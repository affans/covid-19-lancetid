export fit_workplace_closure

"""
        apply_workplace_closure_fit!(p)

This function modifies `p`, where `p` is a named tuple of the parameters, such that work_closure_rate, work_opening_rate, max_workplace_closure, and ε_W, match the google mobility data. 

See `fit_workplace_closure` for the actual model fitting code.
"""
function apply_workplace_closure_fit!(p)
    work_close_interval, _ = COVID_model.get_close_intervals(p)
    opening_rate,closing_rate,ε_W = fit_workplace_closure(p,ontario_data.workplace_closure)
    max_closure = maximum([work_closure(float(i),work_close_interval,opening_rate,closing_rate,ε_W,0.0,1.0) for i=0.0:0.1:200.0])
    display((opening_rate,closing_rate,max_closure))
    return merge(p,(work_closure_rate = closing_rate,work_opening_rate = opening_rate,ε_W = ε_W, max_workplace_closure = max_closure ))
end

"""
    fit_workplace_closure(params,workplace_data::Vector{<:Real})

Using NLOpt we fit a sigmoid function to `workplace_data`.
"""
function fit_workplace_closure(params,workplace_data::Vector{<:Real})            
    work_close_interval, _ = get_close_intervals(params)
    x_0 = [
        params.work_opening_rate,
        params.work_closure_rate,
        params.ε_W
    ]

    lb = [0.001,0.001,0.001]
    ub = [1.0,1.0,maximum(workplace_data)]
    function obj(x)
        tot_loss = eltype(x)(0.0)
        for t = 2:length(workplace_data) 
            tot_loss += (work_closure(float(t),work_close_interval,x[1],x[2],x[3],0.0,1.0) - workplace_data[t])^2
        end
        return tot_loss
    end
    function obj_w_grad(x,grad)
        if length(grad) > 0
            ForwardDiff.gradient!(grad,obj,x)
        end
        return obj(x)
    end
    println(obj_w_grad(x_0,[0.0,0.0,0.0]))
    local_optimizer = NLopt.Opt(:LD_MMA,3) 
    local_optimizer.min_objective = obj_w_grad
    local_optimizer.lower_bounds = lb
    local_optimizer.upper_bounds = ub
    local_optimizer.maxtime = 300
    local_optimizer.ftol_rel = 1e-10
    local_optimizer.ftol_abs = 1e-10
    (optf,best,ret) = NLopt.optimize(local_optimizer,x_0)
    println(ret,best)
    return best
end


"""
    unpack(init_params,l,variables)

Automate the process of reassigning the model variables for arbitrary namedtuples of variables, used in `bayesian_model_selection` since the `KissABC` accepts only an AbstractVector of parameters.
"""
function unpack(init_params,l,variables)
    new_p = merge(init_params, NamedTuple{variables}(ntuple(i -> l[i],length(variables))))
    return new_p
end


"""
    bayesian_model_selection(init_params, nparticles::Int64, threshold,opt_params_dict; test = false)

#Arguments 

- init_params: a set of reasonable initial parameters, used only to debug the loss function more easily
- nparticles: number of particles to use for the ABC computation
- threshold: threshold to pass to the ABCDE algorithm
- opt_params: a dict mapping the parameter names to bounds on those parameters, see `get_opt_params` for examples

This function defines the loss used in the fitting as a closure. It assumes that all priors are uniform distributions bounded by the ranges given in `opt_params_dict`
"""
function bayesian_model_selection(init_params, nparticles::Int64, threshold,opt_params_dict; test = false)
    variables =  ntuple(i -> collect(keys(opt_params_dict))[i],length(keys(opt_params_dict)))

    priors = vcat([Uniform(opt_params_dict[k][1],opt_params_dict[k][2]) for k in keys(opt_params_dict)])
    
    @unpack start_date, infection_data_cumulative,infection_data_incident, distancing_data = ontario_data

    data_length = length(infection_data_cumulative[:,1])
    model_time = (0.0,(Float64)(data_length))
    june_30 = Date(2020,6,30)
    
    prevalence_at_june_30 = 0.011
    prevalence_check_days = Dates.value(june_30 - start_date)
    next_gen_matrices = create_next_gen_matrix(init_params)
    threaded_ngm_list = [deepcopy(next_gen_matrices) for i =1:Threads.nthreads()]

    #loss function, takes a solution to the model
    function loss(sol)
        @views begin
            tot_loss = 0.0
            for i ∈ 2:data_length                    
                tot_loss += (1/population)*(sum(sol[i].C .- sol[i-1].C) - sum(infection_data_cumulative[i,:] .- infection_data_cumulative[i - 1,:]))^2 
                tot_loss += 2*(sol[i].x - distancing_data[i])^2 

            end
            for j = 1:n_g
                tot_loss += 0.01*(1/population_demographics[j])*(sol[data_length].C[j]  - infection_data_cumulative[data_length,j])^2
            end
            return tot_loss + (sol[data_length].x- distancing_data[data_length])^2 + 0.1*(sum(sol[prevalence_check_days].R)/population - prevalence_at_june_30)^2
        end
    end

    #cost function defined for the ABCDE algorithm
    function simulate(opt)
        new_p = unpack(init_params,opt,variables)
        ngm = threaded_ngm_list[Threads.threadid()]
        set_sigma_gamma!(ngm,new_p)
        alloc = allocate_realization(new_p,model_time,fill(0.0,n_g); next_gen_matrices = ngm)
        try
            sol,_ = run_model(alloc; ascertainment_flag = true)
            l = loss(sol)
            return l 
        catch e
            display("inf")
            return Inf64
        end
    end

    prior_dist = Factored(priors...)
    output = []
    res,_ = KissABC.ABCDE(prior_dist,simulate,threshold; nparticles = nparticles,generations = 0.5e3,earlystop = false,verbose=true,parallel = true)

    #this bit reformats the output of the ABCDE process (which is in MonteCarloMeasurements particles) to an easier form to work with, and precomputes the r value, using the next gen matrix, and adds it to each parameter tuple.
    for particle in zip(res...)
        nt = unpack(deepcopy(init_params),particle,variables)
        r = get_r(nt.R_0,init_params.γ_a,nt.γ_s,nt.σ_1,nt.σ_0,init_params.η)
        push!(output,merge(nt,(r = r,)))
    end
    set_r(p) = merge(p,(;r = get_r(p.R_0,p.γ_a,p.γ_s,p.σ_1,p.σ_0,p.η),))

    mean_particle = mean.(res)
    mean_particle_params = unpack(deepcopy(init_params),mean_particle,variables)
    
    return map(p->set_r(p),output),set_r(mean_particle_params)
end