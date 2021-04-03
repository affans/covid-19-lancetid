abstract type ModelOutput end


struct ModelAllocation{V,T2,MOType<:ModelOutput}
    alloc_tuples::V
    constant_params::T2
    model_output_data::MOType
end
struct ModelFitData{V,V2} <: ModelOutput
    incident_infections::V
    distancing::V
    work_shutdown::V
    school_shutdown::V
    ascertained_cumulative_infections::V2
    recovered::V2
end
struct TimeseriesData{V,V2} <: ModelOutput
    incident_infections::V
    distancing::V
    vaccine_refusal::V
    work_shutdown::V
    susceptible::V
    school_shutdown::V
    recovered::V
    vaccinated::V
    cumulative_infections_from_t_vac::V2
end
struct HeatmapData{V,V2} <: ModelOutput
    recovered::V
    cumulative_infections_at_end::V2
end
function ModelFitData(nparticles::Int64,len::Int64)
    return ModelFitData(zeros(nparticles,len),zeros(nparticles,len), zeros(nparticles,len),zeros(nparticles,len), zeros(nparticles,len,n_g), zeros(nparticles,len,n_g))
end

function TimeseriesData(nparticles::Int64,len::Int64)
    return TimeseriesData(zeros(nparticles,len),zeros(nparticles,len), zeros(nparticles,len),zeros(nparticles,len),zeros(nparticles,len),zeros(nparticles,len),zeros(nparticles,len),zeros(nparticles,len),zeros(nparticles))
end
function HeatmapData(nparticles::Int64,len::Int64)
    return HeatmapData(zeros(nparticles,len),zeros(nparticles))
end
function get_u_0(element_type::Type, initial_I_data::Array{<:Real,1}, p_tuple)::(T2 where T2 <:AbstractVector{<:element_type})
    u_0 = @LArray fill((element_type)(0.0),n_g*9 + 3) (
        S = 1:n_g,
        S_v_1 = n_g+1:n_g*2,
        E = n_g*2 + 1:n_g*3,
        P = n_g*3 + 1:n_g*4,
        I_a = n_g*4 + 1:n_g*5,
        I_s = n_g*5 + 1:n_g*6,
        R = n_g*6 + 1:n_g*7,
        C = n_g*7 + 1:n_g*8,
        V = n_g*8 + 1:n_g*9,
        x = n_g*9 + 1,
        D = n_g*9 + 2,
        x_vac = n_g*9 + 3,
    )
    u_0.S .= population_demographics .-  (p_tuple.E_multiplier + p_tuple.P_multiplier + p_tuple.I_multiplier) .* initial_I_data
    u_0.E .= p_tuple.E_multiplier .* initial_I_data
    u_0.P .= p_tuple.P_multiplier .*initial_I_data
    u_0.I_a .= p_tuple.I_multiplier .* p_tuple.η .* initial_I_data
    u_0.I_s .= p_tuple.I_multiplier .* (1-p_tuple.η) .* initial_I_data
    u_0.C .= (p_tuple.E_multiplier + p_tuple.P_multiplier + p_tuple.I_multiplier) .* initial_I_data
    u_0.x = p_tuple.x_0
    u_0.x_vac = p_tuple.x_vac_0
    return u_0
end
function smoothstep(x,threshold,k)
    return tanh(k*(x - threshold))*0.5 + 0.5
end
@views function work_closure(t::T2,work_close_interval::Tuple{Int64,Int64}, work_opening_rate,work_closing_rate,ε_W, D, max_workplace_closure) where {T2<:Real}
    accum = 0.0
    accum +=  ε_W*(tanh(work_opening_rate*(t - work_close_interval[1])) - tanh(work_closing_rate*(t - work_close_interval[2])))*0.5
    accum = min(1.0, accum +  ε_W*D*smoothstep(t,243,0.1)) #TODO: don't hardcode this lol, end of fitting window
    return clamp(accum,0.0, ε_W)
end
@views function school_closure(t::T2,school_close_interval::Tuple{Int64,Int64}, k::T2, D) where {T2<:Real}
    accum = 0.0 
    accum += (tanh(0.5*(t- school_close_interval[1]))- tanh(k*(t- school_close_interval[2])))*0.5
    accum = min(1.0, accum + D*smoothstep((t ),243,0.1) ) #TODO: don't hardcode this lol, end of fitting window
    return min(accum, 1.0)
end
function seasonality(t, seasonality_phase, seasonality_modifier)
    return sin((2*pi/365)*(t - seasonality_phase) - pi/2)*seasonality_modifier + 1.0 
end

function ascertainment(t,ascertainment_vector_1,ascertainment_vector_2,steepness,t_switch)
    return  t > t_switch ? ascertainment_vector_2 : ((t_switch- t)/t_switch) .* ascertainment_vector_1 .+ (t/t_switch).*ascertainment_vector_2#(smoothstep(t,steepness,0.05) * ( 1- begin_factor) + begin_factor)#smoothstep(t,t_1,steepness) .* ascertainment_vector_2 .+ (1 .- smoothstep(t,t_1,steepness)) .* ascertainment_vector_1#
end

function cases_from_travel_proportion(t)
    index = Int(round(t)) + 1
    return index <= ontario_data.data_length ? ontario_data.cases_from_travel[index,:]./population_demographics : ontario_data.cases_from_travel[end,:]./population_demographics #assume cases_from_travel for sim is equal to most recent 7-day average
end
function x_boundary_function(x,p_ul)
    return p_ul*(1 - 2*x)#x<0.5 ? p_ul*(1 - 2x)^p_e : -p_ul*(abs(1 - 2x))^p_e #
end
function contacts!(t,matrix_cache,x_t,p,school_closure_t,work_closure_t)
    @. matrix_cache = ( (1 - p.ε_P_1 * x_t)*(contact_matrices.other_contacts + contact_matrices.home_contacts) + (1 - school_closure_t)*contact_matrices.school_contacts + (1 - work_closure_t)*contact_matrices.work_contacts )
end

"""
This is the RHS of the main ODE system. It has a lot of cache-swapping but that seems to improve performance by about 4x.
"""
function system(du,u,p_tuple,t)
    vector_cache_1,vector_cache_2,matrix_cache,r,closed_intervals,p,ascertainment_vector_1,ascertainment_vector_2 = p_tuple
    @inbounds @views begin
        κ = smoothstep(t,p.t_switch,p.a_m)*p.κ_2 + (1 - smoothstep(t,p.t_switch,p.a_m)) * p.κ_1
        c = smoothstep(t,p.t_switch,p.a_m)*p.c_2 + (1 - smoothstep(t,p.t_switch,p.a_m)) * p.c_1
        work_closure_t = work_closure((t-28),closed_intervals.work,p.work_opening_rate,p.work_closure_rate,p.ε_W,u.D,1.0)
        school_closure_t = school_closure((t-28),closed_intervals.school,0.5,u.D)
        contacts!(t,matrix_cache,u.x,p,school_closure_t,work_closure_t)
        @. vector_cache_2 = u.I_a + u.I_s
        vector_cache_1 .= ascertainment(t,ascertainment_vector_1,ascertainment_vector_2,p.a_m,p.t_switch)
        du.x = κ * u.x * (1.0-u.x) *((1/population) * dot(vector_cache_2,vector_cache_1) - c*u.x) + x_boundary_function(u.x,p.p_ul)
        du.x_vac = p.κ_vac * u.x_vac * (1.0 - u.x_vac) *((1/population) * dot(vector_cache_2,vector_cache_1) - p.c_vac)
        du.D = dot(vector_cache_2,vector_cache_1)>p.dynamic_threshold ? 0.05*(1 - u.D) : (-(p.work_closure_rate)* u.D)
        
        
        @. vector_cache_1 = (vector_cache_2 + u.P) / population_demographics


        mul!(vector_cache_2, matrix_cache,vector_cache_1)
        vector_cache_1 .= (r .* seasonality(t,p.ϕ,p.s) .* vector_cache_2 .+ cases_from_travel_proportion(t))
        
        @. du.S =  -1 * vector_cache_1 * u.S
        @. du.S_v_1 = -1 * vector_cache_1 * u.S_v_1
        @. du.E = vector_cache_1 * (u.S .+ u.S_v_1) - p.σ_0 * u.E
        @. du.P = p.σ_0 * u.E - p.σ_1 * u.P
        @. du.I_a = p.η * p.σ_1 * u.P - p.γ_a* u.I_a
        @. du.I_s = (1 - p.η)* p.σ_1 * u.P - p.γ_s * u.I_s
        @. du.R =  p.γ_a * u.I_a + p.γ_s * u.I_s
        @. du.C = p.σ_1 * u.P
       
    end
end

function get_close_intervals(p)
    c_i = (
        work = (Dates.value(workplace_shutdown_date - p.infection_start),Dates.value( stage_2_begin_date - p.infection_start)),
        school = (Dates.value(school_shutdown_date - p.infection_start),Dates.value( school_reopening_date - p.infection_start))
    )
    return c_i
end
function create_ascertainment_vectors(parameter_tuple)
    ascertainment_vector_1 = SVector{n_g}(vcat(fill(parameter_tuple.a_1_1,4),fill(parameter_tuple.a_1_2,8),fill(parameter_tuple.a_1_3,4)))
    ascertainment_vector_2 = SVector{n_g}(vcat(fill(parameter_tuple.a_2_1,4),fill(parameter_tuple.a_2_2,8),fill(parameter_tuple.a_2_3,4)))
    return ascertainment_vector_1,ascertainment_vector_2
end
function create_susceptibility_vector(r,parameter_tuple)
    susceptibility_vector = SVector{n_g}(vcat(fill(parameter_tuple.r_1,4),fill(parameter_tuple.r_2,8),fill(parameter_tuple.r_3,4)))
    return r .* susceptibility_vector
end

"""
Allocate everything we need to create a model solution, without actually solving it. 
"""
function allocate_realization(parameter_tuple,
    t_span,
    vaccine_distribution;
    next_gen_matrices = nothing
)

    ascertainment_vector_1,ascertainment_vector_2 = create_ascertainment_vectors(parameter_tuple)
    if :r in keys(parameter_tuple)
        susceptibility_vector =  create_susceptibility_vector(parameter_tuple.r, parameter_tuple)
    else
       susceptibility_vector = create_susceptibility_vector( get_r(parameter_tuple.R_0,next_gen_matrices) ,parameter_tuple)
    end
    parameter_types = eltype.(filter(x -> typeof(x) <: Real,values(parameter_tuple)))
    element_type = promote_type(parameter_types...,eltype(vaccine_distribution))
    vector_cache_1 = SizedVector{n_g}(zeros(element_type,n_g))
    vector_cache_2 = SizedVector{n_g}(zeros(element_type,n_g))
    matrix_cache = zeros(element_type,(n_g,n_g))
    u_0 = get_u_0(element_type,parameter_tuple.I_0 ./ ascertainment(0,ascertainment_vector_1,ascertainment_vector_2,parameter_tuple.a_m,parameter_tuple.t_switch),parameter_tuple)
    t_vac = float(Dates.value(parameter_tuple.vaccination_begin_date - parameter_tuple.infection_start))
    stop_times = collect(t_vac:1:t_span[2])

    function fadeout!(integrator)
        integrator.u.I_a .= 0.0
        integrator.u.I_s .= 0.0
    end
    cb_fadeout = DiscreteCallback(
        (u,t,integrator)-> sum(u.I_a .+ u.I_s)<1.0,
        fadeout!,
        save_positions=(false,false)
    )
    transfer_rates = SavedValues(element_type,Vector{element_type})

    function saving_func(u,t,int)
        return (u.S .+ u.S_v_1)./(u.S .+ u.S_v_1 .* ((1-int.p.p_tuple.vaccination_efficacy_else)/(1-int.p.p_tuple.vaccination_efficacy_else*int.p.p_tuple.vac_protection_ratio)))
    end

    saving_cb = SavingCallback(
        saving_func,
        transfer_rates;
        saveat = collect(1.0:1:t_span[2])
    ) 
    vaccination_efficacy =SVector{n_g}(vcat(fill(parameter_tuple.vaccination_efficacy_else,n_g - 4), fill(parameter_tuple.vaccination_efficacy_elderly,4)))

    if !all(vaccine_distribution .== 0)
        cb_vac = PresetTimeCallback(stop_times,
            int->vaccinate_callback!(
                int.u.S,int.u.V,int.u.S_v_1,
                parameter_tuple.ψ_0, vaccine_distribution,vaccination_efficacy,int.u.x_vac
            ),
            save_positions=(false,false)
        )
        cb_set = [cb_vac,cb_fadeout,saving_cb]
    else
        cb_set = [cb_fadeout,saving_cb]
    end
    
    prob = ODEProblem(system,u_0,t_span,(c1 = vector_cache_1,
        c2 = vector_cache_2,c3 = matrix_cache,
        r = susceptibility_vector,c_I = get_close_intervals(parameter_tuple),
        p_tuple = parameter_tuple,a_v_1 = ascertainment_vector_1,a_v_2 = ascertainment_vector_2),
        saveat = 1.0, callback = CallbackSet(cb_set...)
    )
    return (prob, ascertainment_vector_1,ascertainment_vector_2,transfer_rates)
end


"""
Do `allocate_realization` for every particle/parameter_set in `parameter_tuple_list` and store them all in a `ModelAllocation` type. 
"""
function allocate_model(parameter_tuple_list,
    t_span,
    vaccine_distribution,model_output_type
    ;
)
    allocation = map(x -> allocate_realization(x,t_span,vaccine_distribution),parameter_tuple_list)
    model_output = model_output_type( length(parameter_tuple_list),Int(trunc(t_span[2]-t_span[1])))
    return ModelAllocation(allocation,parameter_tuple_list[1], model_output)
end

"""
Solve the model given by `alloc` (output from `allocate_realization`), and return the solution. Used in `bayesian_model_selection`.
"""
function run_model(alloc; ascertainment_flag = false)
    prob,ascertainment_vector_1,ascertainment_vector_2,mortality_adjustment = alloc
    odesol = solve(prob,Tsit5(),maxiters = 5e4,verbose = false)
    if any((s.retcode != :Success for s in odesol))
         error("solution error")
    end
    u = odesol.u
    tlist = odesol.t
    if ascertainment_flag
        @views for j=2:length(u)
            odesol.u[j].C .*= ascertainment(j,ascertainment_vector_1,ascertainment_vector_2,prob.p.p_tuple.a_m,prob.p.p_tuple.t_switch)
        end
    end
    return odesol.u,mortality_adjustment.saveval
end


"""
Solve the models given by `allocations` (output from `allocate_model`), and store the solution in `allocations.model_output_data`. Used to create most of the actual model output that we plot.
"""
function run_model!(allocations::ModelAllocation{V,T,MOType};stiff = false) where {V,T,MOType}
    @unpack alloc_tuples,constant_params,model_output_data = allocations
    for (i,alloc) in enumerate(alloc_tuples)
        prob,ascertainment_vector_first,ascertainment_vector_second,mortality_adjustment = alloc
        odesol = solve(prob,Tsit5(),maxiters = 1e4,verbose = false)
        if any((s.retcode != :Success for s in odesol))
            println("TRYING WITH HIGHER TOL")
            odesol = solve(prob,Tsit5(),reltol = 1e-4,verbose = false)
            if any((s.retcode != :Success for s in odesol))
                println("NOPE GIVING UP")
                error("solution too unstable")
            end
        end
        build_model_output!(i,odesol,model_output_data,constant_params,ascertainment_vector_first,ascertainment_vector_second,prob.p.p_tuple.a_m,prob.p.p_tuple.t_switch,mortality_adjustment.saveval)
    end
    return model_output_data
end

function build_model_output!(iteration_num,sol::ODESolution,model_output_struct::ModelFitData,constant_params,ascertainment_vector_1,ascertainment_vector_2,steepness,t2,mortality_adjustment)
    i = iteration_num
    u = sol.u
    tlist = sol.t
    w_c_intervals,s_c_intervals = get_close_intervals(constant_params)
    @views for j=2:length(u)-1
        u_t = u[j]
        t = tlist[j]
        model_output_struct.incident_infections[i,j] = sum((u_t.C -  u[j-1].C) .* ascertainment(j,ascertainment_vector_1,ascertainment_vector_2,steepness,t2))
        model_output_struct.distancing[i,j] = u_t.x
        model_output_struct.recovered[i,j,:] .= u_t.R
        model_output_struct.work_shutdown[i,j] = COVID_model.work_closure(t,w_c_intervals,constant_params.work_opening_rate,constant_params.work_closure_rate,constant_params.ε_W, u_t.D,1.0)
        model_output_struct.school_shutdown[i,j] = COVID_model.school_closure(t, s_c_intervals,0.5,u_t.D)
        model_output_struct.ascertained_cumulative_infections[i,j,:] .= u_t.C .* ascertainment(j,ascertainment_vector_1,ascertainment_vector_2,steepness,t2)
    end
end
function build_model_output!(iteration_num,sol::ODESolution,model_output_struct::TimeseriesData,constant_params,ascertainment_vector_1,ascertainment_vector_2,steepness,t2,mortality_adjustment)
    i = iteration_num
    u = sol.u
    tlist = sol.t

    w_c_intervals,s_c_intervals = get_close_intervals(constant_params)
    t_vac = Dates.value(constant_params.vaccination_begin_date - constant_params.infection_start)
    @views for j=2:length(u)-1
        u_t = u[j]
        t = tlist[j]
        model_output_struct.incident_infections[i,j] =  sum((u_t.C -  u[j-1].C) .* ascertainment(j,ascertainment_vector_1,ascertainment_vector_2,steepness,t2))
        model_output_struct.distancing[i,j] = u_t.x
        model_output_struct.work_shutdown[i,j] = COVID_model.work_closure(t - 28,w_c_intervals,constant_params.work_opening_rate,constant_params.work_closure_rate,constant_params.ε_W, u_t.D,1.0)
        model_output_struct.school_shutdown[i,j] = COVID_model.school_closure(t - 28, s_c_intervals,0.5,u_t.D)
        model_output_struct.recovered[i,j] = sum(u_t.R)
        model_output_struct.vaccinated[i,j] = sum(u_t.V)
        model_output_struct.susceptible[i,j] = sum(u_t.S)
        model_output_struct.vaccine_refusal[i,j] = u_t.x_vac
        if j > t_vac
            model_output_struct.cumulative_infections_from_t_vac[i] += sum((u_t.C -  u[j-1].C) .* case_age_mortality .* mortality_adjustment[j])
        end
    end
  
end
function build_model_output!(iteration_num,sol::ODESolution,model_output_struct::HeatmapData,constant_params,ascertainment_vector_1,ascertainment_vector_2,steepness,t2,mortality_adjustment)
    i = iteration_num
    u = sol.u
    tlist = sol.t
    w_c_intervals,s_c_intervals = get_close_intervals(constant_params)
    @views for j=2:length(u)-1
        u_t = u[j]
        t = tlist[j]
        model_output_struct.recovered[i,j] = sum(u_t.R)
        model_output_struct.cumulative_infections_at_end[i] += sum((u_t.C -  u[j-1].C).* case_age_mortality  .* mortality_adjustment[j])
    end
end
function fadeout!(integrator)
    integrator.u.I_a .= 0.0
    integrator.u.I_s .= 0.0
end

"""
Defines the callback function called by the ODE solver each day after `t_vac`, which modifies the ODE system state corresponding to vaccinating as per the paper description.
"""
function vaccinate_callback!(susceptible, vaccinated_state, vaccinated_not_immunized_1, ψ_0, vaccine_distribution, vaccine_efficacy, x_vac)
    total_susceptible = susceptible .+ vaccinated_not_immunized_1  
    vaccination_rate_modifier = (total_susceptible ./ (population_demographics .- vaccinated_state) )
    function vaccination_step!(susceptible,vaccinated_not_immunized,vaccine_distribution,num_vaccines,vaccine_efficacy,vaccination_rate_modifier,x_vac,population_demographics)
        vaccinated = min.( vaccination_rate_modifier .* num_vaccines .* vaccine_distribution, susceptible .* x_vac)
        #vector of vaccinated people in each vaccination step is the number of vaccines going to susceptible people (vaccination_rate_modifier .* num_vaccines) 
        #divided according to the vac strategy (vaccine_distribution)
        #or, the number of remaining susceptible people who are willing to be vaccinated S_i.* x_vac
        #whichever is less
        susceptible .= susceptible .- vaccinated
        vaccinated_state .+= vaccine_efficacy  .* vaccinated
        vaccinated_not_immunized .+=  (1 .- vaccine_efficacy)  .* vaccinated
        return num_vaccines - sum(vaccinated ./vaccination_rate_modifier)
    end
    remaining_vaccines = vaccination_step!(susceptible,vaccinated_not_immunized_1,vaccine_distribution,ψ_0,vaccine_efficacy,vaccination_rate_modifier,x_vac,population_demographics)
    while remaining_vaccines >= 1.0

        compartments_not_empty = count(susceptible .> ((population_demographics .* (1 - x_vac)) .+ 1e-2)) 
        if compartments_not_empty == 0
            break
        else
            remaining_vaccines = vaccination_step!(susceptible,vaccinated_not_immunized_1,(1/compartments_not_empty),remaining_vaccines,vaccine_efficacy,vaccination_rate_modifier,x_vac,population_demographics)
        end
    end
end

"""
This function calls `run_model!` on `particles_used` for every parameter set in parameter_plane[k].ranges, for every index k. 
"""
function gen_solution_arrays(folder,model_time,vars_list,filenames,parameter_planes,particles_used,types)
    for (parameter_plane, filename,vars,output_data_type) in zip(parameter_planes,filenames,vars_list,types)
            pp_data = ThreadsX.map(collect(Iterators.product(parameter_plane.ranges...))) do p_tuple
                nt = NamedTuple{vars}(p_tuple)
                particles_used_new = map(p->merge(p,nt),particles_used)
                alloc = allocate_model(particles_used_new,(0.0,model_time),nt.vaccination_distribution, output_data_type)
                return run_model!(alloc)
            end
            axes_labels = NamedTuple{vars}(parameter_plane.names)
            serialize("$folder/$filename",KeyedArray(pp_data; axes_labels...))
    end
end