
function midpoints(l)
    return [e[1] + (e[2]-e[1])/2 for e in l]
end

function setfields!(s,nt)
    for v in keys(nt)
        setfield!(s,v,nt[v])
    end
end
function rebin_data(in_bins, data, out_bins)
    data_size_adjusted = [d * (b[2] - b[1]) for (d, b) in zip(data, in_bins)]
    
    cumulative_data = [sum(data_size_adjusted[1:i]) for i in 1:length(data)]
    bin_centres = midpoints(in_bins)
    approx_func = CubicSpline(cumulative_data, bin_centres)

    binned_data = [(approx_func(e[2]) - approx_func(e[1])) / (e[2] - e[1]) for e in out_bins]
    return clamp.(binned_data, 0.0, Inf)
end
function interpolate_to_num_compartments(vec, lower_bound, upper_bound)
    interpolate_points = vcat(midpoints_age_bins[begin],[midpoints_age_bins[i] for i = 3:2:n_g -3],midpoints_age_bins[end])
    interp_func = QuadraticInterpolation(vec,interpolate_points)
    return max.([interp_func(e) for e in midpoints_age_bins],0.001)
end
function update_scenario(filename)
    old_scenario = deserialize(filename)
    set_p(p,v,x) = merge(p,NamedTuple{tuple(v)}((x,),))
    default_params = get_default_model_parameters()
    scenario_keys = keys(old_scenario.mean_particle)
    symbols_to_add = Symbol[]
    for v in keys(default_params)
        if !(v in scenario_keys)
            push!(symbols_to_add,v)
        end
    end
    values_to_add = [default_params[v] for v in symbols_to_add]
    nt_to_add = namedtuple(symbols_to_add,values_to_add)
    new_mean_particle = merge(old_scenario.mean_particle,nt_to_add)
    new_particles = map(p -> merge(p,nt_to_add),old_scenario.particles)
    new_scenario = ManuscriptScenario(old_scenario.model_time, new_mean_particle, new_particles)
end

function update_scenario_attr(filename, attributes,values)
    old_scenario = deserialize(filename)
    set_p(p) = merge(p,NamedTuple{attributes}(values,))
    new_mean_particle =set_p(old_scenario.mean_particle)
    new_particles = map(p -> set_p(p),old_scenario.particles)
    new_scenario = ManuscriptScenario(old_scenario.model_time, new_mean_particle, new_particles) 
end

using LinearAlgebra
struct NextGenMatrices{A<:AbstractMatrix,B<:AbstractMatrix}
    ??::A
    next_gen_matrix::B
end

function make_T(N,C)
    T = zeros(Float64,n_g*4,n_g*4)
    T[1:n_g,n_g+1:n_g*2] .= C
    for i in 1:n_g
        for j in n_g+1:n_g*2
            T[i,j] *= (N[i]/N[j-n_g])
            T[i,j + n_g] = T[i,j]
            T[i,j + n_g*2] = T[i,j]
        end
    end
    return T
end

function set_Sigma!(??,??_a,??_s,??_1,??_0,??)
    ??[1:n_g,1:n_g] .= Diagonal(fill(-1*??_0,n_g))
    ??[n_g+1:n_g*2,n_g+1:n_g*2] .= Diagonal(fill(-1*??_1,n_g))
    ??[n_g*2  + 1 : n_g*3,n_g*2  + 1 : n_g*3] .= Diagonal(fill(-1*??_a,n_g))
    ??[3*n_g  + 1 : n_g*4, 3*n_g  + 1 : n_g*4] .= Diagonal(fill(-1*??_s,n_g))

    ??[n_g+1: 2*n_g, 1:n_g].= Diagonal(fill( ??_0,n_g))
    ??[n_g*2+1: 3*n_g, n_g + 1:n_g*2].= Diagonal(fill( ??*??_1,n_g))
    ??[3*n_g+1:end, n_g + 1:n_g*2].= Diagonal(fill((1 - ??)*??_1,n_g))
end

function get_r(R_0,??_a,??_s,??_1,??_0,??)
    ?? = zeros(Float64,n_g*4, n_g*4)
    set_Sigma!(??,??_a,??_s,??_1,??_0,??)
    ngm = -1 .* T * inv(??)
    evals,evecs = eigen(ngm)
    vals = real.(evals)
    leading_eigenvector = real.(evecs[1:n_g,end])
    return R_0/maximum(vals)
end


@views function create_next_gen_matrix(params_tuple)::NextGenMatrices
    ??_a = params_tuple.??_a
    ??_s = params_tuple.??_s
    ??_1 = params_tuple.??_1

    ??_0 = params_tuple.??_0
    ?? = params_tuple.??

    ?? = zeros(Float64,n_g*4, n_g*4)
    set_Sigma!(??,??_a,??_s,??_1,??_0,??)
    ngm = -1 .* T * inv(??)

    return NextGenMatrices(??, -1 .* T * inv(??))
end
function leading_eigevector_of_NGM(params_tuple,C,demographics)


    ??_a = params_tuple.??_a
    ??_s = params_tuple.??_s
    ??_0 = params_tuple.??_0
    ??_1 = params_tuple.??_1
    ?? = params_tuple.??

    ?? = zeros(Float64,n_g*4, n_g*4)
    set_Sigma!(??,??_a,??_s,??_1,??_0,??)
    T = make_T(demographics,C)
    ngm = -1 .* T * inv(??)
    vals,vecs = eigen(ngm)
    leading_eigenvector = real.(vecs[1:n_g,end])
    # @show leading_eigenvector
    return leading_eigenvector/sum(leading_eigenvector)
end
@views function set_sigma_gamma!(next_gen_matrices::NextGenMatrices, params_tuple)
    ??_a = params_tuple.??_a
    ??_s = params_tuple.??_s
    ??_1 = params_tuple.??_1

    ??_0 = params_tuple.??_0
    ?? = params_tuple.??
    set_Sigma!(next_gen_matrices.??,??_a,??_s,??_1,??_0,??)
    next_gen_matrices.next_gen_matrix .=  -1 .* T * inv(next_gen_matrices.??)
end

function get_r(R_0::S,next_gen_matrices::NextGenMatrices)::S where S<:Real
    vals = real.(eigen(next_gen_matrices.next_gen_matrix).values)
    return R_0/maximum(vals)
end