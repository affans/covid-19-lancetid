module COVID_model

export n_g,get_default_model_parameters,run_model,ModelParameters,bayesian_model_selection,age_bins, create_simulation_scenarios,do_model_plots,create_model_output

using OrdinaryDiffEq
using LinearAlgebra
using ForwardDiff
using DiffEqCallbacks
using StaticArrays
using AxisKeys
import NLopt
using Statistics
using StatsBase
using Distributions
using KissABC
using DataInterpolations
using LabelledArrays
using DelimitedFiles
using Serialization
using DataInterpolations
using Dates
using RData
using UnPack
using DataFrames
using ThreadsX
using StatsPlots
using LaTeXStrings
using UnPack
using DataStructures
using NamedTupleTools

const PACKAGE_FOLDER = dirname(dirname(pathof(COVID_model)))


include("utils.jl")
include("data.jl")
const n_g = 16
const age_bins = vcat([(i,i+5) for i in 0.0:5:74.0],[(75.0,95.0)])
const midpoints_age_bins = midpoints(age_bins)
const population = 14.57e6
const case_age_mortality = SVector{n_g}(rebin_data(get_canada_case_fatality()..., age_bins))
const contact_matrices = load_contact_matrices()
const population_demographics = SVector{n_g}(rebin_data(get_canada_demographic_distribution()..., age_bins) .* population)
const ontario_data = get_ontario_data()
const T = make_T(population_demographics,contact_matrices.contact_matrices_sum)
const school_shutdown_date = Date(2020, 3, 14)
const workplace_shutdown_date =  Date(2020, 3, 17) 
const stage_2_begin_date = Date(2020, 6, 12)
const school_reopening_date = Date(2020, 9, 8) 

include("model.jl")
include("parameter_fitting.jl")
include("plots.jl")
include("scenario_data.jl")

"""
This type acts as a container for parameter fits.
"""
struct ManuscriptScenario{T1,T2}
    model_time::Float64
    mean_particle::T1
    particles::T2
end

"""
        create_simulation_scenarios()

This function fits the model to the data present in /data/, for the specific scenarios used in the paper and specified by `get_opt_params`. The model fits are placed into the /parameter_fits/ directory.

See `get_opt_params()` for the ranges specified for the priors, and `bayesian_model_selection` for the parameter fitting code itself.
"""
function create_simulation_scenarios() 
    model_time = 1825.0

    init_params = apply_workplace_closure_fit!(get_default_model_parameters())
    display(init_params)
    const_distancing_parameters = merge(deepcopy(init_params),(x_0 = 0.99, κ_1 = 0.0, κ_2 = 0.0, c_1 = 0.0, c_2 = 0.0, p_ul = 0.0))
    baseline_opt_params,opt_params_R_0_2_5,opt_const_distancing = get_opt_params()
    fit_scenarios_parameters = (
        baseline = (baseline_opt_params,init_params, 7.0),
        baseline_tight_thresholds = (baseline_opt_params,init_params, 7.0),
        const_distancing = (opt_const_distancing,const_distancing_parameters,252.0),
    )
    for (name,(opt_params,init_params,threshold)) in pairs(fit_scenarios_parameters)
        display(name)
        particles,mean_particle = bayesian_model_selection(init_params,400,threshold,opt_params)
        scenario = ManuscriptScenario(model_time,mean_particle ,sample(particles,200;replace = false))
        serialize("parameter_fits/$name",scenario)
    end
end

function get_vaccination_rates(vaccinations_per_week) 
    return (vaccinations_per_week ./ 7.0).* population
end

function get_capacity_settings(capacity_fractions_of_first_wave,mean_particle)
    model_alloc = allocate_realization(mean_particle,(0.0,200.0),fill(0.0,n_g))
    prob,ascertainment_vector_1,ascertainment_vector_2 = model_alloc
    model_sol,_ = run_model(model_alloc)
    max_active_infections = maximum([sum((u_t.I_a .+ u_t.I_s) .* ascertainment(t,ascertainment_vector_1,ascertainment_vector_2,prob.p.p_tuple.a_m,prob.p.p_tuple.t_switch)) for (t,u_t) in enumerate(model_sol)])
    return capacity_fractions_of_first_wave .* max_active_infections
end

using BenchmarkTools

"""
        create_model_output()

This function creates all the model data based on the parameter fits in /parameter_fits/, and saves it into the directory /output/
# Optional keyword arguments
- `num_particles` allows one to specify a number of particles less than the number specified in the parameter_fit file. Used for debugging or test, when you don't want to wait for the full model to simulate.

- `file_list` is by default, the list of filenames that are used for the simulations in the paper. However, you can pass any vector of filenames present in /parameter_fits/
"""
function create_model_output(;num_particles = -1,file_list = [
                                        "baseline",
                                        "baseline_tight_thresholds",
                                        "const_distancing",
                                    ]
)
    if !isdir(joinpath(PACKAGE_FOLDER,"output"))
        mkdir(joinpath(PACKAGE_FOLDER,"output"))
    end
    map(file_list) do fname
        display(fname)
        dirname = joinpath(PACKAGE_FOLDER,"output",fname)
        scenario = deserialize(joinpath(PACKAGE_FOLDER,"parameter_fits/$fname"))
        thresholds,(vars_list,parameter_planes,filenames,types) = get_scenario_ranges(scenario)
        if !isdir(dirname)
            mkdir(dirname)
        end
        if num_particles > 1
            particles_used = scenario.particles[1:num_particles] #for debugging/tests
        else
            particles_used = scenario.particles
        end
        allocs = COVID_model.allocate_model(particles_used,(0.0,float(length(ontario_data.infection_data_incident[:,1]))),fill(0.0,n_g),COVID_model.ModelFitData)
        fitting_output = COVID_model.run_model!(allocs);
        serialize(joinpath(PACKAGE_FOLDER,dirname,"fitting.dat"),fitting_output)
        

        particles_with_shutdown = map(p-> merge(p,(dynamic_threshold =only(thresholds),)),particles_used)
        allocs = COVID_model.allocate_model(particles_with_shutdown,(0.0,scenario.model_time),fill(0.0,n_g),COVID_model.ModelFitData)
        fitting_output = COVID_model.run_model!(allocs);
        serialize(joinpath(PACKAGE_FOLDER,dirname,"fitting_with_threshold.dat"),fitting_output)

        
        gen_solution_arrays(dirname,scenario.model_time,vars_list,filenames,parameter_planes, particles_used,types)
    end
    return true
end
"""
        do_model_plots()

Creates all the plots in the paper from the output located in /output/. If you generate a subset of the scenarios, then you will need to modify the `file_list` variable in this function accordingly.

However, it works fine with a subset of the particles. 

This function uses a lot of RAM, since we don't precompute any of the things used for plotting other than the raw model data, and it's all loaded into memory at once, because mmapping arrays with custom types is hard and I have a lot of ram. 
"""
function do_model_plots()
    file_list = readdir(joinpath(PACKAGE_FOLDER,"parameter_fits"))
    file_list = [
        "baseline",
        "baseline_tight_thresholds",
        "const_distancing",
        # "higher_kappa",
        # "lower_kappa",
        # "higher_ascertainment",
        # "lower_ascertainment",       
        # "high_R_0",
    ]
  
    for fname in file_list
        println(fname)
        plotdir = joinpath(PACKAGE_FOLDER,"plots",fname)
        plotting_fit_dir = joinpath(plotdir,"fitting_output")
        plotting_fit_with_threshold_dir = joinpath(plotdir,"fitting_output_threshold")
        dirs = [plotdir,plotting_fit_dir,plotting_fit_with_threshold_dir]
        for dirname in dirs
            if !isdir(dirname)
                mkdir(dirname)
            end
        end

        dates_to_use = [1,2]
        scenario = deserialize(joinpath(PACKAGE_FOLDER,"parameter_fits/$fname"))
        fitting_output = deserialize(joinpath(PACKAGE_FOLDER,"output/$fname/fitting.dat"))
        fitting_output_with_threshold = deserialize(joinpath(PACKAGE_FOLDER,"output/$fname/fitting_with_threshold.dat"))
        big_parameter_plane = deserialize(joinpath(PACKAGE_FOLDER,"output/$fname/big_parameter_plane.dat"))
        display(size(big_parameter_plane[1].recovered))

        bivariate_heatmaps(plotdir, big_parameter_plane,ontario_data.start_date,dates_to_use)
  
        smaller_plane = deserialize(joinpath(PACKAGE_FOLDER,"output/$fname/small_parameter_plane.dat"))
       
        plot_vaccination_timeseries(plotdir,smaller_plane,ontario_data.start_date,ontario_data)
        plot_distribution_plane(plotdir,smaller_plane,dates_to_use)
        
        vaccination_date_parameter_plane = deserialize(joinpath(PACKAGE_FOLDER,"output/$fname/vaccination_date_parameter_plane.dat"))
        R_0_parameter_plane = deserialize(joinpath(PACKAGE_FOLDER,"output/$fname/R_0_parameter_plane.dat"))
        efficacy_plane_output = deserialize(joinpath(PACKAGE_FOLDER,"output/$fname/efficacy_plane.dat"))
        vac_refusal_data = deserialize(joinpath(PACKAGE_FOLDER,"output/$fname/vaccine_refusal.dat"))

        plot_model(plotting_fit_dir,fitting_output,ontario_data,scenario.particles)
        violin_plots(plotting_fit_dir,scenario.particles,get_vaccination_strategies(scenario.mean_particle)[2][5])

        plot_model(plotting_fit_with_threshold_dir,fitting_output_with_threshold,ontario_data,scenario.particles)

        plot_x_vac_timeseries(plotdir,vac_refusal_data,ontario_data.start_date,ontario_data)
        bivariate_heatmaps_x_vac(plotdir,vac_refusal_data,ontario_data.start_date,dates_to_use)
        bivariate_heatmaps_vac_efficacy(plotdir,efficacy_plane_output,ontario_data.start_date,dates_to_use)
        vaccination_analysis_plots(plotdir,big_parameter_plane,dates_to_use)
        R_0_plots(plotdir,R_0_parameter_plane,dates_to_use)
        
        vac_date_plots(plotdir,vaccination_date_parameter_plane,ontario_data.start_date)


        v_T_data = deserialize(joinpath(PACKAGE_FOLDER,"output/$fname/v_T_plane.dat"))
        v_T_plots(plotdir,v_T_data,dates_to_use,scenario.mean_particle)
    end
end



end