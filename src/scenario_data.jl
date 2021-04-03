const vaccination_date_settings =[Date(2021,1,1),Date(2021,9,1)]
const extended_vaccination_date_settings = collect(Date(2020,9,1):Month(1):Date(2021,10,1))
const extended_vaccination_rate_vector = collect(range(0.005,stop = 0.05, step = 0.005))#
const dynamic_threshold_fractions= collect(range(0.5,stop = 3.0, step = 0.25))
const R_0_list = collect(range(1.5,2.5,length = 10))
const small_dynamic_threshold = dynamic_threshold_fractions[1:2:end-1]
const small_ψ_0 =  [0.005,0.01,0.015,0.025,0.05] 
const primary_threshold = 2.0
const efficacy_vectors_elderly = collect(0.4:0.1:0.9) #change
const efficacy_vectors_else = collect(0.4:0.1:0.9)
const kappa_vac_range = collect(range(0.5e2, stop =8e2, step = 2e2))
const c_vac_range = collect(range(2.0e-5, stop =45e-5, step = 9e-5))
const vac_ratio_range = collect(LinRange(1.0,1.27,5))
function get_vaccination_strategies(mean_particle)
    vaccination_strategies_labels = [
        "No vaccination",
        "Oldest first",
        "Youngest first",
        "Uniform",
        "Contact-based"
    ]
    vaccination_strategies = [
        fill(0.0,n_g),
        vcat(fill(0.0,n_g-4),fill(1/4,4)),
        vcat(fill(1/4,4),fill(0.0,n_g - 4)),
        vcat(fill(1/n_g,16)),
        leading_eigevector_of_NGM(mean_particle,contact_matrices.contact_matrices_sum,population_demographics),
        # create_susceptibility_vector(mean_particle.r,mean_particle).*leading_eigevector_of_NGM(mean_particle,contact_matrices.contact_matrices_sum,population_demographics)
    ]
    return (vaccination_strategies_labels,vaccination_strategies )
end
function get_scenario_ranges(scenario)
    @unpack model_time, mean_particle = scenario
    vaccination_strategies_labels,vaccination_strategies = get_vaccination_strategies(mean_particle)
    shutdown_levels = get_capacity_settings(small_dynamic_threshold,mean_particle)
    

    small_parameter_plane = (
        ranges = (
            vaccination_strategies,
            get_vaccination_rates(small_ψ_0),
            vaccination_date_settings,
            get_capacity_settings(small_dynamic_threshold,mean_particle)
        ),
        names = (
            vaccination_strategies_labels,
            small_ψ_0,
            vaccination_date_settings,
            small_dynamic_threshold 
        )
    )

    big_ol_parameter_plane = (
        ranges = (
            vaccination_strategies,
            get_vaccination_rates(extended_vaccination_rate_vector),
            vaccination_date_settings,
            get_capacity_settings(dynamic_threshold_fractions,mean_particle)
        ),
        names = (
            vaccination_strategies_labels,
            extended_vaccination_rate_vector,
            vaccination_date_settings,
            dynamic_threshold_fractions 
        )
    )

    vaccination_date_analysis = (
        ranges = (
            vaccination_strategies,
            get_vaccination_rates(small_ψ_0),
            extended_vaccination_date_settings,
            get_capacity_settings([primary_threshold,5.0],mean_particle)
        ),
        names = (
            vaccination_strategies_labels,
            small_ψ_0,
            extended_vaccination_date_settings,
            [primary_threshold,5.0] 
        ) 
    )
    R_0_analysis = (
        ranges = (
            vaccination_strategies,
            get_vaccination_rates([small_ψ_0[2]]),
            vaccination_date_settings,
            get_capacity_settings([primary_threshold],mean_particle),
            map(R_0 -> get_r(R_0, mean_particle.γ_a,mean_particle.γ_s, mean_particle.σ_1,mean_particle.σ_0, mean_particle.η),R_0_list)
        ),
        names = (
            vaccination_strategies_labels,
            [small_ψ_0[2]],
            vaccination_date_settings,
            [primary_threshold],
            R_0_list
        ) 
    )
    efficacy_analysis = (
        ranges = (
            vaccination_strategies,
            get_vaccination_rates([small_ψ_0[4]]),
            vaccination_date_settings,
            get_capacity_settings([primary_threshold],mean_particle),
            efficacy_vectors_else,
            efficacy_vectors_elderly
        ),
        names = (
            vaccination_strategies_labels,
            [small_ψ_0[4]],
            vaccination_date_settings,
            [primary_threshold],
            efficacy_vectors_else,
            efficacy_vectors_elderly
        ) 
    )
    x_vac_0 = 0.67
    vaccine_refusal_analysis = (
        ranges = (
            [x_vac_0],
            vaccination_strategies,
            get_vaccination_rates([small_ψ_0[2]]),
            vaccination_date_settings,
            get_capacity_settings([primary_threshold],mean_particle),
            c_vac_range,
            kappa_vac_range
        ),
        names = (
            [x_vac_0],
            vaccination_strategies_labels,
            [small_ψ_0[2]],
            vaccination_date_settings,
            [primary_threshold],
            c_vac_range,
            kappa_vac_range
        ) 
    )

    immunization_ratio_plane = (
        ranges = (
            vaccination_strategies,
            get_vaccination_rates([small_ψ_0[4]]),
            vaccination_date_settings,
            get_capacity_settings([primary_threshold],mean_particle),
            vac_ratio_range,
        ),
        names = (
            vaccination_strategies_labels,
            [small_ψ_0[4]],
            vaccination_date_settings,
            [primary_threshold],
            vac_ratio_range
        ) 
    )



    filenames = ("small_parameter_plane.dat","big_parameter_plane.dat","vaccination_date_parameter_plane.dat","R_0_parameter_plane.dat","efficacy_plane.dat","vaccine_refusal.dat","v_T_plane.dat")
    parameter_planes = (small_parameter_plane,big_ol_parameter_plane,vaccination_date_analysis,R_0_analysis,efficacy_analysis,vaccine_refusal_analysis,immunization_ratio_plane)
    vars_list = (
        (:vaccination_distribution, :ψ_0, :vaccination_begin_date,:dynamic_threshold),
        (:vaccination_distribution, :ψ_0, :vaccination_begin_date,:dynamic_threshold),
        (:vaccination_distribution, :ψ_0, :vaccination_begin_date,:dynamic_threshold),
        (:vaccination_distribution, :ψ_0, :vaccination_begin_date,:dynamic_threshold,:r),
        (:vaccination_distribution, :ψ_0, :vaccination_begin_date,:dynamic_threshold, :vaccination_efficacy_else,:vaccination_efficacy_elderly),
        (:x_vac_0,:vaccination_distribution, :ψ_0, :vaccination_begin_date,:dynamic_threshold, :c_vac,:κ_vac),
        (:vaccination_distribution, :ψ_0, :vaccination_begin_date,:dynamic_threshold, :vac_protection_ratio),
    )
    types = (
        TimeseriesData,
        HeatmapData,
        HeatmapData,
        HeatmapData,
        HeatmapData,
        TimeseriesData,
        TimeseriesData,
    )

    return get_capacity_settings([primary_threshold],mean_particle),(vars_list,parameter_planes,filenames,types)
end
function get_opt_params()
    opt_params_baseline =  OrderedDict(
        :σ_0 => (0.0,2.0),
        :σ_1 => (0.0,2.0),
        :γ_s => (0.0,2.0),
        :γ_a => (0.0,2.0),
        :η => (0.0,0.5),
        :κ_1 => (0.0,4e4),
        :κ_2 => (0.0,4e4),
        :t_switch => (120,220),
        :ϕ => (-110,0.0),
        :R_0 => (1.0,2.3),
        :c_1 => (1e-8,1e-3),
        :c_2 => (1e-8,1e-2),
        :ε_P_1 => (0.0,0.5),
        :I_multiplier  => (0.0,12.0),
        :P_multiplier  => (0.0,12.0),
        :E_multiplier  => (0.0,4.0),
        :s => (-0.5,0.0),
        :a_1_1 => (0.01,1.0),
        :a_1_2 => (0.01,1.0),
        :a_1_3 => (0.2,1.0),
        :a_2_1 => (0.01,1.0),
        :a_2_2 => (0.01,1.0),
        :a_2_3 => (0.2,1.0),
        :r_1 => (0.5,2.0),
        :r_2 => (0.5,2.0),
        :r_3 => (0.5,2.0),
    )

    opt_params_R_0_2_5 =  delete!(deepcopy(opt_params_baseline),:R_0)

    opt_params_const_distancing =  OrderedDict(
        :σ_0 => (0.0,2.0),
        :σ_1 => (0.0,2.0),
        :γ_s => (0.0,2.0),
        :γ_a => (0.0,2.0),
        :η => (0.0,0.5),
        :t_switch => (120,220),
        :ϕ => (-110,0.0),
        :R_0 => (1.0,2.3),
        :ε_P_1 => (0.0,0.5),
        :I_multiplier  => (0.0,12.0),
        :P_multiplier  => (0.0,12.0),
        :E_multiplier  => (0.0,4.0),
        :s => (-0.5,0.0),
        :a_1_1 => (0.01,1.0),
        :a_1_2 => (0.01,1.0),
        :a_1_3 => (0.2,1.0),
        :a_2_1 => (0.01,1.0),
        :a_2_2 => (0.01,1.0),
        :a_2_3 => (0.2,1.0),
        :r_1 => (0.5,2.0),
        :r_2 => (0.5,2.0),
        :r_3 => (0.5,2.0),
    )
    return (opt_params_baseline,opt_params_R_0_2_5,opt_params_const_distancing)
end