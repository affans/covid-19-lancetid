const color_palette = palette(:seaborn_pastel)
using Printf
pgfplotsx()
default(dpi = 300)
default(framestyle = :box)
function get_stats(data)
    median_ts = dropdims(mapslices(x -> median(x),data,dims = [1]),dims = 1)
    max_ts = dropdims(mapslices(x -> quantile(x,0.95),data,dims = [1]),dims = 1)
    min_ts = dropdims(mapslices(x -> quantile(x,0.05),data,dims = [1]), dims = 1)
    return (median = median_ts,ribbon = (median_ts .- min_ts, max_ts .- median_ts))
end
function plot_contact_matrices()
    p_list = [plot() for p in 1:4]
    contact_matrices_keys = filter(k -> k != :contact_matrices_sum, keys(contact_matrices))
    titles = [
        "Other",
        "School",
        "Work",
        "Home",
    ]
    for (C,p,title) in zip(contact_matrices_keys,p_list,titles)
        heatmap!(p,midpoints_age_bins,midpoints_age_bins,contact_matrices[C];  xlabel = "Age of individual", c=cgrad(palette(:Blues)),title = title)
    end
    # add_subplot_letters!(p_list)
    plot!(p_list[begin]; ylabel = "Age of contact")
    plot!(p_list[end]; ylabel = "Age of contact")
    
    p = plot(p_list..., layout = (1,4), size=(800,200))
    savefig("plots/contact_matrices.pdf")
end
function seasonality_stats(particles,l)
    seasonality_for_each_particle = hcat(map(p->[COVID_model.seasonality(t,p.ϕ,p.s) for t in range(0.0,step = 1.0,stop = l)],particles)...)
    return get_stats(transpose(seasonality_for_each_particle))
end
function ascertainment_stats(particles,l)
    ascertainment_by_age = []
    v_list = zip(create_ascertainment_vectors.(particles),particles)
    # display([ascertainment(t,p.ascertainment_vector_1,p.ascertainment_vector_2) for t in range(0.0,step = 1.0,stop = l)])
    ascertainment_for_each_particle = transpose(hcat(map(v_i->[ascertainment(t,v_i[1][1],v_i[1][2],v_i[2].a_m,v_i[2].t_switch) for t in range(0.0,step = 1.0,stop = l)],v_list)...))

    for i = 1:n_g
        # display(hcat(ascertainment_for_each_particle[i]...))
        push!(ascertainment_by_age,get_stats(map(x->x[i],ascertainment_for_each_particle)))
    end
    return ascertainment_by_age
end
function cases_from_travel(particles,l)
    cases_from_travel_by_age = []
    cases_from_travel = transpose(hcat(map(p->[COVID_model.cases_from_travel_proportion(t) for t in range(0.0,step = 1.0,stop = l)],particles)...))
    for i = 1:n_g
        # display(hcat(ascertainment_for_each_particle[i]...))
        push!(cases_from_travel_by_age,get_stats(map(x->x[i],cases_from_travel)))
    end
    return cases_from_travel_by_age
end

#plot the model fit
function plot_model(folder::String,sol_data,ontario_data,particles) where V<:AbstractVector
    @unpack start_date,infection_data_cumulative, distancing_data,infection_data_incident,infection_data_incident_raw,workplace_closure = ontario_data
    display(start_date)
    infection_ts = get_stats(sol_data.incident_infections)
    distancing_ts = get_stats(sol_data.distancing)
    work_shutdown_ts = get_stats(sol_data.work_shutdown)
    school_shutdown_ts = get_stats(sol_data.school_shutdown)
    data_length = length(infection_data_cumulative[:,end])
    date_x_points = collect(start_date:Day(1):(start_date + Day(2000))) #lol
    x_points_with_labels= collect(ontario_data.start_date:Month(2):(ontario_data.start_date + Day(999)))
    plot_x_ticks_labels = map(s-> "\\textrm{" * s* "}",Dates.format.(x_points_with_labels,"uuu d, yyyy"))

    p_I = plot(date_x_points[1:length(infection_ts.median)],infection_ts.median,ribbon = infection_ts.ribbon, color = color_palette[1], label="Model, 95\\% confidence",title = "",ylabel = "Incident cases")
    p_x = plot(date_x_points[1:length(distancing_ts.median)],distancing_ts.median,ribbon =  distancing_ts.ribbon, color = color_palette[2], label = "x(t), 95\\% confidence", ylabel = "Fraction of contacts")

    plot!(p_x,date_x_points[1:length(work_shutdown_ts.median)], work_shutdown_ts.median, color = color_palette[3],label = "Work shutdown")
    plot!(p_x,date_x_points[1:length(school_shutdown_ts.median)], school_shutdown_ts.median, color = color_palette[5],label = "School shutdown")
    scatter!(p_x,date_x_points[1:length(workplace_closure)],workplace_closure, label = "Contacts in workplace data",  markersize = 1.5,markerstrokewidth = 0.6,  markerstrokealpha = 0.7, color = color_palette[3])
    scatter!(p_x,date_x_points[1:length(distancing_data)],distancing_data,color = color_palette[2], label = "Distancing data",  markersize = 2,markerstrokewidth = 0.6,  markerstrokealpha = 0.7)
    scatter!(p_I,date_x_points[1:length(sum.(eachrow(infection_data_incident)))],sum.(eachrow(infection_data_incident)),label = "Case notification data",  markersize = 2, markerstrokewidth = 0.6, markerstrokealpha = 0.7, color = color_palette[1])

    #scatter!(p_active,Dates.value(active_cases_data_start_date - start_date):(length(active_cases) + Dates.value(active_cases_data_start_date - start_date) - 1), active_cases, label = "active cases data")
    seasonality_ts = seasonality_stats(particles,float(data_length))
    # p_3 = plot(date_x_points[1:data_length+1],
    # seasonality_ts.median,
    #     ribbon = seasonality_ts.ribbon,
    #     label = "Seasonality curve", ylabel = "Susceptibility factor",color = color_palette[1])
    # seasonality_list = [COVID_model.seasonality(t,constant_params.ϕ,constant_params.s) for t in range(0.0,step = 1.0,stop = 365.0)]
    # println("max = $(findmax(seasonality_list))")

    plot_list = (p_I,p_x)
    add_subplot_letters!(plot_list)
    plot(plot_list...,layout = (2,1),size = (400,400),xlabel = "time (date)",legend = :outerright,xticks = (x_points_with_labels, plot_x_ticks_labels))
    savefig("$folder/plot_model.pdf")
    
    labels = midpoints_age_bins
    # display(data_length)#sum.(eachrow(sol_data.ascertained_cumulative_infections[:,data_length,:])))
    # display(repeat(labels,inner = n_g))
    # display(vcat(collect(eachrow(sol_data.ascertained_cumulative_infections[:,data_length,:]))...))
    p = StatsPlots.violin(repeat(labels,outer = length(particles)),vcat(collect(eachrow(sol_data.ascertained_cumulative_infections[:,data_length,:]))...); seriescolor = color_palette[1], label = "Model cumulative infections on  $(particles[1].infection_start + Day(data_length))",size = (400,300),xrotation = 60)
    scatter!(p,labels,infection_data_cumulative[data_length,:],size = (600,300), markersize = 4,markerstrokewidth = 0.2,  markerstrokealpha = 0.7, color = color_palette[2],label = "Totals cases observed up to $(particles[1].infection_start + Day(data_length))",ylabel = "Cumulative infections on $(particles[1].infection_start + Day(data_length))")
    savefig(p, "$folder/cumulative_infections_by_age.pdf")
    

    ascertainment_ts_list = ascertainment_stats(particles,data_length)
    
    ascertainment_plots = [plot() for i=1:n_g]
    for i=1:n_g
        plot!(ascertainment_plots[i],date_x_points[1:data_length+1],ascertainment_ts_list[i].median,ribbon = ascertainment_ts_list[i].ribbon, title = "Ages: $(age_bins[i])",ylim = [0.0,1.0], label = "Case ascertainment")
    end
    plot(ascertainment_plots...,size = (1800, 800),xlabel = "time (days)")
    savefig("$folder/ascertainment_over_time.pdf")

    # cases_ts_list = mortality_stats(particles,data_length)
    
    # mortality_plots = [plot() for i=1:n_g]
    # for i=1:n_g
    #     plot!(mortality_plots[i],date_x_points[1:data_length+1],cases_ts_list[i].median,ribbon = cases_ts_list[i].ribbon, title = "Ages: $(age_bins[i])")
    # end
    # plot(mortality_plots...,layout = (n_g,1),size = (300,2600),xlabel = "time (days)")
    # savefig("$folder/mortality_over_time.pdf")

    cases_from_travel_list = cases_from_travel(particles,data_length)
    
    cases_from_travel_plots = [plot() for i=1:n_g]
    for i=1:n_g
        plot!(cases_from_travel_plots[i],date_x_points[1:data_length+1],cases_from_travel_list[i].median,ribbon = cases_from_travel_list[i].ribbon, title = "Ages: $(age_bins[i])", label = "Observed cases from travel")
    end
    plot(cases_from_travel_plots...,layout = (n_g,1),size = (500,2600),xlabel = "time (days)")
    savefig("$folder/cases_from_travel_plots_over_time.pdf")


    june_30 = Date(2020,6,30)
    prevalence_check_days = Dates.value(june_30 - start_date)
    recovered_per_day = dropdims(mapslices(sum,sol_data.recovered,dims = (3)),dims = 3) / population
    seroprevalence_ts = get_stats(recovered_per_day)
    seroprevalence_plot = plot(date_x_points[1:length(seroprevalence_ts.median)],seroprevalence_ts.median; ribbon = seroprevalence_ts.ribbon, yformatter = y-> string(y.*100)* "\\%", xlabel = "time(days)",ylabel = "\\% recovered", label = "Model seroprevalence")
    scatter!(seroprevalence_plot,[june_30],[0.011];yerr = ([0.011 - 0.008],[0.013 - 0.011]), markersize = 4,markerstrokewidth = 1.0,  markerstrokealpha = 0.7,
    color = color_palette[2],label = "Observed seroprevalence")
    savefig(seroprevalence_plot,"$folder/seroprevalence_total.pdf")

    seroprevalence_data_labels = ["0-20", "20-60", "60+"]
    seroprevalence_data = [0.008,0.01,0.016] #from june 30 health ontario report
    seroprevalence_CI = (seroprevalence_data .- [0.003,0.007,0.011],[0.014,0.013,0.021] .- seroprevalence_data)
    display(seroprevalence_CI)
    prevalence_check_days = Dates.value(june_30 - start_date)
    seroprevalence_violins_youngest = dropdims(mapslices(sum,sol_data.recovered[:,prevalence_check_days,1:n_g-12],dims = (2)),dims = 2) ./ sum(population_demographics[1:n_g-12])
    seroprevalence_violins_middle = dropdims(mapslices(sum,sol_data.recovered[:,prevalence_check_days,n_g-12:n_g-4],dims = (2)),dims = 2)  ./ sum(population_demographics[n_g-12:n_g-4])
    seroprevalence_violins_oldest = dropdims(mapslices(sum,sol_data.recovered[:,prevalence_check_days,n_g-4:n_g],dims = (2)),dims = 2)  ./ sum(population_demographics[n_g-4:n_g])
    seroprevalence_violins_vector = vcat(seroprevalence_violins_youngest,seroprevalence_violins_middle,seroprevalence_violins_oldest)
    p = StatsPlots.violin(repeat(seroprevalence_data_labels,outer = length(particles)),seroprevalence_violins_vector; seriescolor = color_palette[1], label =  "Model seroprevalence",size = (400,300),xlabel = "Age groups",yformatter = y-> string(y)* "\\%")
    scatter!(p,seroprevalence_data_labels,seroprevalence_data; size = (400,600),yerr = seroprevalence_CI, markersize = 4,markerstrokewidth = 1.0,  markerstrokealpha = 0.7,
     color = color_palette[2],label =  "Observed seroprevalence",yformatter = y-> string(y*100)* "\\%",ylabel = "Seroprevalence on $june_30")
    savefig(p, "$folder/seroprevalence_by_age.pdf")


    # seroprevalence_plots = [plot() for i=1:n_g]
    # for i=1:n_g
    #     plot!(seroprevalence_plots[i],(recovered_ts.median[:,i] ./ population_demographics[i]) * 100,title = "Ages: $(age_bins[i])",xlabel = "time(days)",ylabel = "\\% recovered")
    # end
    # plot(seroprevalence_plots...,layout = (n_g,1),size = (300,2600))
    # savefig("$folder/seroprevalence_by_age.pdf")

end

#helper function to add letters to each subplot
function add_subplot_letters!(plot_list; pos = :top)
    for (i,sp) in enumerate(plot_list)
        letter = string(Char(i+96))
        if pos == :top
            annotate!(sp,xlims(sp)[1] + 0.02*((xlims(sp)[2]  - xlims(sp)[1])),ylims(sp)[2] - 0.11*((ylims(sp)[2]  - ylims(sp)[1])), Plots.text("$letter)", :left, 18))
        elseif pos == :bottom
            annotate!(sp,xlims(sp)[1] + 0.02*((xlims(sp)[2]  - xlims(sp)[1])),ylims(sp)[1] + 0.11*((ylims(sp)[2]  - ylims(sp)[1])), Plots.text("$letter)", :left, 18))
        end
    end
end
#this function plots the posterior distributions of the particle filtering
function violin_plots(folder,params_list::Vector,leading_eigenvector_strategy)
    #Violin plots of fitting output
    fit_params = (:ε_P_1,:κ_1,:κ_2,:c_1,:c_2,:γ_a,:R_0,:γ_s,:σ_0,:σ_1,:I_multiplier, :s,:ϕ,:P_multiplier, :E_multiplier,:η,:t_switch)
    titles = (L"\epsilon_{P_1}", L"\kappa_1",L"\kappa_2",L"c_1",L"c_2",L"\gamma_a",L"\\R_0",L"\gamma_s",L"\sigma_0",L"\sigma_1",L"\\I_0", L"\\s",L"\psi",L"\\P_0",L"\\E_0",L"\eta",L"t_{switch}")#L"\\beta",L"\\omega",L"\\delta", L"a_m")
    plot_objs = [plot() for p in fit_params]
    for (i,(p,title)) in enumerate(zip(fit_params,titles))
        p_dist = map(t->getfield(t,p),params_list)
        plot_objs[i] = StatsPlots.violin(p_dist,xlabel = "$title posterior",seriescolor = color_palette[1],legend =false,label = "Posterior distribution") #marker=(0.0,:blue,stroke(0))
    end 
    plot!(plot_objs[end], legend = true)
    #add_subplot_letters!(plot_objs)

    plot(plot_objs..., layout = (1,length(fit_params)),xticks = false, size = (2500,400))
    savefig("$folder/violins.pdf")
    # display(map(p->create_ascertainment_vectors(p)[1], params_list))
    ascertainment_vectors_1 = hcat(map(p->Array(create_ascertainment_vectors(p)[1]), params_list)...)
    #display(ascertainment_vectors_1)
    ascertainment_vectors_2 = hcat(map(p-> Array(create_ascertainment_vectors(p)[2]), params_list)...)
    labels = midpoints_age_bins
    
    p2 =  StatsPlots.violin(repeat(labels,inner = length(params_list)), vcat(collect(eachrow(ascertainment_vectors_1))...); seriescolor = color_palette[1], ylabel = "Ascertainment",rotation = 60, size = (400,300))
    
    savefig(p2,"$folder/ascertainment_1.pdf")

    r_vectors = hcat(map(p-> Array(create_susceptibility_vector(1.0,p)), params_list)...)
    display("means")
    display(mapslices(mean,r_vectors,dims = [2]))
    display("median")
    display(mapslices(median,r_vectors,dims = [2]))
    CIs = zip(mapslices(x -> quantile(x,0.05),r_vectors,dims = [2]),mapslices(x -> quantile(x, 0.95),r_vectors,dims = [2])) |> collect
    display("CIs")
    display(CIs)
    # display(create_susceptibility_vector(1.0,params_list[1]))
    p2 = StatsPlots.violin(repeat(labels,inner = length(params_list)), vcat(collect(eachrow(r_vectors))...),xlabel = "Age",ylabel = L"\rho_i",seriescolor = color_palette[1], rotation = 60, size = (400,300), label = "Susceptibility factor")
    savefig(p2,"$folder/r_vector.pdf")

    p3 = plot(labels, leading_eigenvector_strategy;xlabel = "Age",ylabel = "Fraction of vaccines",  title = "Contact-based vaccination strategy", size = (300,200), label = "fraction vaccinated under strategy")
    savefig(p3,"$folder/leading_eigenvector.pdf")
end
function plot_distribution_plane(folder,data_matrix,dates_to_use) #this function creates violin plots of the mortality, applied to TimeseriesData
 
    mortality_matrix = map(x->x.cumulative_infections_from_t_vac, data_matrix) 
    # display(data_matrix[1,1,1,1].cumulative_infections_from_t_vac)
    y_format_sci(x) = string(x) * "\\%"
    dim_to_axis = [t[1] => t[2] for t in zip(dimnames(data_matrix),axiskeys(data_matrix))] |> Dict
    threshold_key = :fixed_closure_value in keys(dim_to_axis) ? :fixed_closure_value : :dynamic_threshold
    threshold_labels =  dim_to_axis[threshold_key]
    scheme_labels = filter(x->x!="No vaccination", dim_to_axis[:vaccination_distribution])
    scheme_colors = color_palette[1:length(scheme_labels)]
    vac_labels = dim_to_axis[:ψ_0]
    vaccination_date_settings = dim_to_axis[:vaccination_begin_date][dates_to_use]
    num_particles = length(mortality_matrix[1,1,1,1][:,1])
    println("num particles = $num_particles")
    #define a function that computes the percentage reduction in mortality from no vaccination
    mortality_reduction(vac_scheme, v,t_vac,T) = [((mortality_matrix("No vaccination", v ,t_vac,T) .- mortality_matrix(vac_scheme, v ,t_vac,T)))[i] / mortality_matrix("No vaccination", v ,t_vac,T)[i] for i in 1:num_particles ]
    
    #print out some of the relevant mortality data for figure captions
    for scheme_label in scheme_labels, vac_rate in [0.01,0.015],t_vac in vaccination_date_settings
        arr = mortality_reduction(scheme_label,vac_rate,t_vac,2.0).*100 #|> AxisKeys.unname |> AxisKeys.keyless
        # display(arr)
        data_median =round(median(arr),sigdigits = 4)
        data_upper_quant = round(quantile(arr,0.95),sigdigits = 4)
        data_lower_quant = round(quantile(arr,0.05),sigdigits = 4)

        println("$scheme_label, ψ0 = $vac_rate, tvac = $t_vac :  $data_median %, ($data_upper_quant %,$data_lower_quant %)")
    end



    for t_vac in vaccination_date_settings, T in threshold_labels #loop over vac start dates and thresholds, printing mortality data for the no vaccination case

        arr = mortality_matrix("No vaccination",0.01 ,t_vac,T)
        data_median =round(median(arr),sigdigits = 6)
        data_upper_quant = round(quantile(arr,0.95),sigdigits = 6)
        data_lower_quant = round(quantile(arr,0.05),sigdigits = 6)

        println("no vaccination (total deaths), T = $(T * 100) %, tvac = $t_vac :  $data_median, ($data_upper_quant,$data_lower_quant)")
    end
    # display(table)
    # display(latexify(df; env = :table, latex = false))
    # m_oldest = mortality_matrix("Oldest first", 0.025,vaccination_date_settings[1],2.0)
    # println("Mortality at T = 200%, psi_0 = 2.5, jan, oldest only. Mean: $(mean(m_oldest)), stddev: $(std(m_oldest))")
    # m_uniform = mortality_matrix("Uniform", 0.025 ,vaccination_date_settings[2],2.0)
    # println("Mortality at T = 200%, psi_0 = 2.5, sept, uniform. Mean: $(mean(m_uniform)), stddev: $(std(m_uniform))")

    

    for t_vac in vaccination_date_settings
        vac_date_formatted = Dates.format(t_vac,"uuu d, yyyy")
        plot_matrix = [plot() for vac_rate in vac_labels]
        
        for (subplot, vac_rate) in zip(plot_matrix,vac_labels)

            labels_and_data = [ 
                (
                    100 * x_group,vac_group,color_label,
                    100* mortality_reduction(vac_group,vac_rate,t_vac,x_group)[i],
                    100 * median(mortality_reduction(vac_group,vac_rate,t_vac,x_group)),
                    100 * quantile(mortality_reduction(vac_group,vac_rate,t_vac,x_group), 0.95),
                    100 * quantile(mortality_reduction(vac_group,vac_rate,t_vac,x_group), 0.05)
                )
                 for  i in 1:num_particles, x_group in threshold_labels, (vac_group,color_label) in zip(scheme_labels,scheme_colors)
            ]
            #due to the strange way in which groupedboxplot/groupedviolinplot accepts data, I first organize all data tuples into an array, and then reshape all at once to ensure pairings are preserved
            label_list = reshape(labels_and_data,(:)) 

            x_groups = [l[1] for l in label_list]
            group_labels = [l[2] for l in label_list]
            color_labels = [l[3] for l in label_list]
            data_reshaped = [l[4] for l in label_list]
            medians = [l[5] for l in label_list]
            upper_quartile = [l[6] for l in label_list]
            lower_quartile = [l[7] for l in label_list]
            #unzip pairings into seperate vectors
           
            #median lines are just one point boxplots

            groupedviolin!(subplot,x_groups,data_reshaped;
                groups = group_labels, color = color_labels,legend = true,xformatter = x-> string(x) * "\\%", yformatter = y_format_sci,
                legendtitle = "$(trunc(vac_rate*100,digits = 2))% population vaccinated per week", ylabel = "\\% decrease in mortality",    
            )
            
            groupedboxplot!(subplot,x_groups,medians;
               groups = group_labels, label = false, color = color_labels,
               xlabel = "T, shutdown threshold (\\% of first wave)")

            groupedboxplot!(subplot,x_groups,upper_quartile;
               groups = group_labels, label = false,ylabel = "\\% decrease in mortality",
               xlabel = "T, shutdown threshold (\\% of first wave)", color = color_labels,linewidth = 1.0,
            )
            groupedboxplot!(subplot,x_groups,lower_quartile;
               groups = group_labels, label = false, ylabel = "\\% decrease in mortality",
               xlabel = "T, shutdown threshold (\\% of first wave)", color = color_labels,linewidth = 1.0,
            )



        end
        plot!(plot_matrix[begin];title = "Vaccine available: $vac_date_formatted")

        add_subplot_letters!(plot_matrix; pos = :top)
        p = plot(plot_matrix...,layout = (length(vac_labels),1),size = ( 500,1200))
        
        savefig("$folder/vaccination_by_mortality_boxplots_$(vac_date_formatted).pdf")
    end

    #this is the same as the above code, except over a subset of the data for the main text figure. Could this be organized better? definitely
    plot_matrix = [plot() for t_vac in vaccination_date_settings]
    for (subplot,t_vac) in zip(plot_matrix,vaccination_date_settings) 
            labels_and_data = [ 
                (
                    100 * x_group,vac_group,color_label,
                    100* mortality_reduction(vac_group,x_group,t_vac,2.0)[i],
                    100 * median(mortality_reduction(vac_group,x_group,t_vac,2.0)),
                    100 * quantile(mortality_reduction(vac_group,x_group,t_vac,2.0), 0.95),
                    100 * quantile(mortality_reduction(vac_group,x_group,t_vac,2.0), 0.05)
                )
                 for  i in 1:num_particles, x_group in vac_labels, (vac_group,color_label) in zip(scheme_labels,scheme_colors)
            ]
            
            label_list = reshape(labels_and_data,(:))
            x_groups = ["$(l[1])\\%" for l in label_list]
            group_labels = [l[2] for l in label_list]
            color_labels = [l[3] for l in label_list]
            data_reshaped = [l[4] for l in label_list]
            medians = [l[5] for l in label_list]
            upper_quartile = [l[6] for l in label_list]
            lower_quartile = [l[7] for l in label_list]


            groupedviolin!(subplot,x_groups,data_reshaped; xformatter = x-> string(x) * "\\%", yformatter = y_format_sci,
            groups = group_labels, color = color_labels,legend = false,  ylabel = "\\% decrease in mortality",
            )
            
            groupedboxplot!(subplot,x_groups,medians;
               groups = group_labels, label = false,legend = false, ylabel = "\\% decrease in mortality", color = color_labels,linewidth = 2,
            )
            groupedboxplot!(subplot,x_groups,upper_quartile;
               groups = group_labels, label = false,legend = false, ylabel = "\\% decrease in mortality", color = color_labels,linewidth = 1.0,
            )
            groupedboxplot!(subplot,x_groups,lower_quartile;
               groups = group_labels, label = false,legend = false, ylabel = "\\% decrease in mortality",
               xlabel = "Vaccination rate (\\% vaccinated per week)", color = color_labels,linewidth = 1.0,
            )

    end
    vac_dates_formatted = map(t_vac->Dates.format(t_vac,"uuu d, yyyy"),vaccination_date_settings)
    plot!(plot_matrix[1];title = "Vaccine available: $(vac_dates_formatted[1])",legend = true)

    plot!(plot_matrix[2];title = "Vaccine available: $(vac_dates_formatted[2])")

    add_subplot_letters!(plot_matrix)
    p = plot(plot_matrix...,layout = (length(vaccination_date_settings),1),size = (600,500))
    savefig("$folder/vaccination_by_mortality_small.pdf")
end

#Plot timeseries of vaccination data, this is applied to matricies of the type TimeseriesData
function plot_vaccination_timeseries(folder,data_matrix,infection_start,ontario_data) 
    color_palette = palette(:seaborn_muted) #color_palette
    dim_to_axis_labels = [t[1] => t[2] for t in zip(dimnames(data_matrix),axiskeys(data_matrix))] |> Dict #create dict of names of dimensions and their corresponding axes labels
    threshold_key = :dynamic_threshold
    threshold_labels =  dim_to_axis_labels[threshold_key]
    scheme_labels = filter(x->x!="No vaccination", dim_to_axis_labels[:vaccination_distribution]) #we are not plotting the "no vaccination" scheme in this figure
    scheme_colors = color_palette[1:length(scheme_labels)]
    vac_labels = dim_to_axis_labels[:ψ_0] #the vaccination labels are the labels corresponding to psi_0
    vaccination_date_settings = dim_to_axis_labels[:vaccination_begin_date]
    plot_folder = joinpath(folder,"vaccination_ts")

    window = 1000 #length of timeseries to plot
    #to plot months, need to specify manually the x points
    date_x_points = collect(ontario_data.start_date:Day(1):(ontario_data.start_date + Day(window - 1))) #get all the days from the infection start to the end of the plotting window, these will be our x points
    x_points_with_labels= collect(ontario_data.start_date:Month(6):(ontario_data.start_date + Day(window - 1))) #the labels for the x axis
    plot_x_ticks_labels = map(s-> "\\textrm{" * s* "}",Dates.format.(x_points_with_labels,"uuu d, yyyy")) #format the labels as strings

    main_text_ts_plot = typeof(plot())[] #setup a vector for the main text timeseries figure, which consists of a subset of the infection panels of the timeseries figure 
    @unpack infection_data_cumulative, distancing_data,infection_data_incident,workplace_closure = ontario_data #unpack all the covid data from the ontario data structure
    println(plot_folder)
    if !isdir(plot_folder) 
        mkdir(plot_folder)
    end
    for vaccination_begin_date in vaccination_date_settings,  ψ_0 in reverse(vac_labels), threshold_label in threshold_labels #we create one timeseries plot for each triple 
        
        plot_list = [plot() for i = 1:4]    
        plot_folder_threshold = joinpath(plot_folder,"threshold_$threshold_label")
        println(plot_folder_threshold)
        if !isdir(plot_folder_threshold) 
            mkdir(plot_folder_threshold)
        end
        sol_data = data_matrix(:, ψ_0, vaccination_begin_date, threshold_label)
        rec_max = []
        for (k,label) in enumerate(scheme_labels)
            cases_ts = get_stats(((sol_data(label).incident_infections[:,1:window]))) #
            distancing_ts = get_stats(sol_data(label).distancing[:,1:window])
            recovered_ts = get_stats((sol_data(label).recovered[:,1:window]))
            vac_ts = get_stats((sol_data(label).vaccinated[:,1:window]))
            work_shutdown_ts = get_stats(sol_data(label).work_shutdown[:,1:window])
            school_shutdown_ts = get_stats(sol_data(label).school_shutdown[:,1:window])
            # y_format(x) =string(Int(round(x,sigdigits=1, base = 10)))
            display(maximum(cases_ts.median .+ cases_ts.ribbon[2]))
            yticks_infections = collect(0.0:2500:ceil(maximum(cases_ts.median .+ cases_ts.ribbon[2])/2500)*2500 )
            yticks_labels_infections = map(string ∘ Int,yticks_infections)
            push!(rec_max,maximum(vac_ts.median))
                
            plot!(plot_list[1],date_x_points,cases_ts.median;  ribbon = cases_ts.ribbon, seriescolor = color_palette[k], label="Ascertained incident cases, $label",ylabel = "\\# individuals",yticks = (yticks_infections,yticks_labels_infections), ylims = (0, maximum(yticks_infections)))
            plot!(plot_list[2],date_x_points,distancing_ts.median; ribbon = distancing_ts.ribbon,seriescolor = color_palette[k], label = "x(t), $label", ylabel = "Fraction of baseline contacts")
            plot!(plot_list[3],date_x_points,work_shutdown_ts.median; seriescolor = color_palette[k], label = "Workplace contact reduction, $label", ylabel = "Fraction of baseline contacts")
            plot!(plot_list[4],date_x_points,recovered_ts.median; seriescolor = color_palette[k], label = "Recovered, $label", xlabel = "Time (date)", ylabel = "\\# individuals")
            plot!(plot_list[4],date_x_points,vac_ts.median; seriescolor = color_palette[k], label = "Immunized, $label", linestyle=:dash, xlabel = "Time (date)", ylabel = "\\# individuals")
            plot!(plot_list[3],date_x_points,school_shutdown_ts.median;   seriescolor = color_palette[k],linestyle = :dot, label = "School contact reduction, $label", ylabel = "Fraction of baseline contacts")
        end
        vac_date_formatted = Dates.format(vaccination_begin_date,"uuu d, yy")
        display(rec_max)
        step = floor(Int,(maximum(rec_max)/6)/250_000)*250_000
        yticks_recovered = collect(0.0:step:ceil(maximum(rec_max)/step)*step)
        yticks_labels_recovered = map(string ∘ Int,yticks_recovered)
        plot!(plot_list[1];title = "Vaccination begins on $vac_date_formatted, shutdown at $(threshold_label*100.0)\\% of first wave")
        plot!(plot_list[4]; ylims = (0, maximum(yticks_recovered)),yticks = (yticks_recovered,yticks_labels_recovered))
        vline!(plot_list[1],[vaccination_begin_date]; line = (:black, :dot), label = "vaccination begins")
        vline!(plot_list[2],[vaccination_begin_date], line = (:black, :dot),label = "vaccination begins")
        vline!(plot_list[3],[vaccination_begin_date], line = (:black, :dot),label = "vaccination begins")
        vline!(plot_list[4],[vaccination_begin_date], line = (:black, :dot),label = "vaccination begins,")
         if threshold_label == 2.0
            if ψ_0 in [0.005,0.015]
                panel_a = plot(deepcopy(plot_list[1]);layout = (1,1),legend = nothing, tickfontsize = 12,xticks = (x_points_with_labels,plot_x_ticks_labels),fillalpha = 0.2,linewidth = 2,title = "Vaccination begins on $vac_date_formatted, $(ψ_0*100.0)\\% of pop. vaccinated per week")
                push!(main_text_ts_plot,panel_a)
                display("push!")
            end
        end
        add_subplot_letters!(plot_list)
        p = plot(plot_list...,layout = (length(plot_list),1),size = (800,1000),legend = :outerright, tickfontsize = 12,xticks = (x_points_with_labels,plot_x_ticks_labels),fillalpha = 0.2,linewidth = 2)
      
        savefig(p,"$plot_folder_threshold/plot_$(ψ_0)_$(vac_date_formatted)_$(threshold_label).pdf")
    end
    plot!(main_text_ts_plot[1]; legend = :topright)
    plot!(main_text_ts_plot[end];xlabel = "Time (date)")

    add_subplot_letters!(main_text_ts_plot)
    
    main_text_plot = plot(main_text_ts_plot[1:2]...;layout = (2,1),size = (700,400), tickfontsize = 12,xticks = (x_points_with_labels,plot_x_ticks_labels),fillalpha = 0.2,linewidth = 2)
    savefig(main_text_plot,"$folder/main_text_ts_1.pdf")

    main_text_plot = plot(main_text_ts_plot...;layout = (length(main_text_ts_plot),1),size = (800,1000), tickfontsize = 12,xticks = (x_points_with_labels,plot_x_ticks_labels),fillalpha = 0.2,linewidth = 2)
    savefig(main_text_plot,"$folder/main_text_ts.pdf")

    plot!(main_text_ts_plot[3]; legend = :topright)
    main_text_plot = plot(main_text_ts_plot[3:4]...;layout = (2,1),size = (700,400), tickfontsize = 12,xticks = (x_points_with_labels,plot_x_ticks_labels),fillalpha = 0.2,linewidth = 2)
    savefig(main_text_plot,"$folder/main_text_ts_2.pdf")

end
function vaccination_analysis_plots(folder,vaccination_analysis_data,dates_to_use)
    dim_to_axis = [t[1] => t[2] for t in zip(dimnames(vaccination_analysis_data),axiskeys(vaccination_analysis_data))] |> Dict
    threshold_key = :fixed_closure_value in keys(dim_to_axis) ? :fixed_closure_value : :dynamic_threshold
    threshold_labels =  dim_to_axis[threshold_key]

    scheme_labels = filter(x->x!="No vaccination", dim_to_axis[:vaccination_distribution])
    vac_labels = dim_to_axis[:ψ_0]
    colors = color_palette[1:length(scheme_labels)]
    vaccination_date_settings = dim_to_axis[:vaccination_begin_date][dates_to_use]
    
    plots = [plot() for v in vaccination_date_settings]

    for (subplot,t_vac) in zip(plots,vaccination_date_settings)
        for (color,vaccination_distribution) in zip(colors,scheme_labels)
            data = map(d -> d.cumulative_infections_at_end, vaccination_analysis_data(vaccination_distribution = vaccination_distribution,ψ_0 = :, dynamic_threshold = threshold_labels[end-1] ,vaccination_begin_date = t_vac))
            no_vac_data = map(d -> d.cumulative_infections_at_end, vaccination_analysis_data(vaccination_distribution = "No vaccination",ψ_0 = :, dynamic_threshold = threshold_labels[end-1] ,vaccination_begin_date = t_vac))
            relative_mortality =  [(n_v_d .- d)./n_v_d * 100 for (d,n_v_d) in zip(data,no_vac_data)]
            plot!(subplot,100 .* vac_labels,mean.(relative_mortality); yerror = std.(relative_mortality),
             yformatter = x-> string(x) * "\\%",
             seriescolor = color,label = vaccination_distribution,
             ylabel = "decrease in mortality",markerstrokecolor=:auto,xformatter = x-> string(x) * "\\%")
        end
    end
    add_subplot_letters!(plots)
    vac_dates_formatted = map(t_vac->Dates.format(t_vac,"uuu d, yyyy"),vaccination_date_settings)
    plot!(plots[1];title = "Vaccine available: $(vac_dates_formatted[1])",legend = true)

    plot!(plots[2];title = "Vaccine available: $(vac_dates_formatted[2])",xlabel = "\\% of population vaccinated per week")
    p = plot(plots...,layout = (length(vaccination_date_settings),1),size = (400,500))
    savefig("$folder/mortality_by_vaccine_availablility.pdf")
end


function R_0_plots(folder,analysis_data,dates_to_use)
    dim_to_axis = [t[1] => t[2] for t in zip(dimnames(analysis_data),axiskeys(analysis_data))] |> Dict
    threshold_key = :fixed_closure_value in keys(dim_to_axis) ? :fixed_closure_value : :dynamic_threshold
    threshold_labels =  dim_to_axis[threshold_key]

    scheme_labels = filter(x->x!="No vaccination", dim_to_axis[:vaccination_distribution])
    vac_labels = dim_to_axis[:ψ_0]
    colors = color_palette[1:length(scheme_labels)]
    vaccination_date_settings = dim_to_axis[:vaccination_begin_date][dates_to_use]

    R_0_settings = dim_to_axis[:r] #I know these don't match
    
    plots = [plot() for v in vaccination_date_settings]

    for (subplot,t_vac) in zip(plots,vaccination_date_settings)
        for (color,vaccination_distribution) in zip(colors,scheme_labels)
            data = map(d -> d.cumulative_infections_at_end, analysis_data(r = :, vaccination_distribution = vaccination_distribution,ψ_0 = only(vac_labels), dynamic_threshold = only(threshold_labels) ,vaccination_begin_date = t_vac))
            no_vac_data = map(d -> d.cumulative_infections_at_end, analysis_data(r = :,vaccination_distribution = "No vaccination",ψ_0 = only(vac_labels), dynamic_threshold = only(threshold_labels) ,vaccination_begin_date = t_vac))
            relative_mortality =  [(n_v_d .- d)./n_v_d * 100 for (d,n_v_d) in zip(data,no_vac_data)]
            plot!(subplot,R_0_settings,mean.(relative_mortality); yerror = std.(relative_mortality),
             seriescolor = color,label = vaccination_distribution,
             yformatter = x-> string(x) * "\\%",
             ylabel = "\\% decrease in mortality",markerstrokecolor=:auto)
        end
    end
    vac_dates_formatted = map(t_vac->Dates.format(t_vac,"uuu d, yyyy"),vaccination_date_settings)
    plot!(plots[1];title = "Vaccine available: $(vac_dates_formatted[1])",legend = true)

    plot!(plots[2];title = "Vaccine available: $(vac_dates_formatted[2])",xlabel = L"\\R_0")
    add_subplot_letters!(plots,pos = :bottom)
    p = plot(plots...,layout = (length(vaccination_date_settings),1),size = (400,500))
    savefig("$folder/mortality_by_R_0.pdf")
end


function vac_date_plots(folder,analysis_data,infection_start_date)
    dim_to_axis = [t[1] => t[2] for t in zip(dimnames(analysis_data),axiskeys(analysis_data))] |> Dict
    threshold_key = :fixed_closure_value in keys(dim_to_axis) ? :fixed_closure_value : :dynamic_threshold
    threshold_labels =  dim_to_axis[threshold_key]

    scheme_labels = filter(x->x!="No vaccination", dim_to_axis[:vaccination_distribution])
    vac_rates = [dim_to_axis[:ψ_0][2],dim_to_axis[:ψ_0][end-1]]
    colors = color_palette[1:length(scheme_labels)]
    vaccination_date_settings = dim_to_axis[:vaccination_begin_date]
    vaccination_days_from_start = Dates.value.(vaccination_date_settings .- infection_start_date)
    display(vaccination_days_from_start)

    for threshold_label in threshold_labels
        plots = [plot() for v in vac_rates]
        for (subplot,vac_rate) in zip(plots,vac_rates)
            recovered_series_data = []
            for (i,(days_from_start,vac_date)) in enumerate(zip(vaccination_days_from_start,vaccination_date_settings))
                data_pt = analysis_data(vaccination_distribution = scheme_labels[2],ψ_0 = vac_rate, dynamic_threshold = threshold_label ,vaccination_begin_date = vac_date)     
                #display(mean.(eachcol(data_pt.recovered[:,vaccination_days_from_start]))./population)      
                push!(recovered_series_data,median(data_pt.recovered[:,days_from_start]))
            end
            display((vac_rate,recovered_series_data))
            for (color,vaccination_distribution) in zip(colors,scheme_labels)
                data = map(d -> d.cumulative_infections_at_end, analysis_data(vaccination_distribution = vaccination_distribution,ψ_0 = vac_rate, dynamic_threshold = threshold_label ,vaccination_begin_date = :))
                no_vac_data = map(d -> d.cumulative_infections_at_end, analysis_data(vaccination_distribution = "No vaccination",ψ_0 = vac_rate, dynamic_threshold = threshold_label ,vaccination_begin_date = :))
                relative_mortality =  [(n_v_d .- d)./n_v_d * 100 for (d,n_v_d) in zip(data,no_vac_data)]
                
                xlabels = map(s-> "\\textrm{" * Dates.format(s[1],"uuu d, yyyy") * "| R: "* string(trunc((s[2]/population)*100,digits =1)) * "\\% pop.}",zip(vaccination_date_settings,recovered_series_data))
                #display(xlabels)
                plot!(subplot,xlabels,mean.(relative_mortality); yerror = std.(relative_mortality), seriescolor = color, 
                label = "Mortality: $vaccination_distribution", xrotation = 45,markerstrokecolor=:auto,
                yformatter = x-> string(x) * "\\%")
                #plot!(subplot,vaccination_date_settings,mean.(recovered); linestyle = :dot,seriescolor = color,label = "Recovered")
            end
            plot!(subplot;title = "$(trunc(vac_rate*100,digits = 2))% of population vaccinated per week",legend = true,ylabel = "\\% decrease in mortality")
        end
        plot!(plots[1]; xticklabels = false,xticks = false)
        add_subplot_letters!(plots,pos = :bottom)
        p = plot(plots...,layout = (length(vac_rates),1),size = (600,400))
        savefig("$folder/mortality_by_vac_date_$threshold_label.pdf")
    end
end

function bivariate_heatmaps(folder,vaccination_output,start_date,dates_to_use)
    dim_to_axis = [t[1] => t[2] for t in zip(dimnames(vaccination_output),axiskeys(vaccination_output))] |> Dict
    threshold_key = :dynamic_threshold
    threshold_labels = dim_to_axis[threshold_key][2:end]
    vac_schemes = filter(x->x!="No vaccination", dim_to_axis[:vaccination_distribution])
    display(vac_schemes)
    vac_rates = dim_to_axis[:ψ_0][2:end]
    vaccination_date_settings =dim_to_axis[:vaccination_begin_date][dates_to_use]
    display(vaccination_date_settings)
    plot_list = [plot() for t in vaccination_date_settings]
    histogram_data = [zeros(0) for i in vac_schemes]
    
    for (subplot,vac_date) in zip(plot_list,vaccination_date_settings)

        vac_date_formatted = Dates.format(vac_date,"uuu d, yyyy")
        heatmap_grid = zeros(Int64,length(threshold_labels),length(vac_rates))
        vaccination_begin_index = Dates.value(vac_date - start_date)
        for (i,threshold_label) in enumerate(threshold_labels), (j,vac_rate) in enumerate(vac_rates)
            mortality_by_scheme = map(d -> d.cumulative_infections_at_end, vaccination_output(vaccination_distribution = vac_schemes, ψ_0 = vac_rate, dynamic_threshold = threshold_label ,vaccination_begin_date = vac_date))
            recovered_by_scheme = map(d -> d.recovered[:,vaccination_begin_index], vaccination_output(vaccination_distribution = vac_schemes, ψ_0 = vac_rate, dynamic_threshold = threshold_label ,vaccination_begin_date = vac_date))
           # display((threshold_label,vac_rate))
            list_of_best_vac_scheme = map(pt -> argmin(pt),zip(mortality_by_scheme...))
        #    display(mean.(mortality_by_scheme))
        #    display(list_of_best_vac_scheme)
            mean_mortality_by_scheme = mean.(mortality_by_scheme)
            best_vac_scheme = argmin(mean_mortality_by_scheme)#mode(list_of_best_vac_scheme)
            heatmap_grid[i,j] = best_vac_scheme
            list_of_recovered_by_scheme = zip(recovered_by_scheme...)
            for (best_vac_scheme_pt,recovered_by_scheme_pt) in zip(list_of_best_vac_scheme,list_of_recovered_by_scheme)
                push!(histogram_data[best_vac_scheme_pt],recovered_by_scheme_pt[best_vac_scheme_pt])
            end
        end


        vac_rates_as_percentage = (100 .* vac_rates)
        threshold_labels_as_percentage = (100 .* threshold_labels)
        if minimum(heatmap_grid) == maximum(heatmap_grid) #need a ridiculous hack for single color heatmap
            grid_palette = cgrad(color_palette[minimum(heatmap_grid):minimum(heatmap_grid)+1])#need two colors
            new_heatmap_grid = fill(10,length(threshold_labels),length(vac_rates) + 1)
            new_heatmap_grid[:,1:end-1] .= heatmap_grid
            new_vac_rates = [vac_rates_as_percentage..., 10] 
            heatmap!(subplot,new_vac_rates,threshold_labels_as_percentage,new_heatmap_grid; seriescolors = grid_palette, xformatter = x-> string(x) * "\\%", yformatter = x-> string(x) * "\\%",
            ylims = (minimum(threshold_labels_as_percentage),maximum(threshold_labels_as_percentage)),xlims = (minimum(vac_rates_as_percentage),maximum(vac_rates_as_percentage)),
            ylabel = "T, shutdown threshold (\\% active cases)", title = "Vaccine available on $vac_date_formatted", colorbar = false)
        else
            grid_palette = cgrad(color_palette[minimum(heatmap_grid):maximum(heatmap_grid)])

            heatmap!(subplot,vac_rates_as_percentage,threshold_labels_as_percentage,heatmap_grid; seriescolors = grid_palette, xformatter = x-> string(x) * "\\%", yformatter = x-> string(x) * "\\%",
            ylims = (minimum(threshold_labels_as_percentage),maximum(threshold_labels_as_percentage)),xlims = (minimum(vac_rates_as_percentage),maximum(vac_rates_as_percentage)),
            ylabel = "T, shutdown threshold  (\\% active cases)", title = "Vaccine available on $vac_date_formatted", colorbar = false )
        end

        # end
    end
    plot!(plot_list[end];xlabel = "Vaccination rate (\\% per week)")
    label_heatmap!(plot_list[begin],color_palette[1:length(vac_schemes)],vac_schemes)
    add_subplot_letters!(plot_list)
    p1 = plot(plot_list...,layout = (length(dates_to_use),1),size = (400,300*length(dates_to_use)))
    savefig(p1,"$folder/bivariate_heatmap.pdf")
    plot_list=[]
    histogram_bins = range(minimum(minimum.(filter(x -> !isempty(x), histogram_data))),stop = maximum(maximum.(filter(x -> !isempty(x), histogram_data))), length = 30)
    total_histogram_points = sum(length.(histogram_data))
    for (histogram_data_for_vac_scheme,label,color) in filter(x -> !isempty(x[1]), collect(zip(histogram_data,vac_schemes,color_palette))) 
        hist = fit(Histogram,histogram_data_for_vac_scheme,histogram_bins)
        histogram_scaled = Histogram(hist.edges, hist.weights ./ total_histogram_points)  #normalize by total number of observations across all vac scheme histograms
        # display(histogram_data_for_vac_scheme)
        hist_mean = median(histogram_data_for_vac_scheme)
        # display(hist_mean)
        p = plot(histogram_scaled,seriescolor = color,title = label,xlabel = "Recovered at vaccination begin date", legend = false,dpi = 300, xformatter = x-> string(trunc(100 * x/population,digits = 2)) * "\\%")
        vline!(p,[hist_mean]; linestyle = :dash, linewidth = 2.0, linecolor = :black)
        # annotate!(hist_mean, ylims(p)[2]*0.9, text("mean",12 ,:left))
        push!(plot_list,p)
    end
    add_subplot_letters!(plot_list)
    histograms_plot = plot(plot_list...,layout = (2 + length(plot_list) % 2,2 - length(plot_list) % 2),size = (350*(2 - length(plot_list) % 2),200*(2 + length(plot_list) % 2)),ylabel = "No. parameter sets")
    savefig(histograms_plot,"$folder/histograms.pdf")
    return p1,histograms_plot
end


function bivariate_heatmaps_vac_efficacy(folder,vaccination_output,start_date,dates_to_use)    
    dim_to_axis = [t[1] => t[2] for t in zip(dimnames(vaccination_output),axiskeys(vaccination_output))] |> Dict
    threshold_labels = dim_to_axis[:dynamic_threshold]
    efficacy_elderly_labels =  dim_to_axis[:vaccination_efficacy_elderly]
    efficacy_else_labels = dim_to_axis[:vaccination_efficacy_else]
    vac_schemes = filter(x->x != "No vaccination", dim_to_axis[:vaccination_distribution])
    display(vac_schemes)
    vac_rates = [dim_to_axis[:ψ_0][end]]
    vaccination_date_settings = dim_to_axis[:vaccination_begin_date][dates_to_use]
    display(vaccination_date_settings)
    plot_list = [plot() for t in vaccination_date_settings]
    histogram_data = [zeros(0) for i in vac_schemes]
    
    for (subplot,vac_date) in zip(plot_list,vaccination_date_settings)

        vac_date_formatted = Dates.format(vac_date,"uuu d, yyyy")
        heatmap_grid = zeros(Int64,length(efficacy_elderly_labels),length(efficacy_else_labels))
        vaccination_begin_index = Dates.value(vac_date - start_date)
        for (i,efficacy_elderly) in enumerate(efficacy_elderly_labels), (j,efficacy_else) in enumerate(efficacy_else_labels)
            mortality_by_scheme = map(d -> d.cumulative_infections_at_end, vaccination_output(
                vaccination_distribution = vac_schemes,
                ψ_0 = only(vac_rates), 
                dynamic_threshold = only(threshold_labels),
                vaccination_begin_date = vac_date, 
                vaccination_efficacy_elderly = efficacy_elderly,
                vaccination_efficacy_else = efficacy_else
            ))
            recovered_by_scheme = map(d -> d.recovered[:,vaccination_begin_index],vaccination_output(
                vaccination_distribution = vac_schemes,
                ψ_0 = only(vac_rates), 
                dynamic_threshold = only(threshold_labels),
                vaccination_begin_date = vac_date, 
                vaccination_efficacy_elderly = efficacy_elderly,
                vaccination_efficacy_else = efficacy_else
            ))
           # display((threshold_label,vac_rate))
           list_of_best_vac_scheme = map(pt -> argmin(pt),zip(mortality_by_scheme...))
           #    display(mean.(mortality_by_scheme))
           #    display(list_of_best_vac_scheme)
               mean_mortality_by_scheme = mean.(mortality_by_scheme)
               best_vac_scheme = argmin(mean_mortality_by_scheme)#mode(list_of_best_vac_scheme)
               heatmap_grid[i,j] = best_vac_scheme
            list_of_recovered_by_scheme = zip(recovered_by_scheme...)
            for (best_vac_scheme_pt,recovered_by_scheme_pt) in zip(list_of_best_vac_scheme,list_of_recovered_by_scheme)
                push!(histogram_data[best_vac_scheme_pt],recovered_by_scheme_pt[best_vac_scheme_pt])
            end
        end


        elderly_efficacy_as_percentage = (100 .* efficacy_elderly_labels)
        else_efficacy_as_percentage = (100 .* efficacy_else_labels)
        if minimum(heatmap_grid) == maximum(heatmap_grid) #need a ridiculous hack for single color heatmap
            grid_palette = cgrad(color_palette[minimum(heatmap_grid):minimum(heatmap_grid)+1])#need two colors
            display(heatmap_grid)
            new_heatmap_grid = fill(10,length(efficacy_elderly_labels),length(efficacy_else_labels) + 1)
            new_heatmap_grid[:,1:end-1] .= heatmap_grid
            new_vac_rates = [else_efficacy_as_percentage..., 100] 
            heatmap!(subplot,new_vac_rates,elderly_efficacy_as_percentage,new_heatmap_grid; seriescolors = grid_palette, xformatter = x-> string(x) * "\\%", yformatter = x-> string(x) * "\\%",
            ylims = (minimum(elderly_efficacy_as_percentage),maximum(elderly_efficacy_as_percentage)),xlims = (minimum(else_efficacy_as_percentage),maximum(else_efficacy_as_percentage)),
            ylabel = "Vaccination efficacy (age >60)", title = "Vaccine available on $vac_date_formatted", colorbar = false)
        else
            grid_palette = cgrad(color_palette[minimum(heatmap_grid):maximum(heatmap_grid)])
            display(heatmap_grid)
            heatmap!(subplot,else_efficacy_as_percentage,elderly_efficacy_as_percentage,heatmap_grid; seriescolors = grid_palette, xformatter = x-> string(x) * "\\%", yformatter = x-> string(x) * "\\%",
            ylims = (minimum(elderly_efficacy_as_percentage),maximum(elderly_efficacy_as_percentage)),xlims = (minimum(else_efficacy_as_percentage),maximum(else_efficacy_as_percentage)),
            ylabel = "Vaccination efficacy (age >60)", title = "Vaccine available on $vac_date_formatted", colorbar = false )
        end

        # end
    end
    plot!(plot_list[end];xlabel = "Vaccination efficacy (age <60)")
    label_heatmap!(plot_list[begin],color_palette[1:length(vac_schemes)],vac_schemes)
    add_subplot_letters!(plot_list)
    p1 = plot(plot_list...,layout = (length(dates_to_use),1),size = (400,300*length(dates_to_use)))
    savefig(p1,"$folder/bivariate_heatmap_efficacy.pdf")
    return p1
end


function bivariate_heatmaps_x_vac(folder,vaccination_output,start_date,dates_to_use)    
    dim_to_axis = [t[1] => t[2] for t in zip(dimnames(vaccination_output),axiskeys(vaccination_output))] |> Dict
    display(axiskeys(vaccination_output))
    threshold_labels = dim_to_axis[:dynamic_threshold]
    kappa_vac_labels =  dim_to_axis[:κ_vac]
    c_vac_labels = dim_to_axis[:c_vac]
    vac_schemes = filter(x->x!="No vaccination", dim_to_axis[:vaccination_distribution])
    display(c_vac_labels)
    display(kappa_vac_labels)
    vac_rates = dim_to_axis[:ψ_0]
    x_vac_0 = only(dim_to_axis[:x_vac_0])
    vaccination_date_settings = dim_to_axis[:vaccination_begin_date][dates_to_use]
    plot_list = [plot() for t in vaccination_date_settings]
    histogram_data = [zeros(0) for i in vac_schemes]
    
    for (subplot,vac_date) in zip(plot_list,vaccination_date_settings)

        vac_date_formatted = Dates.format(vac_date,"uuu d, yyyy")
        heatmap_grid = zeros(Int64,length(kappa_vac_labels),length(c_vac_labels))
        vaccination_begin_index = Dates.value(vac_date - start_date)
        for (i,κ_vac) in enumerate(kappa_vac_labels), (j,c_vac) in enumerate(c_vac_labels)
            # display(c_vac)
            mortality_by_scheme = map(d -> d.cumulative_infections_from_t_vac, vaccination_output(
                vaccination_distribution = vac_schemes,
                ψ_0 = only(vac_rates), 
                dynamic_threshold = only(threshold_labels),
                vaccination_begin_date = vac_date, 
                κ_vac = κ_vac,
                c_vac = c_vac,
                x_vac_0 = x_vac_0,
            ))
            recovered_by_scheme = map(d -> d.recovered[:,vaccination_begin_index],vaccination_output(
                vaccination_distribution = vac_schemes,
                ψ_0 = only(vac_rates), 
                dynamic_threshold = only(threshold_labels),
                vaccination_begin_date = vac_date, 
                κ_vac = κ_vac,
                c_vac = c_vac,
                x_vac_0 = x_vac_0,
            ))
            list_of_best_vac_scheme = map(pt -> argmin(pt),zip(mortality_by_scheme...))
            #    display(mean.(mortality_by_scheme))
            #    display(list_of_best_vac_scheme)
                mean_mortality_by_scheme = mean.(mortality_by_scheme)
                best_vac_scheme = argmin(mean_mortality_by_scheme)#mode(list_of_best_vac_scheme)
                heatmap_grid[i,j] = best_vac_scheme
            list_of_recovered_by_scheme = zip(recovered_by_scheme...)
            for (best_vac_scheme_pt,recovered_by_scheme_pt) in zip(list_of_best_vac_scheme,list_of_recovered_by_scheme)
                push!(histogram_data[best_vac_scheme_pt],recovered_by_scheme_pt[best_vac_scheme_pt])
            end
        end
        display(heatmap_grid)
        if minimum(heatmap_grid) == maximum(heatmap_grid) #need a ridiculous hack for single color heatmap
            grid_palette = cgrad(color_palette[minimum(heatmap_grid):minimum(heatmap_grid)+1])#need two colors
            new_heatmap_grid = fill(100,length(kappa_vac_labels),length(c_vac_labels) + 1)
            new_heatmap_grid[:,1:end-1] .= heatmap_grid
            new_vac_rates = [c_vac_labels..., 100] 
            heatmap!(subplot,new_vac_rates,kappa_vac_labels,new_heatmap_grid; seriescolors = grid_palette,
            ylims = (minimum(kappa_vac_labels),maximum(kappa_vac_labels)),xlims = (minimum(c_vac_labels),maximum(c_vac_labels)),
        ylabel =L"\kappa_{vac}"*", vacc. social learning rate", title = "Vaccine available on $vac_date_formatted", colorbar = false)
        else
            grid_palette = cgrad(color_palette[minimum(heatmap_grid):maximum(heatmap_grid)])
            heatmap!(subplot,string.(c_vac_labels),string.(kappa_vac_labels),heatmap_grid; seriescolors = grid_palette,
            #ylims = (minimum(kappa_vac_labels),maximum(kappa_vac_labels)),#xlims = (minimum(c_vac_labels),maximum(c_vac_labels)),
            ylabel = L"\kappa_{vac}"*", vacc. social learning rate", title = "Vaccine available on $vac_date_formatted", colorbar = false )
        end

        # end
    end
    plot!(plot_list[end];xlabel = L"c_{vac}" *", incentive not to vaccinate" )
    label_heatmap!(plot_list[begin],color_palette[1:length(vac_schemes)],vac_schemes)
    add_subplot_letters!(plot_list)
    p1 = plot(plot_list...,layout = (length(dates_to_use),1),size = (600,300*length(dates_to_use)))
    savefig(p1,"$folder/bivariate_heatmap_x_vac.pdf")
end




function plot_x_vac_timeseries(folder,data_matrix,infection_start,ontario_data)
    color_palette = palette(:seaborn_muted)
    dim_to_axis_labels = [t[1] => t[2] for t in zip(dimnames(data_matrix),axiskeys(data_matrix))] |> Dict
    threshold_key = :dynamic_threshold
    threshold_label = only(dim_to_axis_labels[threshold_key])
    scheme_labels = filter(x->x!="No vaccination", dim_to_axis_labels[:vaccination_distribution])
    scheme_colors = color_palette[1:length(scheme_labels)]
    vac_label = only(dim_to_axis_labels[:ψ_0])
    vaccination_date_settings = dim_to_axis_labels[:vaccination_begin_date]
    plot_folder = joinpath(folder,"vaccine_refusal_ts")
    # display(dim_to_axis_labels)
    kappa_vac_labels =  dim_to_axis_labels[:κ_vac]
    c_vac_labels = dim_to_axis_labels[:c_vac]

    x_vac_0 = only(dim_to_axis_labels[:x_vac_0])
    window = 1000

    date_x_points = collect(ontario_data.start_date:Day(1):(ontario_data.start_date + Day(window-1)))
    x_points_with_labels= collect(ontario_data.start_date:Month(6):(ontario_data.start_date + Day(window-1)))
    plot_x_ticks_labels = map(s-> "\\textrm{" * s* "}",Dates.format.(x_points_with_labels,"uuu d, yyyy"))

    @unpack infection_data_cumulative, distancing_data,infection_data_incident,workplace_closure = ontario_data
    println(plot_folder)
    if !isdir(plot_folder) 
        mkdir(plot_folder)
    end
    for kappa_vac_label in kappa_vac_labels, vaccination_begin_date in vaccination_date_settings, c_vac_label in c_vac_labels
        vac_date_formatted = Dates.format(vaccination_begin_date,"uuu d, yyyy")
        
        plot_list = [plot() for i = 1:4]    
        plot_folder_threshold = joinpath(plot_folder,"threshold_$threshold_label")
        println(plot_folder_threshold)
        if !isdir(plot_folder_threshold) 
            mkdir(plot_folder_threshold)
        end
        sol_data = data_matrix(x_vac_0,:, vac_label, vaccination_begin_date, threshold_label,c_vac_label,kappa_vac_label)
        rec_max = []
        for (k,label) in enumerate(scheme_labels)
            cases_ts = get_stats(((sol_data(label).incident_infections[:,1:window]))) #
            distancing_ts = get_stats(sol_data(label).distancing[:,1:window])
            recovered_ts = get_stats((sol_data(label).recovered[:,1:window]))
            vac_ts = get_stats((sol_data(label).vaccinated[:,1:window]))
            work_shutdown_ts = get_stats(sol_data(label).work_shutdown[:,1:window])
            school_shutdown_ts = get_stats(sol_data(label).school_shutdown[:,1:window])
            vac_refusal_ts = get_stats(sol_data(label).vaccine_refusal[:,1:window])
            yticks_infections = collect(0.0:2500:ceil(maximum(cases_ts.median .+ cases_ts.ribbon[2])/2500)*2500 )
            yticks_labels_infections = map(string ∘ Int,yticks_infections)
            push!(rec_max,maximum(vac_ts.median))
            plot!(plot_list[2],date_x_points,vac_refusal_ts.median, seriescolor = color_palette[k],linestyle = :dot, label = "y(t), fraction of vaccinators, $label")
            plot!(plot_list[1],date_x_points,cases_ts.median;  ribbon = cases_ts.ribbon, seriescolor = color_palette[k], label="Ascertained incident cases, $label",ylabel = "\\# individuals",ylim = (0, maximum(yticks_infections)), yticks = (yticks_infections,yticks_labels_infections))
            plot!(plot_list[2],date_x_points,distancing_ts.median; ribbon = distancing_ts.ribbon,seriescolor = color_palette[k], label = "x(t), $label", ylabel = "Fraction of baseline contacts")
            plot!(plot_list[3],date_x_points,work_shutdown_ts.median; seriescolor = color_palette[k], label = "Workplace contact reduction, $label", ylabel = "Fraction of baseline contacts")
            plot!(plot_list[4],date_x_points,recovered_ts.median; seriescolor = color_palette[k], label = "Recovered, $label", xlabel = "Time (date)", ylabel = "\\# individuals")
            plot!(plot_list[4],date_x_points,vac_ts.median; seriescolor = color_palette[k], label = "Immunized, $label", linestyle=:dash, xlabel = "Time (date)", ylabel = "\\# individuals")
            plot!(plot_list[3],date_x_points,school_shutdown_ts.median;   seriescolor = color_palette[k],linestyle = :dot, label = "School contact reduction, $label", ylabel = "Fraction of baseline contacts")
        end

        # scatter!(plot_list[2],workplace_closure, label = "Workplace closure data",  markersize = 2,markerstrokewidth = 1,  markerstrokealpha = 0.7, color = color_palette[3])
        # scatter!(plot_list[2],distancing_data,color = color_palette[2], label = "Distancing data",  markersize = 2,markerstrokewidth = 1,  markerstrokealpha = 0.7)
        # scatter!(plot_list[1],sum.(eachrow(infection_data_incident)),label = "Case notification data",  markersize = 2, markerstrokewidth = 1, markerstrokealpha = 0.7, color = color_palette[1])
        vac_date_formatted = Dates.format(vaccination_begin_date,"uuu d, yy")

        step = floor(Int,(maximum(rec_max)/6)/250_000)*250_000
        yticks_recovered = collect(0.0:step:ceil(maximum(rec_max)/step)*step)
        yticks_labels_recovered = map(string ∘ Int,yticks_recovered)
        t_vac_index = vaccination_begin_date
        plot!(plot_list[1]; title = "Vaccination begins on $vac_date_formatted, shutdown at $(threshold_label*100.0)% of first wave", yformatter = x -> string(x)*"\\%")
        plot!(plot_list[4]; yticks = (yticks_recovered,yticks_labels_recovered), ylim = maximum(rec_max))
        vline!(plot_list[1],[t_vac_index]; line = (:black, :dot), label = "vaccination begins")
        vline!(plot_list[2],[t_vac_index], line = (:black, :dot),label = "vaccination begins")
        vline!(plot_list[3],[t_vac_index], line = (:black, :dot),label = "vaccination begins")
        vline!(plot_list[4],[t_vac_index], line = (:black, :dot),label = "vaccination begins")
        add_subplot_letters!(plot_list)
        p = plot(plot_list...,layout = (length(plot_list),1),size = (800,1000),legend = :outerright, tickfontsize = 12,xticks = (x_points_with_labels,plot_x_ticks_labels))

        savefig(p,"$plot_folder_threshold/plot_$(c_vac_label)_$(vac_date_formatted)_$(kappa_vac_label).pdf")
    end
end

function label_heatmap!(p,colors,labels; legend_title = nothing)
    for (color,label) in zip(colors,labels)
        plot!(p,[0.001, 0.001],fillrange =0.001,seriescolor = color,
        label = label, widen = false,legend = true, legendtitle = legend_title)
    end
end


function v_T_plots(folder,analysis_data,dates_to_use,mean_particle)
    dim_to_axis = Dict(t[1] => t[2] for t in zip(dimnames(analysis_data),axiskeys(analysis_data)))
    threshold_key = :fixed_closure_value in keys(dim_to_axis) ? :fixed_closure_value : :dynamic_threshold
    threshold_labels =  dim_to_axis[threshold_key]

    scheme_labels = filter(x->x!="No vaccination", dim_to_axis[:vaccination_distribution])
    vac_labels = dim_to_axis[:ψ_0]
    colors = color_palette[1:length(scheme_labels)]
    vaccination_date_settings = dim_to_axis[:vaccination_begin_date][dates_to_use]

    v_D_settings = mean_particle.vaccination_efficacy_elderly .* dim_to_axis[:vac_protection_ratio]
    
    plots = [plot() for v in vaccination_date_settings]

    for (subplot,t_vac) in zip(plots,vaccination_date_settings)
        for (color,vaccination_distribution) in zip(colors,scheme_labels)
            data = map(d -> d.cumulative_infections_from_t_vac, analysis_data(vac_protection_ratio = :, vaccination_distribution = vaccination_distribution,ψ_0 = only(vac_labels), dynamic_threshold = only(threshold_labels) ,vaccination_begin_date = t_vac))
            no_vac_data = map(d -> d.cumulative_infections_from_t_vac, analysis_data(vac_protection_ratio = :,vaccination_distribution = "No vaccination",ψ_0 = only(vac_labels), dynamic_threshold = only(threshold_labels) ,vaccination_begin_date = t_vac))
            display(size(data[1]))
            relative_mortality =  [(n_v_d .- d)./n_v_d * 100 for (d,n_v_d) in zip(data,no_vac_data)]
            plot!(subplot,v_D_settings,mean.(relative_mortality); yerror = std.(relative_mortality),
             seriescolor = color,label = vaccination_distribution,
             yformatter = x-> string(x) * "\\%",
             ylabel = "\\% decrease in mortality",markerstrokecolor=:auto)
        end
    end
    vac_dates_formatted = map(t_vac->Dates.format(t_vac,"uuu d, yyyy"),vaccination_date_settings)
    plot!(plots[1];title = "Vaccine available: $(vac_dates_formatted[1])",legend = true)

    plot!(plots[2];title = "Vaccine available: $(vac_dates_formatted[2])",xlabel = L"v_D")
    add_subplot_letters!(plots,pos = :bottom)
    p = plot(plots...,layout = (length(vaccination_date_settings),1),size = (400,500))
    savefig("$folder/mortality_by_v_D.pdf")
end
