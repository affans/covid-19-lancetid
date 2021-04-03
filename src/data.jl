export get_ontario_data

struct OntarioData{A,B,C}
    data_length::Int
    start_date::A
    infection_data_cumulative::B 
    workplace_closure::C
    distancing_data::C
    infection_data_incident::B
    infection_data_incident_raw::B
    cases_from_travel::B
    mortality_over_time::B
end

function get_distancing_data_sets()::Tuple{Vector{Float64},Vector{Float64}}
    f = readdlm(joinpath(PACKAGE_FOLDER,"data/csv/distancing_data.csv"), ',')
    df = filter(row -> row[:country_region_code] == "CA" && row[:sub_region_1] == "Ontario", DataFrame([f[1,i] => f[2:end, i] for i = 1:length(f[1,:])]))
    distancing_data = clamp.((-1 .* filter(x -> typeof(x) <: Real, df[:,:retail_and_recreation_percent_change_from_baseline])) / 100, 0.0, 1.0)
    workplace_data = clamp.((-1 .* filter(x -> typeof(x) <: Real, df[:,:workplaces_percent_change_from_baseline])) / 100, 0.0, 1.0)
    return (workplace_data, distancing_data)
end
function get_canada_demographic_distribution()::Tuple{Vector{Tuple{Float64,Float64}},Vector{Float64}}
    f = readdlm(joinpath(PACKAGE_FOLDER,"data/csv/demographic_data.csv"), ',')
    df = DataFrame([f[1,i] => f[2:end, i] for i = 1:length(f[1,:])])
    return map(parse_string_as_float_pair, df[:,:demographic_data_bins]), df[:,:demographic_data]
end
function get_canada_case_fatality()::Tuple{Vector{Tuple{Float64,Float64}},Vector{Float64}}
    f = readdlm(joinpath(PACKAGE_FOLDER,"data/csv/case_fatality_data.csv"), ',')
    df = DataFrame([f[1,i] => f[2:end, i] for i = 1:length(f[1,:])])
    return  map(parse_string_as_float_pair, df[:,:case_fatality_bins]), df[:,:case_fatality_data]
    # https://www.publichealthontario.ca/-/media/documents/ncov/epi/covid-19-severe-outcomes-ontario-epi-summary.pdf?la=en
end
function parse_string_as_float_pair(s)::Tuple{Float64,Float64}
    parsed = strip(s, ['(',')']) |> s -> split(s, ",")
    if length(parsed) > 2
        error("input tuple too long")
    end
    return (parse(Float64, parsed[begin]), parse(Float64, parsed[end]))
end

function parse_cases_data()::NamedTuple{(:incident_cases_raw,:incident_cases, :cumulative_cases, :cases_from_travel, :mortality_over_time),NTuple{5,Array{Float64,2}}}
    f = readdlm(joinpath(PACKAGE_FOLDER,"data/csv/COVID_ontario_data.csv"), ',')
    start_date = Date("2020-01-01")
    last_date = Date(2020, 1, 1)
    for line in eachrow(f[2:end,:])
        if Date(line[2]) > last_date 
            last_date = Date(line[2])
        end
    end
    end_date = last_date - Day(7) # date to which data is likely to change according to Ontario government
    bins = [(0.0, 20.0),(20.0, 30.0),(30.0, 40.0),(40.0, 50.0),(50.0, 60.0),(60.0, 70.0),(70.0, 80.0),(80.0, 90.0),(90.0, 125.0)]
    infection_data = zeros(Float64,(end_date - start_date |> Dates.value) + 1, length(bins))
    deaths = zeros(Float64,(end_date - start_date |> Dates.value) + 1, length(bins))
    cases_from_travel = zeros(Float64,(end_date - start_date |> Dates.value) + 1, length(bins))
    
    
    for (i,line) in enumerate(eachrow(f[2:end,:]))
        date = Date(line[2])
        if start_date <= date < end_date
            bin = line[6]
            if length(bin)>0 && bin != "UNKNOWN"
                if '<' in bin
                    cleaned_bin = 1
                else
                    cleaned_bin = parse(Int, strip(bin, ['0','s']))
                end
                if line[9] == "Fatal"
                    deaths[(date - start_date |> Dates.value) + 1 , cleaned_bin] +=  1.0
                end

                if occursin("travel", lowercase(line[8]))
                    cases_from_travel[(date - start_date |> Dates.value) + 1 , cleaned_bin] +=  1.0
                else
                    infection_data[(date - start_date |> Dates.value) + 1 , cleaned_bin] += 1.0
                end
            else #if age is not specified, assume age bin follows distribution up until then
                if line[9] == "Fatal"
                    deaths[(date - start_date |> Dates.value) + 1,:] .+=  sum.(eachcol(deaths))./i
                end

                if occursin("travel", lowercase(line[8]))
                    cases_from_travel[(date - start_date |> Dates.value) + 1,:] .+=  sum.(eachcol(cases_from_travel))./i
                else
                    infection_data[(date - start_date |> Dates.value) + 1,:] .+= sum.(eachcol(infection_data))./i
                end
            end
        end
    end
    moving_total(arr,len) = [sum(arr[max(1,i-len):i,j]) for i in 1:length(arr[:,1]), j in 1:length(bins)] 
    moving_average(arr,len) = [sum(arr[max(1,i-len):i,j])/len for i in 1:length(arr[:,1]), j in 1:length(bins)]
  
    mortality_over_time::Array{Float64,2} =  moving_total(deaths,30) ./ moving_total(infection_data.+cases_from_travel,30)
    incident_cases = moving_average(infection_data,7)

    return (
        incident_cases_raw = infection_data,
        incident_cases = incident_cases,
        cumulative_cases = [sum(incident_cases[1:i,j]) for i in 1:length(incident_cases[:,1]), j in 1:length(bins)],
        cases_from_travel = moving_average(cases_from_travel,7),
        mortality_over_time = mortality_over_time
    )
end
function get_default_model_parameters()
    infection_start = ontario_data.start_date
    I_0 = ontario_data.infection_data_incident[begin,:]
    parameters = (
        κ_vac = 0.0,
        c_vac = 0.0,
        p_ul = 0.01,
        σ_0 = 1/2,
        σ_1 = 1/2,
        γ_s = 1/7,
        γ_a = 1.7,
        κ_1 = 1.1e4,
        κ_2 = 10.0,
        R_0 = 2.0,
        P_multiplier = 3.0,
        E_multiplier = 0.5,
        r_modifier = 0.0,
        c_1 = 5e-5,
        c_2 = 5e-5,
        ε_P_1 = 0.3,
        ε_W =0.5,
        s = -0.4,
        a_m = 0.05,
        ϕ = -80.0,
        I_multiplier = 2.0,
        t_switch = 160,
        ψ_0 = 0.0,
        x_0 = 0.01,
        x_vac_0 = 1.0,
        η = 0.15,
        vac_protection_ratio = 1.0,
        work_closure_rate = 0.1,
        work_opening_rate = 0.1,
        dynamic_threshold = population - 1,
        infection_start = infection_start,
        vaccination_begin_date = Date(3000, 1, 1),
        vaccination_efficacy_elderly = 0.75,
        vaccination_efficacy_else = 0.75, 
        r_1 = 1.0,
        r_2 = 1.0,
        r_3 = 1.0,
        a_1_1 = 0.5203184244350092, #start with good parameters so that successive runs of
        a_1_2 = 0.3316526579948613, #parameter fitting are easier
        a_1_3 = 0.7993083727994479, 
        a_2_1 = 0.9344602582238042,
        a_2_2 = 0.2655124202811558,
        a_2_3 = 0.1134815835693088,
        I_0 = I_0
    )
    return parameters
end
function load_contact_matrices()::(NamedTuple{Sym,NTuple{5,T}} where {T<:AbstractArray,Sym})
    data_files = (load(joinpath(PACKAGE_FOLDER,"data/contact_matrices/contact_others.rdata"))["contact_others"],
        load(joinpath(PACKAGE_FOLDER,"data/contact_matrices/contact_school.rdata"))["contact_school"],
        load(joinpath(PACKAGE_FOLDER,"data/contact_matrices/contact_work.rdata"))["contact_work"],
        load(joinpath(PACKAGE_FOLDER,"data/contact_matrices/contact_home.rdata"))["contact_home"])
    contact_matrices = ([a["CAN"] for a in data_files])
    
    return (
        other_contacts  = contact_matrices[1],
        school_contacts = contact_matrices[2],
        work_contacts = contact_matrices[3],
        home_contacts = contact_matrices[4],
        contact_matrices_sum = contact_matrices[1] .+ contact_matrices[2] .+ contact_matrices[3] .+ contact_matrices[4]
    )
end

function get_ontario_data()
    infected_ts_bins = [(0.0, 20.0),(20.0, 30.0),(30.0, 40.0),(40.0, 50.0),(50.0, 60.0),(60.0, 70.0),(70.0, 80.0),(80.0, 90.0),(90.0, 125.0)]
    infected_ts_raw,infected_ts, infected_ts_cumulative, cases_from_travel, mortality_data = parse_cases_data()
    distancing_data_start_date = Date(2020, 02, 15)
    distancing_data = [vcat(zeros(Dates.value(distancing_data_start_date - Date(2020, 1, 1))), l) for l in get_distancing_data_sets()]
    infected_ts_bins = [(0.0, 20.0),(20.0, 30.0),(30.0, 40.0),(40.0, 50.0),(50.0, 60.0),(60.0, 70.0),(70.0, 80.0),(80.0, 90.0),(90.0, 125.0)]
   
    jan1 = Date(2020, 1, 1)
    infected_ts_binned_cumulative = mapslices(l -> rebin_data(infected_ts_bins, l, age_bins), infected_ts_cumulative, dims=[2])
    infected_ts_binned = mapslices(l -> rebin_data(infected_ts_bins, l, age_bins), infected_ts, dims=[2]) 
    infected_ts_raw_binned = mapslices(l -> rebin_data(infected_ts_bins, l, age_bins), infected_ts_raw, dims=[2]) 
    
    # infected_ts_binned_smoothed = mapslices(l->moving_average(l),infected_ts_binned,dims = [1])
    infection_start = jan1 + Day(findfirst(x -> x >= 50.0, sum.(eachrow(infected_ts_binned))))
    infection_start_value = Dates.value(infection_start - jan1)

    cases_from_travel_binned = mapslices(l -> rebin_data(infected_ts_bins, l, age_bins), cases_from_travel, dims=[2])
    mortality_over_time_binned = mapslices(l -> rebin_data(infected_ts_bins, l, age_bins), mortality_data, dims=[2])

    return OntarioData(
        length(infected_ts_binned_cumulative[infection_start_value:end,1]),
        infection_start,
        infected_ts_binned_cumulative[infection_start_value:end,:],
        distancing_data[1][infection_start_value:end],
        distancing_data[2][infection_start_value:end],
        infected_ts_binned[infection_start_value:end,:],
        infected_ts_raw_binned[infection_start_value:end,:],
        cases_from_travel_binned[infection_start_value:end,:],
        mortality_over_time_binned[infection_start_value:end,:]
    )
end
