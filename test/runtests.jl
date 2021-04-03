using COVID_model
using Test
using Dates
using AxisKeys
using DataStructures
using Serialization
using StatsBase
#testing could use work, writing good tests for research code is hard!
file_list = ["baseline",]
@test create_model_output(;num_particles = 2 ,file_list = file_list)

@testset "scenario $folder" for folder in file_list
    @testset "output $fname" for fname in readdir(joinpath(COVID_model.PACKAGE_FOLDER,"output",folder))
        if occursin("fitting",fname)
            test_data = deserialize(joinpath(COVID_model.PACKAGE_FOLDER,"test/output_test",folder,fname)) 
            data = deserialize(joinpath(COVID_model.PACKAGE_FOLDER,"output",folder,fname))
            display(fieldnames(typeof(data)))
            for field in fieldnames(typeof(data))
                @test all(getfield(data,field) .≈ getfield(test_data,field))
            end
        else
            test_data = deserialize(joinpath(COVID_model.PACKAGE_FOLDER,"test/output_test",folder,fname)) |> AxisKeys.unname |> AxisKeys.keyless
            data = deserialize(joinpath(COVID_model.PACKAGE_FOLDER,"output",folder,fname))|> AxisKeys.unname |> AxisKeys.keyless
            for index in eachindex(test_data)
                for field in fieldnames(typeof(data[index]))
                    @test all(getfield(data[index],field) .≈ getfield(test_data[index],field))
                end
            end
        end
    end
end


# ontario_covid_data = COVID_model.get_ontario_data()
# default_parameters = COVID_model.get_default_model_parameters(ontario_covid_data.start_date,  ontario_covid_data.infection_data_incident[1,:])

# init_params = COVID_model.apply_workplace_closure_fit(default_parameters,ontario_covid_data)

# @testset "workplace fit $v" for v in (:work_opening_rate,:work_closure_rate,:ε_W)
#     @test abs(init_params[v] - mean_particle[v])/mean_particle[v] < 0.05
# end;
# opt_params,_ = COVID_model.get_opt_params()
# baseline_testing,mean_particle_testing,converged_testing = COVID_model.bayesian_model_selection(init_params,500,7.0,opt_params,ontario_covid_data; test =true)
# display(mean_particle_testing)
# display(mean_particle)
# @test converged_testing

# for v in keys(opt_params)
#     println("$v save: $(mean_particle[v]), testing: $(mean_particle_testing[v])")
# end
# @testset "parameter fitting" begin
#     @testset "parameter fit $v" for v in keys(opt_params)
#         particle_list = []
#         for particle in baseline
#             push!(particle_list,particle[v])
#         end
#         measured_std = std(particle_list)
#         println("std: $measured_std")
#         @test abs(mean_particle_testing[v] - mean_particle[v]) < measured_std
#     end
# end