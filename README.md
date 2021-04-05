https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(21)00057-8/fulltext#seccestitle80

This repo contains the code for the paper "Prioritising COVID-19 vaccination in changing social and epidemiological landscapes"
by Peter C. Jentsch, Madhur Anand, Chris T Bauch


Provided you have Julia 1.5.x or above, generating the paper figures is (in theory) as simple as
```
git clone https://git.uwaterloo.ca/pjentsch/interactions-between-social-distancing-and-vaccination.git
cd interactions-between-social-distancing-and-vaccination
julia --project make_paper.jl
```

and then waiting 1-10 days, depending on how big your computer is.
It will use however many threads you give Julia.

You can ensure all tests pass by running
```
]test COVID_model
```
in the Julia REPL after you have added the package.

This package has only been tested on linux, although it should work on anything else with full support for Julia 1.5.x.

The package exports three functions, they are:


- `create_simulation_scenarios()`: does parameter fitting, saves to files in /parameter_fits/
- `create_model_output()`: computes all the output needed for the plotting, saves to /output/
- `do_model_plots()`: makes the plots from the data in /output/

The script `fetch_data.sh` should fetch the most recent COVID-19 and distancing data for Ontario, Canada, and save it to the right place, but google and Public Health Ontario break it pretty often. The repo currently contains the data used in whatever the most recent version of the paper is.

# Summary of Code Structure

There are four steps in the model:
- loading/parsing cases and mobility data located in /data/
- fitting model parameters to cases and mobility data (begins from `create_simulation_scenarios()`)
- running model over several multidimensional arrays of parameter values, generated in `scenario_data.jl`, this process begins from `create_model_output()`
- plotting model output, begins from `do_model_plots()`, contained in `plots.jl`


## Model Structure

The individual infection model is a compartmental ODE model (see paper for more details), defined in `model.jl`. 

Parameters are stored and passed to the ODE solver in the form of NamedTuples. The default parameter NamedTuple can be obtained by calling `get_default_model_parameters()`.

The model allocation and computing the model solution are split into two steps, given by `allocate_realization(...)` and `run_model(...)`, or `run_model!(...)`.

We also preallocate the space used for the model solution, in the form of three types, `ModelFitData`, `HeatmapData`, and `TimeseriesData`. `run_model!(...)` calls the method `build_model_output!` which dispatches on these types, and given a set of model solutions fills the type with the appropriate data from them. 

The method `run_model(...)` is used during the model fitting and only generates solutions from a single model realization, as opposed to a list of them. 

In hindsight, I think there are better ways to solve these problems, but this was my first real project in Julia.

There is certainly a lot more documentation to be written, which I will work on as I have the time, but definitely email me if you have questions about the code. Even if you don't have questions, if you find this code useful let me know. 

Peter
