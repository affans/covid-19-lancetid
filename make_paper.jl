using Pkg
Pkg.instantiate()

using COVID_model
using Plots
pgfplotsx()
default(dpi = 300)
default(framestyle = :box)
create_model_output()
do_model_plots()