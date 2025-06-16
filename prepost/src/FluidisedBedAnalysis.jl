module FluidisedBedAnalysis

# Based module for fluidised bed pre/post processing

include("curve_plots.jl")
include("model_analysis.jl")

using .CurvePlots
using .BoModels

export FluidisedBed,
       plot_pressure,
       plot_voidfrac,
       sim_params,
       overshoot_model,
       dhr_model,
       hyst_model,
       intrinsic_bond_num,
       model_summary


end