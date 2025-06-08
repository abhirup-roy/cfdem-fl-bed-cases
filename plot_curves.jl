import Pkg

Pkg.activate(".")
Pkg.instantiate()

using FluidisedBedAnalysis

pressure_path = "CFD/postProcessing/cuttingPlane/"
velcfg_path = "prepost/velcfg.txt"
plots_dir = "plots/"

flbed = FluidisedBed(
    presure_path=pressure_path,
    n_probes=5,
    dump2csv=true,
    velcfg_path=velcfg_path,
    plots_dir="plots/"
)

plot_pressure(
    flbed,
    x_var="velocity",
    slice_dirn='z',
    use_slices=true,
)

plot_pressure(
    flbed,
    slice_dirn="z",
    x_var="time",
    png_name="pressure_time_plot_z",
    use_slices=True
)