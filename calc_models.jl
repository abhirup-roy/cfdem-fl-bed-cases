import Pkg

Pkg.activate(".")
Pkg.instantiate()

using FluidisedBedAnalysis


flbed = FluidisedBed(
    presure_path="CFD/postProcessing/cuttingPlane/",
    n_probes=5,
    dump2csv=false,
    velcfg_path="prepost/velcfg.txt",
    plots_dir="plots/"
)

sim_params(flbed,
    p_diameter=150e-6,
    rho_p=2700,
    cg_factor=2.44,
    poisson_ratio=0.25,
    youngs_modulus=5.4e6,
    ced=10^4)
println("Simulation parameters set")

summary = model_summary(flbed)

println("Model Summary:")
println(summary)
