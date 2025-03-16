#!/usr/bin/env python3

from postprocessing.bondno.model_analysis import ModelAnalysis
from postprocessing.plot_P.probe_analysis import ProbeAnalysis

"""
COMPARING MODEL RESULTS
"""

model = ModelAnalysis(
    pressure_path='CFD/postProcessing/cuttingPlane/',
    nprobes=5,
    velcfg_path='postprocessing/plot_P/velcfg.txt',
    dump2csv=False,
    plots_dir='plots/'
)

model.define_params(
    diameter=150e-6,
    rho_p=2700,
    cg_factor=2.44
)

print(model.model_summary())


"""
PLOTTING PRESSURE AND VOID FRACTION VS VELOCITY
"""

pressure_path = '../../CFD/postProcessing/cuttingPlane/'
velcfg_path = 'velcfg.txt'


probe_cfdem_slices = ProbeAnalysis(
    pressure_path=pressure_path,
    nprobes=5,
    velcfg_path=velcfg_path,
    dump2csv=False
)

probe_cfdem_slices.plot_pressure(slice_dirn="z", 
    x_var="velocity", 
    png_name="pressure_vel_plot_z",
    use_slices=True
)

probe_cfdem_slices.plot_voidfrac(
    slice_dirn="y", 
    x_var="velocity",
    png_name="voidfrac_time_plot_y"
)