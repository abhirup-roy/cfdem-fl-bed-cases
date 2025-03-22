#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Abhirup Roy
# license: MIT


"""
Plotting the fluidisation curves
"""

from prepost.plotting import ProbeAnalysis
import pandas as pd


pressure_path = 'CFD/postProcessing/cuttingPlane/'
velcfg_path = 'prepost/velcfg.txt'


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

probe_cfdem_slices.plot_pressure(slice_dirn="z", 
    x_var="time", 
    png_name="pressure_time_plot_z",
    use_slices=True
)

probe_cfdem_slices.plot_contactarea()