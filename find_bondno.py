#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculates bond number for the simulation data using different models
"""

from prepost.model_analysis import ModelAnalysis

pressure_path = 'CFD/postProcessing/cuttingPlane/'
velcfg_path = 'prepost/velcfg.txt'

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
    bed_mass=0.0058,
    cg_factor=2.44
)

print(model.model_summary())