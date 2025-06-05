#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculates bond number for the simulation data using different models
"""
import json
from prepost.model_analysis import ModelAnalysis

if __name__ == "__main__":
        
    pressure_path = 'CFD/postProcessing/cuttingPlane/'
    velcfg_path = 'prepost/velcfg.txt'

model = ModelAnalysis(
    pressure_path=pressure_path,
    nprobes=5,
    velcfg_path=velcfg_path,
    dump2csv=False,
    plots_dir='plots/'
)

model.define_params(
    diameter=150e-6,
    rho_p=2700,
    bed_mass=0.0058,
    cg_factor=2.44
)

summary = model.model_summary()
print(summary)

# Save the summary to a JSON file
with open('pyoutputs/model_summary.json', 'w') as f:
    json.dump(summary, f, indent=4)
print(f"Model summary saved to pyoutputs/model_summary.json")

