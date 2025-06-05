#usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Abhirup Roy
# license: MIT

"""
Templating for CFDEM simulations
"""

from prepost.jinja2cfdem import LIGGGHTSTemplatePopulator

if __name__ == "__main__":
    write_dir = 'DEM'
    template_dir = 'templates'

    ltp = LIGGGHTSTemplatePopulator(
        write_dir = write_dir,
        template_dir = template_dir,
        auto_cg = False,
        radius = 0.000183,
        density = 1109,
        bed_mass = 0.0849,
        contact_dumpstep = 2645
    )

    # ltp.populate_jkr_template(autocomp_workofadhesion = True, surface_energy = 0.057, young_mod = 5.4e6, poisson_ratio = 0.25, contact_dumpstep = 2645, dump_params = True)
    ltp.populate_sjkr_template(ced = 1e2, dump_params = True)