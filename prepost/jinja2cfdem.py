#usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Populates LIGGGHTS templating scripts for the JKR and SJKR models
"""

import os
import jinja2 as jj
import numpy as np
import json

__author__ = "Abhirup Roy"
__credits__ = ["Abhirup Roy"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Abhirup Roy"
__email__ = "axr154@bham.ac.uk"
__status__ = "Development"

class LIGGGHTSTemplatePopulator:

    def __init__(self, write_dir:str, template_dir:str, auto_cg:bool, **kwargs:float):
        """
        Constructor for the LIGGGHTSTemplatePopulator class.

        ========================================================================
        ARGS:
        ========================================================================
        STRING write_dir:       The directory where the populated templates will be written
        STRING template_dir:    The directory where the templates are stored
        BOOL auto_cg:           Whether or not the sim inputs are to be coarse-grained

        ========================================================================
        KWARGS:
        ========================================================================
        FLOAT radius:           The particle radius (in meters)
        INT density:            The density of the material (in kg/m^3)
        FLOAT bed_mass:         The mass of the bed (in kg)
        INT contact_dumpstep:   The dumpstep for contact data
        FLOAT cg_factor:        A custom coarse-graining factor (if auto_cg is set to True). Not needed if auto_cg is set to False or if you coarse-grain parameters yourself
        """

        self.write_dir = write_dir
        self.template_dir = template_dir
        # Check if the template directory exists
        if not os.path.isdir(template_dir):
            raise ValueError('Template directory does not exist.')
        
        # Note default values are for from preliminary work on Eskal150 simulations  
        self.R = kwargs.get('radius', 0.000183)
        self.density = kwargs.get('density', 1109)
        self.bed_mass = kwargs.get('bed_mass', 0.0849)
        self.contact_dumpstep = kwargs.get('contact_dumpstep', 2645)

        # Default time steps
        self.timestep_init:float = 1e-6
        self.timestep_run:float = 5e-6
        self.n_settle_steps:int = 10000

        if auto_cg:
            if 'cg_factor' not in kwargs:
                raise ValueError('cg_factor must be provided for coarse-grained simulations. Otherwise cg must be set to False or coarsegrain factors yourself.')
            else:
                cg_factor = kwargs.get('cg_factor')
            # Apply the coarse-graining factor [Model B from Che et al. (2023)]
            R *= cg_factor
            density *= 1/cg_factor
            bed_mass *= cg_factor**3
    
    def _dump_params(self, params, target_dir:str, save_as:str, filename:str):
        """
        Helper function to dump the parameters to a JSON file.

        ARGS:
            DICT/ARRAY  params: The parameters to be dumped
            STRING target_dir: The directory where the parameters will be dumped
        """
        if not os.path.isdir(target_dir):
            raise ValueError(f'Dump target directory {target_dir} does not exist.')
        if save_as == 'txt':
            np.savetxt(f'{target_dir}/{filename}.txt', params)
        
        elif save_as == 'json':
            with open(f'{target_dir}/{filename}.json', 'w') as f:
                json.dump(params, f)
        else:
            raise ValueError('Invalid file format. Must be "txt" or "json".')

    def set_timestep(self, timestep:float, kind:str):
        """
        Set the time step for the simulation for the init and/or run phase.

        ========================================================================
        ARGS:
        ========================================================================
        FLOAT timestep: The time step for the simulation
        STRING kind: The kind of time step to set ("init" or "run" or "all")
        """
        if kind == 'init':
            self.timestep_init = timestep
        elif kind == 'run':
            self.timestep_run = timestep
        elif kind == 'all':
            self.timestep_init = timestep
            self.timestep_run = timestep
        else:
            raise ValueError('Invalid time step kind. Must be "init", "run" or "all".')

    def set_init_time(self, time:int):
        """
        Set the time for the init simulation i.e. the settling time. MUST be run after set_timestep (if being used).
        ========================================================================
        ARGS:
        ========================================================================
        INT time: The settling time for the init simulation
        """
        self.n_settle_steps = int(time/self.timestep_init)


    def _compute_workofadhesion(self, surface_energy:float, radius:float, young_mod:float, poisson_ratio:float):
        """
        Helper function to calculate the work of adhesion of a material using the JKR model (Thornton and Ning, 1998):

        ARGS:
            FLOAT surface_energy: The surface energy of the material (in J/m^2)
            FLOAT radius: The particle radius (in meters)
            FLOAT young_mod: The Young's modulus of the material (in Pa)
            FLOAT poisson_ratio: The Poisson's ratio of the material
        """
        # Assuming monodisperse particles
        E_eq = 0.5 * young_mod/(1-poisson_ratio**2)
        R_eq = radius / 2
        workofadhesion = 7.09 * (surface_energy**5 * R_eq**4 / E_eq**2)**(1/3)
        return workofadhesion

    def populate_sjkr_template(self, ced:int, dump_params:bool, **kwargs):
        """
        Populate the template for the SJKR model:
        ========================================================================
        ARGS
        ========================================================================
        INT ced: The cohesion energy density value for the simulation
        BOOL dump_params: Whether or not to dump the parameters to a text file

        ========================================================================
        KWARGS
        ========================================================================
        STRING dump_filename: The name of the file to dump the parameters to
        STRING dump_filetype: The file format to dump the parameters to (txt or json)
        """

        jkr_jinja_env = jj.Environment(loader=jj.FileSystemLoader(f'{self.template_dir}/sjkr'))

        init_script_template = jkr_jinja_env.get_template('in.liggghts_init')
        init_script_contxt = {
            'CED': ced,
            'RADIUS': self.R,
            'DENSITY': self.density,
            'BED_MASS': self.bed_mass,
            'TIMESTEP_INIT': self.timestep_init,
            'SETTLETIME': self.n_settle_steps
        }
        init_script_rendered = init_script_template.render(init_script_contxt)

        run_script_template = jkr_jinja_env.get_template('in.liggghts_run')
        run_script_contxt = {
            'CED': ced,
            'DUMPSTEP_CONTACT': self.contact_dumpstep,
            'TIMESTEP_RUN': self.timestep_run
        }
        run_script_rendered = run_script_template.render(run_script_contxt)

        with open(f'{self.write_dir}/in.liggghts_init', 'w') as f:
            f.write(init_script_rendered)
        with open(f'{self.write_dir}/in.liggghts_run', 'w') as f:
            f.write(run_script_rendered)

        if dump_params:
            dump_filename:str = kwargs.get('dump_filename', 'params')
            dump_filetype:str = kwargs.get('dump_filetype', 'json')
            if dump_filetype == 'txt':
                params = [
                    f'radius={self.R}',
                    f'ced={ced}',
                    f'density={self.density}',
                    f'bed_mass={self.bed_mass}',
                    f'timestep={self.timestep_run}'
                    f'model=sJKR'
                ]
                with open(f'{self.write_dir}/{dump_filename}.txt', 'w') as f:
                    f.writelines(params)
            elif dump_filetype == 'json':
                params = {
                    "radius": self.R,
                    "ced": ced,
                    "density": self.density,
                    "bed_mass": self.bed_mass,
                    "timestep": self.timestep_run,
                    "model": "sJKR"
                }

            self._dump_params(params=params, target_dir='pyoutputs', save_as=dump_filetype, filename=dump_filename)

        

    def populate_jkr_template(self, autocomp_workofadhesion:bool, dump_params:bool, young_mod:float, 
                              poisson_ratio:float, contact_dumpstep:int, **kwargs:float):
        """
        Populate the template for the SJKR model:
        ========================================================================
        ARGS
        ========================================================================
        BOOL autocomp_workofadhesion: Whether or not to calculate the work of adhesion 
                                      automatically for script generation
        FLOAT young_mod: The Young's modulus of the material (in Pa) - REQUIRED
        FLOAT poisson_ratio: The Poisson's ratio of the material - REQUIRED
        INT contact_dumpstep: The dumpstep for contact data
        
        ========================================================================
        KWARGS
        ========================================================================
        FLOAT surface_energy: The surface energy of the material (in J/m^2) - REQUIRED if autocomp_workofadhesion is True
        FLOAT workofadhesion: The work of adhesion of the material (in J/m^2) - REQUIRED if autocomp_workofadhesion is False
        """

        sjkr_jinja_env = jj.Environment(loader=jj.FileSystemLoader(f'{self.template_dir}/jkr'))
        
        if autocomp_workofadhesion:
            if 'surface_energy' not in kwargs:
                raise ValueError('surface_energy must be provided for automatic work of adhesion calculation.')
            else:
                workofadhesion = self._compute_workofadhesion(
                    surface_energy = kwargs.get('surface_energy'),
                    radius = self.R, 
                    young_mod = young_mod,
                    poisson_ratio = poisson_ratio
                )
        else:
            workofadhesion = kwargs.get('workofadhesion', 0.116)

        init_script_template = sjkr_jinja_env.get_template('in.liggghts_init')
        init_script_contxt = {
            'RADIUS': self.R,
            'YOUNG_MOD': young_mod,
            'POISSON_RATIO': poisson_ratio,
            'WORKOFADHESION': workofadhesion,
            'DENSITY': self.density,
            'BED_MASS': self.bed_mass,
            'SETTLETIME': self.n_settle_steps,
            'TIMESTEP_INIT': self.timestep_init
        }
        init_script_rendered = init_script_template.render(init_script_contxt)

        run_script_template = sjkr_jinja_env.get_template('in.liggghts_run')
        run_script_contxt = {
           'DUMPSTEP_CONTACT': contact_dumpstep,
            'WORKOFADHESION': workofadhesion,
            'YOUNG_MOD': young_mod,
            'POISSON_RATIO': poisson_ratio,
            'TIMESTEP_RUN': self.timestep_run
        }
        run_script_rendered = run_script_template.render(run_script_contxt)

        with open(f'{self.write_dir}/in.liggghts_init', 'w') as f:
            f.write(init_script_rendered)
        with open(f'{self.write_dir}/in.liggghts_run', 'w') as f:
            f.write(run_script_rendered)
        if dump_params:
            dump_filename:str = kwargs.get('dump_filename', 'params')
            dump_filetype:str = kwargs.get('dump_filetype', 'json')
            if dump_filetype == 'txt':
                params = [
                    f'radius={self.R}',
                    f'young_mod={young_mod}',
                    f'poisson_ratio={poisson_ratio}',
                    f'workofadhesion={workofadhesion}',
                    f'density={self.density}',
                    f'bed_mass={self.bed_mass}',
                    f'timestep={self.timestep_run}'
                    f'model=JKR'
                ]
                with open(f'{self.write_dir}/{dump_filename}.txt', 'w') as f:
                    f.writelines(params)
            elif dump_filetype == 'json':
                params = {
                    "radius": self.R,
                    "young_mod": young_mod,
                    "poisson_ratio": poisson_ratio,
                    "workofadhesion": workofadhesion,
                    "density": self.density,
                    "bed_mass": self.bed_mass,
                    "timestep": self.timestep_run,
                    "model": "JKR"
                }

            self._dump_params(params=params, target_dir='pyoutputs', save_as=dump_filetype, filename=dump_filename)





if __name__ == '__main__':
    ltp = LIGGGHTSTemplatePopulator(
        write_dir = 'DEM',
        template_dir = 'templates',
        auto_cg = False,
        radius = 0.000183,
        density = 1109,
        bed_mass = 0.0849,
        contact_dumpstep = 2645
    )

    # ltp.populate_jkr_template(autocomp_workofadhesion = True, surface_energy = 0.057, young_mod = 5.4e6, poisson_ratio = 0.25, contact_dumpstep = 2645, dump_params = True)
    # ltp.populate_sjkr_template(ced = 0, dump_params = True)

