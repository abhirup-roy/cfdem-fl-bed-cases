#usr/bin/env python3

import os
import jinja2 as jj
import numpy as np

class LIGGGHTSTemplatePopulator:

    def __init__(self, write_dir:str, template_dir:str, auto_cg:bool, **kwargs:float):
        """
        Constructor for the LIGGGHTSTemplatePopulator class:
        
        STRING write_dir:       The directory where the populated templates will be written
        STRING template_dir:    The directory where the templates are stored
        BOOL auto_cg:           Whether or not the sim inputs are to be coarse-grained
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
        self.R = kwargs.get('radius', 0.00183)
        self.density = kwargs.get('density', 1109)
        self.bed_mass = kwargs.get('bed_mass', 0.0849)
        self.contact_dumpstep = kwargs.get('contact_dumpstep', 2645)

        if auto_cg:
            if 'cg_factor' not in kwargs:
                raise ValueError('cg_factor must be provided for coarse-grained simulations. Otherwise cg must be set to False or coarsegrain factors yourself.')
            else:
                cg_factor = kwargs.get('cg_factor')
            # Apply the coarse-graining factor [Model B from Che et al. (2023)]
            R *= cg_factor
            density *= 1/cg_factor
            bed_mass *= cg_factor**3
    
    def _dump_params(self, params:list, target_dir:str):
        """
        Helper function to dump the parameters to a text file:
            LIST params: The parameters to be dumped
            STRING target_dir: The directory where the parameters will be dumped
        """
        param_arr = np.array(params)
        np.savetxt(f'{target_dir}/params.txt', param_arr)

    def _compute_workofadhesion(self, surface_energy:float, radius:float, young_mod:float, poisson_ratio:float):
        """
        Helper function to calculate the work of adhesion of a material using the JKR model (Thornton and Ning, 1998):
            FLOAT surface_energy: The surface energy of the material (in J/m^2)
            FLOAT radius: The particle radius (in meters)
            FLOAT young_mod: The Young's modulus of the material (in Pa)
            FLOAT poisson_ratio: The Poisson's ratio of the material
        """
        # Assuming monodisperse particles
        E_eq = young_mod/(1-poisson_ratio**2)
        R_eq = radius / 2
        workofadhesion = 7.09 * (surface_energy**5 * R_eq**4 / E_eq**2)**(1/3)
        return workofadhesion

    def populate_sjkr_template(self, ced:int, dump_params:bool):
        """
        Populate the template for the SJKR model:
        INT ced: The cohesion energy density value for the simulation
        BOOL dump_params: Whether or not to dump the parameters to a text file
        """
        jkr_jinja_env = jj.Environment(loader=jj.FileSystemLoader(f'{self.template_dir}/sjkr'))

        init_script_template = jkr_jinja_env.get_template('in.liggghts_init')
        init_script_contxt = {
            'CED': ced,
            'RADIUS': self.R,
            'DENSITY': self.density,
            'BED_MASS': self.bed_mass
        }
        init_script_rendered = init_script_template.render(init_script_contxt)

        run_script_template = jkr_jinja_env.get_template('in.liggghts_run')
        init_script_contxt = {
            'CED': ced,
            'DUMPSTEP_CONTACT': self.contact_dumpstep
        }
        run_script_rendered = run_script_template.render(init_script_contxt)

        with open(f'{self.write_dir}/in.liggghts_init', 'w') as f:
            f.write(init_script_rendered)
        with open(f'{self.write_dir}/in.liggghts_run', 'w') as f:
            f.write(run_script_rendered)

        if dump_params:
            params = [self.R, self.density, self.bed_mass, ced]
            self._dump_params(params=params, target_dir='postprocessing')

        

    def populate_jkr_template(self, autocomp_workofadhesion:bool, dump_params:bool, young_mod:float, 
                              poisson_ratio:float, contact_dumpstep:int, **kwargs:float):
        """
        Populate the template for the SJKR model:
        BOOL autocomp_workofadhesion: Whether or not to calculate the work of adhesion automatically for script generation
        FLOAT young_mod: The Young's modulus of the material (in Pa) - REQUIRED
        FLOAT poisson_ratio: The Poisson's ratio of the material - REQUIRED
        INT contact_dumpstep: The dumpstep for contact data
        
        KWARGS:
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
            'BED_MASS': self.bed_mass
        }
        init_script_rendered = init_script_template.render(init_script_contxt)

        run_script_template = sjkr_jinja_env.get_template('in.liggghts_run')
        init_script_contxt = {
           'DUMPSTEP_CONTACT': contact_dumpstep,
            'WORKOFADHESION': workofadhesion,
            'YOUNG_MOD': young_mod,
            'POISSON_RATIO': poisson_ratio
        }
        run_script_rendered = run_script_template.render(init_script_contxt)

        with open(f'{self.write_dir}/in.liggghts_init', 'w') as f:
            f.write(init_script_rendered)
        with open(f'{self.write_dir}/in.liggghts_run', 'w') as f:
            f.write(run_script_rendered)
        if dump_params:
            params = [self.R, young_mod, poisson_ratio, workofadhesion, self.density, self.bed_mass]
            self._dump_params(params=params, target_dir='postprocessing')





if __name__ == '__main__':
    ltp = LIGGGHTSTemplatePopulator(
        write_dir = 'DEM',
        template_dir = 'templates',
        auto_cg = False,
        radius = 0.00183,
        density = 1109,
        bed_mass = 0.0849,
        contact_dumpstep = 2645
    )

    # ltp.populate_jkr_template(autocomp_workofadhesion = True, surface_energy = 0.2e-3, young_mod = 5.4e6, poisson_ratio = 0.25, contact_dumpstep = 2645, dump_params = True)
    # ltp.populate_sjkr_template(ced = 0, dump_params = True)