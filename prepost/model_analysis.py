#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculates bond number for the simulation data using different models
"""

from plotting import ProbeAnalysis
import pandas as pd
import numpy as np

__author__ = "Abhirup Roy"
__credits__ = ["Abhirup Roy"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Abhirup Roy"
__email__ = "axr154@bham.ac.uk"
__status__ = "Development"

class ModelAnalysis(ProbeAnalysis):
    def __init__(self, **kwargs):
        """
        Class to calculate the bond number using different models

        ========================
        Parameters
        ========================
        STRING  pressure_path: Path to the pressure data. If using slices, point to cuttingPlane directory. 
        INT     nprobes: Number of probes
        STRING  velcfg_path: Path to the vel_cfg file
        BOOL    dump2csv: Save the probe data to a csv file
        STRING  plots_dir: Directory to save the plots
        """
        
        pressure_path:str = kwargs.get("pressure_path", 'CFD/postProcessing/cuttingPlane/')
        nprobes:int = kwargs.get("nprobes", 5)
        velcfg_path:str = kwargs.get("velcfg_path", 'postprocessing/plot_P/velcfg.txt')
        dump2csv:bool = kwargs.get("dump2csv", False)
        plots_dir:str = kwargs.get("plots_dir", 'plots/')

        super().__init__(pressure_path, nprobes, velcfg_path, dump2csv=dump2csv, plots_dir=plots_dir)
        
        self._store_data()


    def _access_pressures(self):
        """
        Access the pressure data and divide into aerated and non-aerated regions 
        """

        pressure_df:pd.DataFrame = super()._probe2df(
            use_slices=True, 
            slice_dirn="z", y_agg=None
        )

        super()._calc_vel(df=pressure_df)
        vel_plot_df = pressure_df.groupby(["direction", "V_z"]).mean()

        vel_up = vel_plot_df[
                vel_plot_df.index.get_level_values(level='direction').isin(["up", "max"])
            ].reset_index('direction', drop=True).sort_index()
        
        vel_down = vel_plot_df[
                vel_plot_df.index.get_level_values(level='direction').isin(["down", "max"])
            ].reset_index('direction', drop=True).sort_index()
        
        return vel_up['Probe 0'], vel_down['Probe 0']

    def _access_voidfrac(self) -> pd.Series:
        """
        Access the void fraction data and divide into aerated and non-aerated regions
        """

        voidfrac_df = super()._read_voidfrac(slice_dirn="y", post_dir="CFD/postProcessing/cuttingPlane/")

        super()._calc_vel(df=voidfrac_df)

        vel_plot_df = voidfrac_df.groupby(["direction", "V_z"]).mean()

        vel_up = vel_plot_df[
            vel_plot_df.index.get_level_values(level='direction').isin(["up", "max"])
        ].reset_index('direction', drop=True).sort_index()
                
        vel_down = vel_plot_df[
            vel_plot_df.index.get_level_values(level='direction').isin(["down", "max"])
        ].reset_index('direction', drop=True).sort_index()

        return vel_up.squeeze(), vel_down.squeeze()
    
    def _access_contactn(self, contact_csv_path='DEM/post/collisions.csv') -> pd.Series:
        """
        Calculate the contact data and divide into aerated and non-aerated regions

        STRING  contact_csv_path: Path to the contact data
        """
        contact_df = pd.read_csv(contact_csv_path, sep="\s+", index_col="time")

        contact_df['contactn'] = contact_df.n_contact / contact_df.n_atoms 
        super()._calc_vel(df=contact_df)

        contact_plot_df = contact_df.groupby(["direction", "V_z"]).mean()
        vel_up = contact_plot_df[
            contact_plot_df.index.get_level_values(level='direction').isin(["up", "max"])
        ].reset_index('direction', drop=True).sort_index()

        vel_down = contact_plot_df[
            contact_plot_df.index.get_level_values(level='direction').isin(["down", "max"])
        ].reset_index('direction', drop=True).sort_index()

        return vel_up['contactn'], vel_down['contactn']
    
    def _store_data(self):
        """
        Store the data in the class. Called in __init__ to prevent repeated calculations
        """
        self.pressure_up, self.pressure_down = self._access_pressures()
        self.contactn_up, self.contactn_down = self._access_contactn()
        self.voidfrac_up, self.voidfrac_down = self._access_voidfrac()

        self.u_mf = self.pressure_up.idxmax()

    def define_params(self, diameter:float, rho_p:float, 
                      bed_mass:float, cg_factor:float=None):
        
        """
        Define the parameters for the model

        FLOAT   diameter: Diameter of the particle
        FLOAT   rho_p: Density of the particle
        FLOAT   bed_mass: Mass of the bed
        FLOAT   cg_factor: Coarse-graining factor. If the simulation is coarse-grained, 
                          provide the factor to scale the parameters
        """

        if cg_factor:
            self.rho_p = rho_p/cg_factor
            self.diameter = diameter*cg_factor
        else:
            self.rho_p = rho_p
            self.diameter = diameter

    def overshoot_model(self)->float:
        """
        Calculate the Bond number overshoot model (Hsu, Huang and Kuo, 2018)
        """

        if not hasattr(self, 'rho_p'):
            raise AttributeError("Define the parameters first using `define_params`")
        elif not hasattr(self, 'diameter'):
            raise AttributeError("Define the parameters first using `define_params`")

        p_1 = self.pressure_up.max()
        p_ss = self.pressure_up.iloc[-1]
        p_over = p_1-p_ss

        avg_contactn = pd.concat([self.contactn_up, self.contactn_down]).mean()
        avg_voidfrac = pd.concat([self.voidfrac_up, self.voidfrac_down]).mean()

        return (6 * p_over)/(avg_contactn**2 * (1-avg_voidfrac) * self.diameter * self.rho_p * 9.81) 

    def dhr_model(self)->float:
        """"
        Calculate the Bond number usin DHR model (Soleimani et al., 2021)"
        """
        voidfrac1 = self.voidfrac_up.loc[self.u_mf]
        voidfrac2 = self.voidfrac_down.loc[self.u_mf]

        return (voidfrac2/voidfrac1)**3 * (1-voidfrac1)/(1-voidfrac2) - 1
    
    def hyst_model(self)->float:
        """
        Calculate the Bond number using hysteresis model (Affleck et al., 2023)
        """
        p_1 = self.pressure_up.max()
        p_ss = self.pressure_up.iloc[-1]
        p_2 = self.pressure_down.loc[self.u_mf]

        delta_k = self.contactn_up.loc[self.u_mf] - self.contactn_down.loc[self.u_mf]

        return (p_1-p_2)/(p_ss*delta_k)

    def model_summary(self)->dict:
        """
        Return a summary Bond number calculated using different models
        """
        overshoot = self.overshoot_model()
        dhr = self.dhr_model()
        hyst = self.hyst_model()

        return {
            "Overshoot": overshoot,
            "DHR": dhr,
            "Hysteresis": hyst
        }

if __name__ == "__main__":
    pass
    # See post_scripts/model_analysis.py for usage