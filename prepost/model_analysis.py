#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculates bond number for the simulation data using different models
"""

from .plotting import FlBedPlot
import pandas as pd
import numpy as np
from typing import Optional

__author__ = "Abhirup Roy"
__credits__ = ["Abhirup Roy"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Abhirup Roy"
__email__ = "axr154@bham.ac.uk"
__status__ = "Development"


class ModelAnalysis(FlBedPlot):
    def __init__(self, **kwargs):
        """
        Initialise object to calculate the bond number using different models

        Args:
          pressure_path:
            Path to the pressure data. If using slices, point to CFD cuttingPlane directory.
          probes:
            Number of probes in simulation. Default is 5
          velcfg_path:
            Path to the vel_cfg file
          dump2csv:
            Whether or not to save the probe data to a csv file
          plots_dir:
            Directory to save the plots
        """

        pressure_path: str = kwargs.get(
            "pressure_path", "CFD/postProcessing/cuttingPlane/"
        )
        nprobes: int = kwargs.get("nprobes", 5)
        velcfg_path: str = kwargs.get("velcfg_path", "prepost/velcfg.txt")
        dump2csv: bool = kwargs.get("dump2csv", False)
        plots_dir: str = kwargs.get("plots_dir", "plots/")

        super().__init__(
            pressure_path=pressure_path,
            nprobes=nprobes,
            velcfg_path=velcfg_path,
            dump2csv=dump2csv,
            plots_dir=plots_dir,
        )

        self._store_data()

    def _access_pressures(self) -> tuple[pd.Series, pd.Series]:
        """
        Helper function to load the pressure data and divide into aerated and non-aerated regions

        Returns:
            A tuple of pd.Series of pressure data with increasing and decreasing velocity
            for probe 0
        """

        pressure_df: pd.DataFrame = super()._probe2df(
            use_slices=True, slice_dirn="z", y_agg=None
        )

        super()._calc_vel(df=pressure_df)
        vel_plot_df = pressure_df.groupby(["direction", "V_z"]).mean()

        vel_up = (
            vel_plot_df[
                vel_plot_df.index.get_level_values(level="direction").isin(
                    ["up", "max"]
                )
            ]
            .reset_index("direction", drop=True)
            .sort_index()
        )

        vel_down = (
            vel_plot_df[
                vel_plot_df.index.get_level_values(level="direction").isin(
                    ["down", "max"]
                )
            ]
            .reset_index("direction", drop=True)
            .sort_index()
        )

        return vel_up["Probe 0"], vel_down["Probe 0"]

    def _access_voidfrac(self) -> tuple[pd.Series, pd.Series]:
        """
        Helper function to load void fraction data and divide into
        increasing and decreasing velocity regions

        Returns:
          A tuple of pd.Series of void fraction data with increasing and decreasing velocity
        """

        voidfrac_df = super()._read_voidfrac(
            slice_dirn="y", post_dir="CFD/postProcessing/cuttingPlane/"
        )

        super()._calc_vel(df=voidfrac_df)

        vel_plot_df = voidfrac_df.groupby(["direction", "V_z"]).mean()

        vel_up = (
            vel_plot_df[
                vel_plot_df.index.get_level_values(level="direction").isin(
                    ["up", "max"]
                )
            ]
            .reset_index("direction", drop=True)
            .sort_index()
        )

        vel_down = (
            vel_plot_df[
                vel_plot_df.index.get_level_values(level="direction").isin(
                    ["down", "max"]
                )
            ]
            .reset_index("direction", drop=True)
            .sort_index()
        )

        squeezed_up = vel_up.squeeze()
        squeezed_down = vel_down.squeeze()

        # Ensure we return Series objects
        if not isinstance(squeezed_up, pd.Series):
            raise TypeError("squeezed_up is not a pd.Series")
        if not isinstance(squeezed_down, pd.Series):
            raise TypeError("squeezed_down is not a pd.Series")

        return squeezed_up, squeezed_down

    def _access_contactn(
        self, contact_csv_path="DEM/post/collisions.csv"
    ) -> tuple[pd.Series, pd.Series]:
        """
        Helper function to read the contact data and divide into increasing and decreasing
        velocity regions

        Args:
          contact_csv_path: Path to the contact data csv file. Default is 'DEM/post/collisions.csv'

        Returns:
          A tuple of pd.Series of contact number data with increasing and decreasing velocity
        """
        contact_df = super()._read_collisions(contact_csv_path, calltype="contactn")
        contact_df.set_index("time", inplace=True)
        contact_df.index -= contact_df.index.min()

        super()._calc_vel(df=contact_df)

        contact_plot_df = contact_df.groupby(["direction", "V_z"]).mean()
        vel_up = (
            contact_plot_df[
                contact_plot_df.index.get_level_values(level="direction").isin(
                    ["up", "max"]
                )
            ]
            .reset_index("direction", drop=True)
            .sort_index()
        )

        vel_down = (
            contact_plot_df[
                contact_plot_df.index.get_level_values(level="direction").isin(
                    ["down", "max"]
                )
            ]
            .reset_index("direction", drop=True)
            .sort_index()
        )

        return vel_up["contactn"], vel_down["contactn"]

    def _store_data(self):
        """
        Store the data in the class. Called in __init__ to prevent repeated calculations
        """
        self.pressure_up, self.pressure_down = self._access_pressures()
        self.contactn_up, self.contactn_down = self._access_contactn()
        self.voidfrac_up, self.voidfrac_down = self._access_voidfrac()

        self.u_mf = self.pressure_up.idxmax()

    def define_params(
        self,
        diameter: float,
        rho_p: float,
        bed_mass: float,
        cg_factor: Optional[float] = None,
    ):
        """
        Define the parameters for the model. Must be called before calculating the Bond number.

        Args:
          diameter:
            Diameter of the particles (in m)
          rho_p:
            Density of the particles (in kg/m^3)
          bed_mass:
            Mass of the bed (in kg)
          cg_factor:
            Coarse-graining factor. If not provided, no coarse-graining is applied.
        """

        if cg_factor:
            self.rho_p = rho_p / cg_factor
            self.diameter = diameter * cg_factor
        else:
            self.rho_p = rho_p
            self.diameter = diameter

    def overshoot_model(self) -> float:
        """
        Calculate the Bond number overshoot model from Hsu, Huang and Kuo (2018)
        """

        if not hasattr(self, "rho_p"):
            raise AttributeError("Define the parameters first using `define_params`")
        elif not hasattr(self, "diameter"):
            raise AttributeError("Define the parameters first using `define_params`")

        p_1 = self.pressure_up.max()
        p_ss = self.pressure_up.iloc[-1]
        p_over = p_1 - p_ss

        avg_contactn = pd.concat([self.contactn_up, self.contactn_down]).mean()
        avg_voidfrac = pd.concat([self.voidfrac_up, self.voidfrac_down]).mean()

        return (6 * p_over) / (
            avg_contactn**2 * (1 - avg_voidfrac) * self.diameter * self.rho_p * 9.81
        )

    def dhr_model(self) -> float:
        """ "
        Calculate the Bond number usin DHR model from Soleimani et al. (2021)"
        """
        voidfrac1 = self.voidfrac_up.loc[self.u_mf]
        voidfrac2 = self.voidfrac_down.loc[self.u_mf]

        return (voidfrac2 / voidfrac1) ** 3 * (1 - voidfrac1) / (1 - voidfrac2) - 1

    def hyst_model(self) -> float:
        """
        Calculate the Bond number using hysteresis model from Affleck et al. (2023)
        """
        p_1 = self.pressure_up.max()
        p_ss = self.pressure_up.iloc[-1]
        p_2 = self.pressure_down.loc[self.u_mf]

        delta_k = np.abs(
            self.contactn_up.loc[self.u_mf] - self.contactn_down.loc[self.u_mf]
        )
        return (p_1 - p_2) / (p_ss * delta_k)

    def model_summary(self) -> dict:
        """
        Return a summary Bond number calculated using different models
        """
        overshoot = self.overshoot_model()
        dhr = self.dhr_model()
        hyst = self.hyst_model()

        return {"Overshoot": overshoot, "DHR": dhr, "Hysteresis": hyst}


if __name__ == "__main__":
    pass
    # See model_analysis.py for usage
