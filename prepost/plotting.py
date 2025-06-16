#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plotting utils for CFDEM fluidised bed simulations.
"""

import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
import pyvista as pv
from pyevtk.hl import pointsToVTK
from warnings import warn

__author__ = "Abhirup Roy"
__credits__ = ["Abhirup Roy"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Abhirup Roy"
__email__ = "axr154@bham.ac.uk"
__status__ = "Development"


class FlBedPlot():

    def __init__(
            self,
            pressure_path: str,
            nprobes: int,
            velcfg_path: str,
            dump2csv: bool = True,
            plots_dir: str = 'plots/'):
        """
        Initialise the FlBedPlot class for plotting pressure and void fraction data from fluidised bed simulations

        Args:
          pressure_path:
            Path to the pressure data. If using slices, point to cuttingPlane directory.
          nprobes:
            Number of probes in simulation
          velcfg_path:
            Path to the vel_cfg file
          dump2csv:
            Save the probe data to a csv file
          plots_dir:
            Directory to save the plots

        Raises:
          FileNotFoundError: If the pressure_path, velcfg_path or plots_dir does not exist
        """
        # Check if the paths exist
        if not os.path.exists(pressure_path):
            raise FileNotFoundError(
                f"Pressure data at {pressure_path} does not exist")

        if not os.path.exists(velcfg_path):
            raise FileNotFoundError(
                f"Velocity config at {velcfg_path} does not exist")
        if not os.path.isdir(plots_dir):
            raise FileNotFoundError(
                f" Plots directory at {plots_dir} does not exist")

        self.pressure_path = pressure_path
        self.nprobes = nprobes
        self.dump2csv = dump2csv
        self.velcfg_path = velcfg_path
        self.plots_dir = plots_dir

        rcParams.update({'font.size': 20})

    def _probe2df(
            self,
            use_slices: bool,
            slice_dirn: str,
            y_agg: str) -> pd.DataFrame:
        """
        Helper function to map the probe data to a pandas dataframe, with data indexed by time.
        Dumps the data to a csv file if dump2csv attribute of class is True.

        Args:
          use_slices:
            Whether to use slices or probes. If True, reads the pressure data from the cuttingPlane directory.
          slice_dirn:
            Direction of the slice. "z" for z-normal slices, "y" for y-normal slices.
          y_agg:
            Aggregation method for y-normal slices. "cdf_median" for median of CDF, "mean" for mean, "median" for median.

        Returns:
          A pandas DataFrame with the pressure data indexed by time.

        Raises:
            ValueError: If use_slices is True and slice_dirn is not "z" or "y".
            ValueError: If y_agg is not "cdf_median", "mean" or "median".
        """
        if use_slices:
            times = os.listdir(self.pressure_path)
            pressure_dict = dict()

            for dir in times:
                if slice_dirn == "z":
                    p_lst = []
                    pdata = pv.read(
                        self.pressure_path + '/' + dir + '/p_zNormal0.vtk')
                    p_arr = pdata.get_array('p')
                    p_lst.append(p_arr.mean().item())

                    pdata = pv.read(
                        self.pressure_path + '/' + dir + '/p_zNormal1.vtk')
                    p_arr = pdata.get_array('p')
                    p_lst.append(p_arr.mean().item())

                    pdata = pv.read(
                        self.pressure_path + '/' + dir + '/p_zNormal2.vtk')
                    p_arr = pdata.get_array('p')
                    p_lst.append(p_arr.mean().item())

                    pdata = pv.read(
                        self.pressure_path + '/' + dir + '/p_zNormal3.vtk')
                    p_arr = pdata.get_array('p')
                    p_lst.append(p_arr.mean().item())

                    pdata = pv.read(
                        self.pressure_path + '/' + dir + '/p_zNormal4.vtk')
                    p_arr = pdata.get_array('p')
                    p_lst.append(p_arr.mean().item())

                    pressure_dict[float(dir)] = p_lst

                    pressure_df = pd.DataFrame.from_dict(
                        pressure_dict, orient='index', columns=[
                            f"Probe {i}" for i in range(
                                self.nprobes)]).sort_index()

                elif slice_dirn == "y":
                    pdata = pv.read(
                        self.pressure_path + '/' + dir + '/p_yNormal.vtk')
                    p_arr = pdata.get_array('p')

                    if y_agg == 'cdf_median':
                        pressure_dict[float(dir)] = p_arr
                        pressure_dict = {
                            k: self.find_cdfmedian(v) for k,
                            v in pressure_dict.items()}
                    elif y_agg == 'mean':
                        pressure_dict[float(dir)] = p_arr.mean().item()
                    elif y_agg == "median":
                        pressure_dict[float(dir)] = np.median(p_arr)
                    else:
                        raise ValueError(
                            "Invalid aggregation method. Choose 'cdf_median' or 'mean'")
                    pressure_df = pd.DataFrame.from_dict(
                        pressure_dict, orient='index', columns=['pressure']).sort_index()
                    pressure_df.to_csv("probe_pressure.csv")

                else:
                    raise ValueError(
                        "Invalid slice direction. Choose 'z' or 'y'")

        else:
            # Make df from the probe data
            headers = ["Probe Time"]
            for i in range(self.nprobes):
                headers.append(f"Probe {i}")

            pressure_df = pd.read_csv(
                self.pressure_path,
                delim_whitespace=True,
                comment='#',
                names=headers,
                header=None,
            ).set_index("Probe Time")

        # Dump to csv if specified
        if self.dump2csv:
            pressure_df.to_csv("probe_pressure.csv")

        return pressure_df

    def _read_probetxt(self):
        """
        Helper function to read the velocity config file and extract the time and corresponding fluid velocities.
        """
        with open(self.velcfg_path, "r") as f:
            probe_text = f.read().splitlines(False)

            # read plot times and corresponding fluid velocities
            self.t = []
            self.v_z = []

            for line in probe_text:
                line_splt = line.replace("(", "").replace(")", "").split()
                self.t.append(float(line_splt[0]))
                self.v_z.append(float(line_splt[-1]))
            print("Selected times: ", self.t)
            print("Corresponding vel: ", self.v_z)

    def _calc_vel(self, df: pd.DataFrame) -> None:
        """
        Helper function to map the velocity to pressure. Reads the velocity config file and maps the velocity to the time-series data.
        This function assumes that the DataFrame is indexed by time and contains pressure data. It adds a new column "V_z" to the DataFrame
        with the corresponding fluid velocities and a "direction" column indicating whether the velocity is increasing, decreasing or at max.

        Args:
          df:
            A pandas DataFrame with the pressure data indexed by time.
        """
        # Initialise bounds and velocity
        bounds = []
        vel = []

        self._read_probetxt()
        # Find the bounds for velocity
        for i in range(len(self.t) - 1):
            if self.v_z[i] == self.v_z[i + 1]:
                bounds.append([self.t[i], self.t[i + 1]])
                vel.append(self.v_z[i])

                if self.v_z[i] == max(self.v_z):
                    max_vel_t1, max_vel_t2 = self.t[i], self.t[i + 1]

            else:
                pass

        # Lower and upper time bounds for each velocity
        lb = [b[0] for b in bounds]
        ub = [b[1] for b in bounds]
        vz_arr = np.zeros_like(
            df.index.to_numpy(),
            dtype=float)  # Initialize as float array

        # Map the velocity to the pressure data
        for i in range(len(bounds)):
            mask = (
                df.index.to_numpy() > lb[i]) & (
                df.index.to_numpy() < ub[i])
            vz_arr[mask] = vel[i]

            if i < len(bounds) - 1:
                gap_mask = (df.index.to_numpy() >= ub[i]) & (
                    df.index.to_numpy() <= lb[i + 1])
                vz_arr[gap_mask] = np.nan
        df["V_z"] = vz_arr

        def _map_direction(x) -> None:
            """
            Simple helper function to map the direction of the velocity. i.e whether it is increasing, decreasing or at max
            """
            if x < max_vel_t1:
                return "up"
            elif x >= max_vel_t1 and x <= max_vel_t2:
                return "max"
            else:
                return "down"
        # Add direction to the dataframe
        df["direction"] = df.index.to_series().apply(_map_direction)

    def find_cdfmedian(self, arr: np.array) -> float:
        """
        Finds the median of the CDF of the given array.

        Args:
          arr:
            A numpy array of values for which the median of the CDF is to be calculated.

        Returns:
            The median of the CDF of the given array as a float.
        """

        x, counts = np.unique(arr, return_counts=True)
        cusum = np.cumsum(counts)
        cdf = cusum / cusum[-1]

        median_idx = cdf.tolist().index(np.percentile(cdf, 50, method='nearest'))
        return x[median_idx].item()

    def plot_pressure(
            self,
            x_var: str,
            png_name: str = None,
            use_slices: bool = True,
            slice_dirn: str = None,
            y_agg=None,
            dump_probe0: bool = True):
        """
        Plot the pressure data from simulation

        Args:
          x_var:
            Variable to plot against. "time" for time, "velocity" for velocity.
          png_name:
            (OPTIONAL) Name of the png file to save the plot. If not specified, the filename is selected automatically
          use_slices:
            (OPTIONAL) Whether to use slices or probes. If True, reads the pressure data from the cuttingPlane directory.
            Default is True.
          slice_dirn:
            (OPTIONAL) Direction of the slice. "z" for z-normal slices, "y" for y-normal slices. If use_slices is False, this argument is ignored.
          y_agg:
            (OPTIONAL) Aggregation method for y-normal slices. "cdf_median" for median of CDF, "mean" for mean, "median" for median.
            If not specified, defaults to None.
          dump_probe0:
            (OPTIONAL) Whether to dump the probe0 data to a numpy file for further analysis. Default is True.

        Raises:
          ValueError: If use_slices is True and slice_dirn is not "z" or "y".
          ValueError: If y_agg is not "cdf_median", "mean" or "median".
        """
        if slice_dirn == "y":
            warn("Aggregating pressure data using y-slices can yield inaccurate results. Use with caution")

        plot_suffix = "slices" if use_slices else "probes"
        pressure_df = self._probe2df(
            use_slices=use_slices,
            slice_dirn=slice_dirn,
            y_agg=y_agg)

        if x_var == "time":
            pressure_df.plot(
                xlabel="Time (s)",
                figsize=(
                    20,
                    10),
                fontsize=20,
                ylabel="Pressure (Pa)",
                title=f"Pressure at {plot_suffix}")
            plt.savefig(
                self.plots_dir +
                png_name +
                ".png") if png_name else plt.savefig(
                self.plots_dir +
                "probe_pressure.png")

        elif x_var == "velocity":
            plt.figure(figsize=[20, 10])
            self._calc_vel(df=pressure_df)

            vel_plot_df = pressure_df.groupby(["direction", "V_z"]).mean()

            # Sort the data for plotting
            vel_up = vel_plot_df[
                vel_plot_df.index.get_level_values(level='direction').isin(["up", "max"])
            ].reset_index('direction', drop=True).sort_index()

            vel_down = vel_plot_df[
                vel_plot_df.index.get_level_values(level='direction').isin(["down", "max"])
            ].reset_index('direction', drop=True).sort_index()

            if slice_dirn == "z":
                for i in range(self.nprobes):
                    plt.plot(vel_up.index,
                             vel_up[pressure_df.columns[i]],
                             label=f"Probe {i} (Up)",
                             color=f'C{i}',
                             marker='o')
                    plt.plot(vel_down.index,
                             vel_down[pressure_df.columns[i]],
                             label=f"Probe {i} (Down)",
                             color=f'C{i}',
                             marker='o',
                             linestyle='dashed')
                if dump_probe0:
                    probe0_up_v = vel_up.index.to_numpy()
                    probe0_up_p = vel_up[pressure_df.columns[0]].to_numpy()
                    probe0_down_v = vel_down.index.to_numpy()
                    probe0_down_p = vel_down[pressure_df.columns[0]].to_numpy()
                    probe0_2d = np.vstack(
                        (probe0_up_v, probe0_up_p, probe0_down_v, probe0_down_p)).T
                    np.save(
                        os.path.join(
                            self.plots_dir,
                            "probe0_plot_voidfrac.npy"),
                        probe0_2d)

            else:

                plt.plot(
                    vel_up.index,
                    vel_up['pressure'],
                    label=r"$V_z$ Increasing",
                    color='C0',
                    marker='o')
                plt.plot(
                    vel_down.index,
                    vel_down['pressure'],
                    label=r"$V_z$ Increasing",
                    color='C0',
                    marker='o',
                    linestyle='dashed')
            plt.xlabel("Velocity (m/s)")
            plt.ylabel("Pressure (Pa)")

        plt.legend()
        plt.title(f"Pressure vs Velocity for {plot_suffix}")

        plt.savefig(
            self.plots_dir +
            f"{png_name}.png") if png_name else plt.savefig(
            self.plots_dir +
            f"pressure_vel_plot_{plot_suffix}.png")

    def _read_voidfrac(self, post_dir: str, slice_dirn: str) -> pd.DataFrame:
        """
        Helper function to read the void fraction data using `pyvista` module.
        Args:
          post_dir:
            Path to the postprocessing directory containing the void fraction data.
          slice_dirn:
            Direction of the slice. "z" for z-normal slices, "y" for y-normal slices.

        Returns:
          A pandas DataFrame with the void fraction data indexed by time.

        Raises:
          ValueError: If slice_dirn is not "z" or "y".
        """

        times = os.listdir(post_dir)

        voidfrac_dict = dict()
        for dir in times:
            voidfrac_lst = []

            if slice_dirn == 'z':
                void_data = pv.read(
                    post_dir + '/' + dir + '/voidfraction_zNormal0.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_lst.append(voidfrac_arr.mean().item())

                void_data = pv.read(
                    post_dir + '/' + dir + '/voidfraction_zNormal1.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_lst.append(voidfrac_arr.mean().item())

                void_data = pv.read(
                    post_dir + '/' + dir + '/voidfraction_zNormal2.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_lst.append(voidfrac_arr.mean().item())

                void_data = pv.read(
                    post_dir + '/' + dir + '/voidfraction_zNormal3.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_lst.append(voidfrac_arr.mean().item())

                void_data = pv.read(
                    post_dir + '/' + dir + '/voidfraction_zNormal4.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_lst.append(voidfrac_arr.mean().item())
                voidfrac_dict[float(dir)] = voidfrac_lst
                voidfrac_df = pd.DataFrame.from_dict(
                    voidfrac_dict, orient='index', columns=[
                        f"Probe {i}" for i in range(
                            self.nprobes)]).sort_index()

            elif slice_dirn == 'y':
                void_data = pv.read(
                    post_dir + '/' + dir + '/voidfraction_yNormal.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_dict[float(dir)] = voidfrac_arr
                voidfrac_dict = {
                    k: self.find_cdfmedian(v) for k,
                    v in voidfrac_dict.items()}
                voidfrac_df = pd.DataFrame.from_dict(
                    voidfrac_dict, orient='index', columns=['void_frac']).sort_index()
            else:
                raise ValueError("Invalid slice direction. Choose 'z' or 'y'")

        return voidfrac_df

    def plot_voidfrac(
            self,
            slice_dirn: str,
            x_var: str,
            post_dir: str = "CFD/postProcessing/cuttingPlane/",
            png_name: str = None,
            dump_probe0: bool = True):
        """
        Plots the void fraction data against time or velocity, depending on the x_var argument.
        If x_var is "time", the plot will show void fraction vs time. If x_var is "velocity",
        the plot will show void fraction vs velocity. The velocity is calculated from the velocity config file.

        Args:
          slice_dirn:
            Direction of the slice. "z" for z-normal slices, "y" for y-normal slices.
          x_var:
            Variable to plot against. "time" for time, "velocity" for velocity.
          post_dir:
            (OPTIONAL) Path to the postprocessing directory containing the void fraction data. Default is "CFD/postProcessing/cuttingPlane/".
          png_name:
            (OPTIONAL) Name of the png file to save the plot. If not specified, the filename is selected automatically
          dump_probe0:
            (OPTIONAL) Whether to dump the probe0 data to a numpy file for further analysis. Default is True.

        Raises:
          ValueError: If slice_dirn is not "z" or "y".
          ValueError: If x_var is not "time" or "velocity".
        """
        plt.figure(figsize=[20, 10])
        voidfrac_df = self._read_voidfrac(
            slice_dirn=slice_dirn, post_dir=post_dir)

        if x_var == "time":
            voidfrac_df.plot(
                xlabel="Time (s)",
                ylabel="Void Fraction (-)",
                title="Void Fraction vs Time")
            plt.savefig(
                self.plots_dir +
                f"{png_name}.png") if png_name else plt.savefig(
                self.plots_dir +
                f"voidfrac_time_plot_{slice_dirn}.png")
            plt.xlabel("Velocity (m/s)")

        elif x_var == "velocity":
            self._calc_vel(df=voidfrac_df)

            vel_plot_df = voidfrac_df.groupby(["direction", "V_z"]).mean()

            # Sort the data for plotting
            vel_up = vel_plot_df[
                vel_plot_df.index.get_level_values(level='direction').isin(["up", "max"])
            ].reset_index('direction', drop=True).sort_index()
            # print("Vel up", vel_up)

            vel_down = vel_plot_df[
                vel_plot_df.index.get_level_values(level='direction').isin(["down", "max"])
            ].reset_index('direction', drop=True).sort_index()

            if slice_dirn == "z":
                for i in range(self.nprobes):
                    plt.plot(vel_up.index,
                             vel_up[voidfrac_df.columns[i]],
                             label=f"Probe {i} (Up)",
                             color=f'C{i}',
                             marker='o')
                    plt.plot(vel_down.index,
                             vel_down[voidfrac_df.columns[i]],
                             label=f"Probe {i} (Down)",
                             color=f'C{i}',
                             marker='o',
                             linestyle='dashed')
            else:
                plt.plot(
                    vel_up.index,
                    vel_up['void_frac'],
                    label=r"$V_z$ Increasing",
                    color='C0',
                    marker='o')
                plt.plot(
                    vel_down.index,
                    vel_down['void_frac'],
                    label=r"$V_z$ Increasing",
                    color='C0',
                    marker='o',
                    linestyle='dashed')
                if dump_probe0:
                    probe0_up_v = vel_up.index.to_numpy()
                    probe0_up_p = vel_up['void_frac'].to_numpy()
                    probe0_down_v = vel_down.index.to_numpy()
                    probe0_down_p = vel_down['void_frac'].to_numpy()
                    probe0_2d = np.vstack(
                        (probe0_up_v, probe0_up_p, probe0_down_v, probe0_down_p)).T
                    np.save(
                        os.path.join(self.plots_dir, "probe0_plot_P.npy"),
                        probe0_2d
                    )

            plt.xlabel("Velocity (m/s)")
            plt.ylabel("Void Fraction (-)")

        plt.legend()
        plt.title("Void Fraction vs Velocity")

        plt.savefig(
            self.plots_dir +
            f"{png_name}.png") if png_name else plt.savefig(
            self.plots_dir +
            f"voidfrac_vel_plot_{slice_dirn}.png")

    def _read_collisions(self, csv_path: str, calltype: str) -> pd.DataFrame:
        """
        Helper function to read the collision data from the DEM simulation

        Args:
          csv_path:
            Path to the csv file containing the collision data.
          calltype:
            Type of data to return. "contactarea" for contact area per atom, "contactn" for contact number per atom.

        Returns:
          A pandas DataFrame with the collision data. If calltype is "contactarea", the DataFrame will contain a_contact_peratom column.
          If calltype is "contactn", the DataFrame will contain a contactn column.
        """
        df = pd.read_csv(csv_path, sep='\\s+')

        if calltype == "contactarea":
            df['a_contact_peratom'] = df.a_contact / df.n_atoms
            return df.drop(columns=['n_atoms'])
        elif calltype == "contactn":
            df['contactn'] = df.n_contact / df.n_atoms
            return df

    def plot_contactarea(
            self,
            csv_path: str = 'DEM/post/collisions.csv',
            png_name: str = None):
        """
        Plot the contact area data from the DEM simulations

        Args:
          csv_path:
            Path to the csv file containing the collision data. Default is 'DEM/post/collisions.csv'.
          png_name:
            (OPTIONAL) Name of the png file to save the plot. If not specified, the filename is selected automatically.
        """
        fig = plt.figure(figsize=[20, 10])

        contact_df = self._read_collisions(
            csv_path=csv_path, calltype="contactarea")
        plt.plot(
            contact_df.time,
            contact_df.a_contact_peratom,
            label="Contact Area per Atom",
            color='C0')
        plt.xlabel("Time (s)")
        plt.ylabel(r"Contact Area per Atom ($m^2$)")
        plt.title("Contact Area vs Time")
        if png_name:
            plt.savefig(self.plots_dir + f"{png_name}.png")
        else:
            plt.savefig(self.plots_dir + "contactarea_time_plot.png")


def liggghts2vtk(
        timestep: float = 5e-6,
        vtk_dir: str = None,
        dump_every: int = None,
        liggghts_dump_dir: str = None,
        file_suffix: str = ".liggghts_run"
):
    """
    Convert LIGGGHTS dump files to VTK format for visualization using the `pyvetk` module.

    Args:
      timestep:
        (OPTIONAL) Time step of the simulation in seconds. Default is 5e-6.
      vtk_dir:
        (OPTIONAL) Directory to save the VTK files. If not specified, defaults to "DEM/post/vtk/".
      dump_every:
        (OPTIONAL) Dump every nth file. If not specified, all files are processed.
      liggghts_dump_dir:
        (OPTIONAL) Directory containing the LIGGGHTS dump files. If not specified, defaults to "DEM/post/".
      file_suffix:
        (OPTIONAL) Suffix of the LIGGGHTS dump files. Default is ".liggghts_run".
    """

    if vtk_dir is None:
        vtk_dir = "DEM/post/vtk/"
    if not os.path.isdir(vtk_dir):
        os.makedirs(vtk_dir, exist_ok=True)

    if not liggghts_dump_dir:
        liggghts_dump_dir = "DEM/post/"
    elif not os.path.isdir(liggghts_dump_dir):
        raise Exception(
            f"LIGGGHTS dump directory {liggghts_dump_dir} does not exist")

    def is_liggghts_file(filename): return filename.endswith(file_suffix)
    liggghts_files = os.listdir(liggghts_dump_dir)
    liggghts_files = list(
        filter(is_liggghts_file, liggghts_files)
    )

    if not liggghts_files:
        raise Exception(f"No LIGGGHTS dump files found in {liggghts_dump_dir}")

    def sort_key(filename):
        match = re.search(r'\d+', filename)
        return int(match.group(0)) if match else 0
    liggghts_files.sort(key=sort_key)

    if dump_every:
        liggghts_files = liggghts_files[::dump_every]

    with open(os.path.join(liggghts_dump_dir, liggghts_files[0]), 'r') as f:
        header = f.read().split('\n')[8]
        header = header.split()[2:]

    dump_diff = sort_key(liggghts_files[1]) - sort_key(liggghts_files[0])
    step_diff = sort_key(liggghts_files[0]) - dump_diff

    for file in liggghts_files:

        vtk_name = file.split(file_suffix)[0]
        dump_df = pd.read_csv(
            os.path.join(liggghts_dump_dir, file),
            skiprows=9, sep=' ', header=None,
            engine='pyarrow'
        ).iloc[:, :-1]
        dump_df.columns = header

        x = np.array(dump_df["x"].values, dtype=np.float64)
        y = np.array(dump_df["y"].values, dtype=np.float64)
        z = np.array(dump_df["z"].values, dtype=np.float64)
        velocity = np.array(np.linalg.norm(
            dump_df[['vx', 'vy', 'vz']].values, axis=1), dtype=np.float64)
        velocity_x = np.array(dump_df['vx'].values, dtype=np.float64)
        velocity_y = np.array(dump_df['vy'].values, dtype=np.float64)
        velocity_z = np.array(dump_df['vz'].values, dtype=np.float64)
        radius = np.array(dump_df['radius'].values, dtype=np.float64)

        time = (sort_key(file) - step_diff) * timestep

        pointsToVTK(
            os.path.join(vtk_dir, vtk_name),
            x, y, z,
            data=dict(
                time=np.full(len(x), time, dtype=np.float64),
                radius=radius,
                velocity=velocity,
                velocity_x=velocity_x,
                velocity_y=velocity_y,
                velocity_z=velocity_z
            )
        )


if __name__ == "__main__":

    # Example usage - usable as a standalone script for default paths

    pressure_path = 'CFD/postProcessing/cuttingPlane/'
    velcfg_path = 'velcfg.txt'

    probe_cfdem_slices = ProbeAnalysis(
        pressure_path=pressure_path,
        nprobes=5,
        velcfg_path=velcfg_path,
        dump2csv=False
    )

    """
    Z-Normal Slices vs Time
    """
    # probe_cfdem_slices.plot_pressure(slice_dirn="z", x_var="time", png_name="pressure_time_plot_z", use_slices=True)
    # probe_cfdem_slices.plot_voidfrac(slice_dirn="z", x_var="time", png_name="voidfrac_time_plot_z")

    """
    Y-Normal Slices vs Velocity
    """
    # probe_cfdem_slices.plot_pressure(slice_dirn="y", x_var="velocity", png_name="pressure_vel_plot_y", use_slices=True, y_agg='median')
    # probe_cfdem_slices.plot_voidfrac(slice_dirn="y", x_var="velocity", png_name="voidfrac_vel_plot_y")

    """
    Y-Normal Slices vs Time
    """
    # probe_cfdem_slices.plot_pressure(slice_dirn="y", x_var="time", png_name="pressure_time_plot_y", use_slices=True, y_agg='median')
    # probe_cfdem_slices.plot_voidfrac(slice_dirn="y", x_var="time", png_name="voidfrac_time_plot_y")

    """
    Z-normal Slices vs Velocity
    """
    # probe_cfdem_slices.plot_pressure(slice_dirn="z", x_var="velocity", png_name="pressure_vel_plot_z", use_slices=True)
    # probe_cfdem_slices.plot_voidfrac(slice_dirn="z", x_var="velocity", png_name="voidfrac_vel_plot_z")

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
