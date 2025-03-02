
#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd

import plotly.express as px
import pyvista as pv

class ProbeAnalysis():

    def __init__(self, pressure_path, nprobes, velcfg_path, dump2csv=True, plots_dir='../../plots/'):
        """
        STRING  pressure_path: Path to the pressure data. If using slices, point to cuttingPlane directory. 
        INT     nprobes: Number of probes
        STRING  velcfg_path: Path to the vel_cfg file
        BOOL    dump2csv: Save the probe data to a csv file
        STRING  plots_dir: Directory to save the plots
        """
        # Check if the paths exist
        if not os.path.exists(pressure_path):
            raise Exception(f"Pressure data at {pressure_path} does not exist")
                
        if not os.path.exists(velcfg_path):
            raise Exception(f"Velocity config at {velcfg_path} does not exist")
        if not os.path.isdir(plots_dir):
            raise Exception(f" Plots directory at {plots_dir} does not exist")

        self.pressure_path = pressure_path
        self.nprobes = nprobes
        self.dump2csv = dump2csv
        self.velcfg_path = velcfg_path
        self.plots_dir = plots_dir

        rcParams.update({'font.size': 20})
    
    def _probe2df(self, use_slices, slice_dirn):
        """
        Convert the probe data to a pandas dataframe, with data indexed by time. Dump if specified
        """
        if use_slices:
            times = os.listdir(self.pressure_path)
            pressure_dict = dict()

            for dir in times:
                if slice_dirn == "z":
                    p_lst = []
                    pdata = pv.read(self.pressure_path + '/' + dir + '/p_zNormal0.vtk')
                    p_arr = pdata.get_array('p')
                    p_lst.append(p_arr.mean().item())
                    
                    pdata = pv.read(self.pressure_path + '/' + dir + '/p_zNormal1.vtk')
                    p_arr = pdata.get_array('p')
                    p_lst.append(p_arr.mean().item())

                    pdata = pv.read(self.pressure_path + '/' + dir + '/p_zNormal2.vtk')
                    p_arr = pdata.get_array('p')
                    p_lst.append(p_arr.mean().item())

                    pdata = pv.read(self.pressure_path + '/' + dir + '/p_zNormal3.vtk')
                    p_arr = pdata.get_array('p')
                    p_lst.append(p_arr.mean().item())

                    pdata = pv.read(self.pressure_path + '/' + dir + '/p_zNormal4.vtk')
                    p_arr = pdata.get_array('p')
                    p_lst.append(p_arr.mean().item())

                    pressure_dict[float(dir)] = p_lst
                
                    pressure_df = pd.DataFrame.from_dict(
                        pressure_dict, orient='index',columns=[f"Probe {i}" for i in range(self.nprobes)]
                    ).sort_index()
 
            
                elif slice_dirn == "y":
                    pdata = pv.read(self.pressure_path+ '/' + dir + '/p_yNormal.vtk')
                    p_arr = pdata.get_array('p')
                    pressure_dict[float(dir)] = p_arr

                    pressure_dict = {k: self.find_cdfmedian(v) for k,v in pressure_dict.items()}
                    pressure_df = pd.DataFrame.from_dict(pressure_dict, orient='index', columns=['pressure']).sort_index()

                else:
                    raise Exception("Invalid slice direction. Choose 'z' or 'y'")
                                  
            
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
        Read the velocity config file
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

    def _calc_vel(self, df):
        """
        Helper function to map the velocity to pressure. Reads the velocity config file and maps the velocity to the time-series data.
        """
        # Initialise bounds and velocity
        bounds = []
        vel =[]
        
        self._read_probetxt()
        # Find the bounds for velocity
        for i in range(len(self.t)-1):
            if self.v_z[i] == self.v_z[i+1]:
                bounds.append([self.t[i], self.t[i+1]])
                vel.append(self.v_z[i])

                if self.v_z[i] == max(self.v_z):
                    max_vel_t1, max_vel_t2 = self.t[i], self.t[i+1]

            else:
                pass
        
        # Lower and upper time bounds for each velocity 
        lb = [b[0] for b in bounds]
        ub = [b[1] for b in bounds]
        vz_arr = np.zeros_like(df.index.to_numpy())

        # Map the velocity to the pressure data
        for i in range(len(bounds)):
            mask = (df.index.to_numpy() > lb[i]) & (df.index.to_numpy() < ub[i])
            vz_arr[mask] = vel[i]

            if i < len(bounds) - 1:
                gap_mask = (df.index.to_numpy() >= ub[i]) & (df.index.to_numpy() <= lb[i + 1])
                vz_arr[gap_mask] = np.nan
            
        df["V_z"] = vz_arr


        def _map_direction(x):
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

    def find_cdfmedian(self, arr):
                    x, counts = np.unique(arr, return_counts=True)
                    cusum = np.cumsum(counts)
                    cdf = cusum/cusum[-1]
        
                    median_idx = cdf.tolist().index(np.percentile(cdf,50,method='nearest'))
                    return x[median_idx].item()
 
    def plot_pressure(self, x_var, png_name=None, use_slices=True, slice_dirn=None):
        """
        Plot the pressure data from simulation.
        STRING x_var: Variable to plot against. "time" for time, "velocity" for velocity.
        STRING slice_dirn (optional): Direction of the slice. "z" for z-normal, "y" for y-normal. MUST be specified if use_slices is True!!!
        STRING png_name (optional): Name of the png file to save the plot. If not specified, the filename is selected automatically
        BOOL use_slices (optional): Use slices or probes. Default is True
        """
        fig = plt.figure(figsize=[20,10])
        plot_suffix = "slices" if use_slices else "probes"
        pressure_df = self._probe2df(use_slices=use_slices, slice_dirn=slice_dirn)

        if x_var == "time":

            pressure_df.plot(xlabel="Time (s)", ylabel="Pressure (Pa)", title=f"Pressure at {plot_suffix}")
            plt.savefig(self.plots_dir+png_name+".png") if png_name else plt.savefig(self.plots_dir+"probe_pressure.png")
        
        elif x_var == "velocity":
            self._calc_vel(df=pressure_df)

            vel_plot_df = pressure_df.groupby(["direction", "V_z"]).mean()
            
            # Sort the data for plotting
            vel_up = vel_plot_df[
                vel_plot_df.index.get_level_values(level='direction').isin(["up", "max"])
            ].reset_index('direction', drop=True).sort_index()
            print("Vel up", vel_up)
            
            vel_down = vel_plot_df[
                vel_plot_df.index.get_level_values(level='direction').isin(["down", "max"])
            ].reset_index('direction', drop=True).sort_index()
            
            if slice_dirn == "z":
                for i in range(self.nprobes):
                    plt.plot(vel_up.index, vel_up[pressure_df.columns[i]], label=f"Probe {i} (Up)", color=f'C{i}', marker='o')
                    plt.plot(vel_down.index, vel_down[pressure_df.columns[i]], label=f"Probe {i} (Down)", color=f'C{i}', marker='o', linestyle='dashed')
            else:
                plt.plot(vel_up.index, vel_up['pressure'], label=r"$V_z$ Increasing", color='C0', marker='o')
                plt.plot(vel_down.index, vel_down['pressure'], label=r"$V_z$ Increasing", color='C0', marker='o', linestyle='dashed')

        plt.xlabel("Velocity (m/s)")
        plt.ylabel("Pressure (Pa)")
        plt.legend()
        plt.title(f"Pressure vs Velocity for {plot_suffix}")

        plt.savefig(self.plots_dir + f"{png_name}.png") if png_name else plt.savefig(self.plots_dir + f"pressure_vel_plot_{plot_suffix}.png")
            
    
    def _read_voidfrac(self, post_dir, slice_dirn):
        """
        Helper function to read the void fraction data using `pyvista`
        """

        times = os.listdir(post_dir)
        
        voidfrac_dict = dict()
        for dir in times:
            voidfrac_lst = []

            if slice_dirn=='z':
                void_data = pv.read(post_dir + '/' + dir + '/voidfraction_zNormal0.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_lst.append(voidfrac_arr.mean().item())
                
                void_data = pv.read(post_dir + '/' + dir + '/voidfraction_zNormal1.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_lst.append(voidfrac_arr.mean().item())

                void_data = pv.read(post_dir + '/' + dir + '/voidfraction_zNormal2.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_lst.append(voidfrac_arr.mean().item())

                void_data = pv.read(post_dir + '/' + dir + '/voidfraction_zNormal3.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_lst.append(voidfrac_arr.mean().item())

                void_data = pv.read(post_dir + '/' + dir + '/voidfraction_zNormal4.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_lst.append(voidfrac_arr.mean().item())
                voidfrac_dict[float(dir)] = voidfrac_lst
                voidfrac_df = pd.DataFrame.from_dict(
                    voidfrac_dict, orient='index', columns=[f"Probe {i}" for i in range(self.nprobes)]
                ).sort_index()


            elif slice_dirn=='y':
                void_data = pv.read(post_dir + '/' + dir + '/voidfraction_yNormal.vtk')
                voidfrac_arr = void_data.get_array('voidfraction')
                voidfrac_dict[float(dir)] = voidfrac_arr
                voidfrac_dict = {k: self.find_cdfmedian(v) for k,v in voidfrac_dict.items()}
                voidfrac_df = pd.DataFrame.from_dict(voidfrac_dict, orient='index', columns=['void_frac']).sort_index()
            else:
                raise Exception("Invalid slice direction. Choose 'z' or 'y'")
            
        return voidfrac_df
    
    def plot_voidfrac(self, slice_dirn, x_var, post_dir="../../CFD/postProcessing/cuttingPlane/", png_name=None):
        """
        Plot the void fraction data
        STRING slice_dirn: Method to plot the void fraction data. "slices" for z-normal slices, "cdf_median" for median of CDF of y-normal slice
        STRING x_var: Variable to plot against. "time" for time, "velocity" for velocity.
        STRING post_dir (optional): Path to the postprocessing directory (Default: CFD/postProcessing/cuttingPlane/)
        STRING png_name (optional): Name of the png file to save the plot. If not specified, the filename is selected automatically

        """
        fig = plt.figure(figsize=[20,10])
        voidfrac_df = self._read_voidfrac(slice_dirn=slice_dirn, post_dir=post_dir)

        if x_var == "time":
            voidfrac_df.plot(xlabel="Time (s)", ylabel="Void Fraction (-)", title="Void Fraction vs Time")
            plt.savefig(self.plots_dir + f"{png_name}.png") if png_name else plt.savefig(self.plots_dir + f"voidfrac_time_plot_{slice_dirn}.png")
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
                    plt.plot(vel_up.index, vel_up[voidfrac_df.columns[i]], label=f"Probe {i} (Up)", color=f'C{i}', marker='o')
                    plt.plot(vel_down.index, vel_down[voidfrac_df.columns[i]], label=f"Probe {i} (Down)", color=f'C{i}', marker='o', linestyle='dashed')
            else:
                plt.plot(vel_up.index, vel_up['void_frac'], label=r"$V_z$ Increasing", color='C0', marker='o')
                plt.plot(vel_down.index, vel_down['void_frac'], label=r"$V_z$ Increasing", color='C0', marker='o', linestyle='dashed')

            plt.xlabel("Velocity (m/s)")
            plt.ylabel("Void Fraction (-)")
            plt.legend()
            plt.title("Void Fraction vs Velocity")
            
        #     for i in range(self.nprobes):
        #         if slice_dirn == "y":
        #             plt.plot(vel_up.index, vel_up['void_frac'], label=f"Probe {i} (Up)", color=f'C{i}', marker='o')
        #             plt.plot(vel_down.index, vel_down['void_frac'], label=f"Probe {i} (Down)", color=f'C{i}', marker='o', linestyle='dashed')
        #         else:
        #             plt.plot(vel_up.index, vel_up[voidfrac_df.columns[i]], label=f"Probe {i} (Up)", color=f'C{i}', marker='o')
        #             plt.plot(vel_down.index, vel_down[voidfrac_df.columns[i]], label=f"Probe {i} (Down)", color=f'C{i}', marker='o', linestyle='dashed')
                

        #         plt.xlabel("Velocity (m/s)")
        #         plt.ylabel("Void Fraction (-)")
        #         plt.legend()
        #         plt.title("Void Fraction vs Velocity (Z-Normal Slices)")
        #         plt.savefig(self.plots_dir + f"{png_name}.png") if png_name else plt.savefig(self.plots_dir + "voidfrac_vel_plot_z0-4.png")
                
        # if slice_dirn=="y":
        #     voidfrac_dict = self._read_voidfrac(post_dir, slice_dirn=slice_dirn)
                
        #     voidfrac_dict = {k: self.find_cdfmedian(v) for k,v in voidfrac_dict.items()}
        #     voidfrac_df = pd.DataFrame.from_dict(voidfrac_dict, orient='index', columns=['void_frac']).sort_index()
        #     if x_var == "time":
        #         voidfrac_df.plot(xlabel="Time (s)", ylabel="Void Fraction (-)", title="Void Fraction vs Time (Y-Normal Slice)")
        #         plt.savefig(self.plots_dir + f"{png_name}.png") if png_name else plt.savefig(self.plots_dir + "voidfrac_time_plot_y.png") 
        #     elif x_var == "velocity":
        #         self._calc_vel(df=voidfrac_df)

        #         vel_plot_df = voidfrac_df.groupby(["direction", "V_z"]).mean()

        #         # Sort the data for plotting
        #         vel_up = vel_plot_df[
        #             vel_plot_df.index.get_level_values(level='direction').isin(["up", "max"])
        #         ].reset_index('direction', drop=True).sort_index()

        #         vel_down = vel_plot_df[
        #             vel_plot_df.index.get_level_values(level='direction').isin(["down", "max"])
        #         ].reset_index('direction', drop=True).sort_index()
                
        #         plt.plot(vel_up.index, vel_up['void_frac'], label=r"$V_z$ Increasing", color='C0', marker='o')
        #         plt.plot(vel_down.index, vel_down['void_frac'], label=r"$V_z$ Increasing", color='C0', marker='o', linestyle='dashed')

        #     plt.xlabel("Velocity (m/s)")
        #     plt.ylabel("Void Fraction (-)")
        #     plt.legend()
        #     plt.title("Void Fraction vs Velocity (Y-Normal Slice)")
        #     plt.savefig(self.plots_dir + f"{png_name}.png") if png_name else plt.savefig(self.plots_dir + "voidfrac_vel_plot_y.png")




if __name__ == "__main__":

    pressure_path = '../../CFD/postProcessing/cuttingPlane/'
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
    probe_cfdem_slices.plot_pressure(slice_dirn="z", x_var="time", png_name="pressure_time_plot_z", use_slices=True)
    probe_cfdem_slices.plot_voidfrac(slice_dirn="z", x_var="time", png_name="voidfrac_time_plot_z")

    """
    Y-Normal Slices vs Velocity
    """
    probe_cfdem_slices.plot_pressure(slice_dirn="y", x_var="velocity", png_name="pressure_vel_plot_y", use_slices=True)
    probe_cfdem_slices.plot_voidfrac(slice_dirn="y", x_var="velocity", png_name="voidfrac_vel_plot_y")
