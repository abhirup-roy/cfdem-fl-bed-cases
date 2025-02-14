#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px

class ProbeAnalysis():

    def __init__(self, pressure_path, nprobes, velcfg_path, dump2csv=True, plots_dir='../../plots/'):
        """
        pressure_path: Path to the pressure data
        nprobes: Number of probes
        velcfg_path: Path to the vel_cfg file
        dump2csv: Save the probe data to a csv file
        plots_dir: Directory to save the plots
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

        self._probe2df()
    
    def _probe2df(self):
        """
        Convert the probe data to a pandas dataframe, with data indexed by time. Dump if specified
        """
        # Make df from the probe data
        headers = ["Probe Time"]
        for i in range(self.nprobes):
            headers.append(f"Probe {i}")

        self.probe_df = pd.read_csv(
            self.pressure_path, 
            delim_whitespace=True, 
            comment='#', 
            names=headers, 
            header=None,
        ).set_index("Probe Time")

        # Dump to csv if specified
        if self.dump2csv:
            self.probe_df.to_csv("probe_pressure.csv")

    def plot_probe(self, plot_backend):
        """
        Plot the pressure data at the probes ported to matplotlib or plotly

        plot_backend: "matplotlib" for png or "plotly" for interactive html
        """
        # Plotly backend
        pd.options.plotting.backend = plot_backend

        if plot_backend == "plotly":
            probe_px = self.probe_df.plot(title="Pressure at Probes", template="simple_white",
                        labels=dict(index="Time (s)", value="Pressure", variable="Probe"))
            probe_px.write_html(self.plots_dir+"probe_pressure.html")

        elif plot_backend == "matplotlib":
            plt.figure(figsize=[30,20])
            self.probe_df.plot(xlabel="Time (s)", ylabel="Pressure (Pa)", title="Pressure at Probes")
            plt.savefig(self.plots_dir+"probe_pressure.png")
    
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

    def _calc_vel(self):
        """
        Map the velocity to pressure
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
        vz_arr = np.zeros_like(self.probe_df.index.to_numpy())

        # Map the velocity to the pressure data
        for i in range(len(bounds)):
            mask = (self.probe_df.index.to_numpy() > lb[i]) & (self.probe_df.index.to_numpy() < ub[i])
            vz_arr[mask] = vel[i]

            if i < len(bounds) - 1:
                gap_mask = (self.probe_df.index.to_numpy() >= ub[i]) & (self.probe_df.index.to_numpy() <= lb[i + 1])
                vz_arr[gap_mask] = np.nan
            
        self.probe_df["V_z"] = vz_arr


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
        self.probe_df["direction"] = self.probe_df.index.to_series().apply(_map_direction)

        


    def plot_pressure_v(self):
        """
        Plot the pressure data at the probes against the velocity
        """
        self._calc_vel()

        # Create the plot
        plt.figure(figsize=[20,10])
        
        # Use groupby to get the mean pressure at each velocity for increasing, decreasing and max velocity
        vel_plot_df = self.probe_df.groupby(["direction", "V_z"]).mean()
        # Set the pressure at zero velocity to zero
        vel_plot_df.loc[vel_plot_df.index.get_level_values('V_z') == 0, :] = 0
        
        # Sort the data for plotting
        vel_up = vel_plot_df[
            vel_plot_df.index.get_level_values(level='direction').isin(["up", "max"])
        ].reset_index('direction', drop=True
        ).sort_index()
        
        vel_down = vel_plot_df[
            vel_plot_df.index.get_level_values(level='direction').isin(["down", "max"])
        ].reset_index('direction', drop=True
        ).sort_index()
        
        for i in range(self.nprobes):
            plt.plot(vel_up.index, vel_up[self.probe_df.columns[i]], label=self.probe_df.columns[i], color=f'C{i}', marker='o')
            plt.plot(vel_down.index, vel_down[self.probe_df.columns[i]], label=self.probe_df.columns[i], color=f'C{i}', marker='o', linestyle='dashed')

        plt.xlabel("Velocity (m/s)")
        plt.ylabel("Pressure (Pa)")

        plt.savefig(self.plots_dir + "pressure_vel_plot.png")
    


if __name__ == "__main__":
    pressure_path = '../../CFD/postProcessing/probes/0/p'
    velcfg_path = 'velcfg.txt'

    probe_cfdem = ProbeAnalysis(
        pressure_path=pressure_path,
        nprobes=5,
        velcfg_path=velcfg_path,
        dump2csv=True
    )
    
    probe_cfdem.plot_probe(plot_backend="plotly")
    probe_cfdem.plot_probe(plot_backend="matplotlib")
    probe_cfdem.plot_pressure_v()