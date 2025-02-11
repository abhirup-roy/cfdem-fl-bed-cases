#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px

class ProbeAnalysis():

    def __init__(self, probe_path, nprobes, probes_text_path, dump2csv=True):
        """
        probe_path: Path to the probe file
        nprobes: Number of probes
        probes_text_path: Path to the probes.txt file
        dump2csv: Save the probe data to a csv file
        """
        self.probe_path = probe_path
        self.nprobes = nprobes
        self.dump2csv = dump2csv
        self.probes_text_path = probes_text_path

        self._probe2df()
    
    def _probe2df(self):
        if not os.path.exists(probe_path):
            raise Exception(f"Path {probe_path} does not exist")

        headers = ["Probe Time"]
        for i in range(self.nprobes):
            headers.append(f"Probe {i}")

        self.probe_df = pd.read_csv(
            probe_path, 
            delim_whitespace=True, 
            comment='#', 
            names=headers, 
            header=None,
        ).set_index("Probe Time")

        if self.dump2csv:
            self.probe_df.to_csv("probe_pressure.csv")
    def plot_probe(self, plot_backend = "plotly"):

        pd.options.plotting.backend = plot_backend

        if plot_backend == "plotly":
            probe_px = self.probe_df.plot(title="Pressure at Probes", template="simple_white",
                        labels=dict(index="Time (s)", value="Pressure", variable="Probe"))
            probe_px.write_html("probe_pressure.html")

        elif plot_backend == "matplotlib":
            plt.figure(figsize=[30,20])
            self.probe_df.plot(xlabel="Time (s)", ylabel="Pressure (Pa)", title="Pressure at Probes")
            plt.savefig("probe_pressure.png")
    
    def _read_probetxt(self):
        with open(self.probes_text_path, "r") as f:
            probe_text = f.read().splitlines(False)

            self.t = []
            self.v_z = []

            for line in probe_text:
                line_splt = line.replace("(", "").replace(")", "").split()
                self.t.append(float(line_splt[0]))
                self.v_z.append(float(line_splt[-1]))
            print("Selected times: ", self.t)
            print("Corresponding vel: ", self.v_z)

    def _calc_vel(self):
        bounds = []
        vel =[]
        
        self._read_probetxt()
        
        for i in range(len(self.t)-1):
            if self.v_z[i] == self.v_z[i+1]:
                bounds.append([self.t[i], self.t[i+1]])
                vel.append(self.v_z[i])
            else:
                pass
        
        # def _vel_mask(df, bounds, vel):
        #     for i in range(len(bounds)):
        #         if df["Probe Time"] < bounds[i][0] and df["Probe Time"] > bounds[i][1]:
        #             df["V_z"] = vel[i]
        #         else: pass
        # self.probe_df = self.probe_df.apply(_vel_mask, args=(bounds, vel), axis=1)

        lb = [b[0] for b in bounds]
        ub = [b[1] for b in bounds]
        try:
            vz_arr = np.zeros_like(self.probe_df.index.to_numpy())
        
            for i in range(len(bounds)):
                mask = (self.probe_df.index.to_numpy() > lb[i]) & (self.probe_df.index.to_numpy() < ub[i])
                vz_arr[mask] = vel[i]
        except:
            print(self.probe_df)

        self.probe_df["V_z"] = vz_arr


    def plot_pressure_v(self):
        self._calc_vel()

        print(self.probe_df)
        plt.figure(figsize=[30,20])
        for i in range(self.nprobes):
            plt.scatter(self.probe_df["V_z"], self.probe_df[f"Probe {i}"], label=f"Probe {i}",)
            
        plt.legend()
        plt.savefig("pressure_vel_plot.png")
    


if __name__ == "__main__":
    probe_path = '../../CFD/postProcessing/probes/0/p'
    probes_text_path = 'probes_txt'

    probe_cfdem = ProbeAnalysis(probe_path=probe_path, nprobes=5, probes_text_path=probes_text_path, dump2csv=True)
    probe_cfdem.plot_probe(plot_backend="plotly")
    probe_cfdem.plot_probe(plot_backend="matplotlib")
    probe_cfdem.plot_pressure_v()