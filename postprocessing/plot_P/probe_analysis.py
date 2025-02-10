#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px

"""
probe_path = '../../CFD/postProcessing/probes/0/p'

# Check if correct path is provided
if not os.path.exists(probe_path):
    raise Exception(f"Path {probe_path} does not exist")

headers = ["Probe Time", "Probe 0", "Probe 1", "Probe 2", "Probe 3", "Probe 4"]

# Extract pressure at probes
probe_df = pd.read_csv(
    probe_path, 
    delim_whitespace=True, 
    comment='#', 
    names=headers, 
    header=None,
).set_index("Probe Time")

probe_df.to_csv("probe_pressure.csv")

#Matplotlib plotting

plt.figure(figsize=[20,10])

probe_df.plot(xlabel="Time (s)", ylabel="Pressure (Pa)", title="Pressure at Probes")
plt.savefig("probe_pressure.png")
"""

# pd.options.plotting.backend = "plotly"
# times_lst = [0,0.5, 0.6, 1.6, 1.7, 2.7, 2.8, 3.8, 3.9, 4.9, 5, 6,
#              6.1, 7.1, 7.2, 8.2, 8.3, 9.3, 9.4, 10.4, 10.5, 11.5,
#              12.5, 12.6, 13.6, 13.7]

# prob_px = px.line(probe_df, title="Pressure at Probes")
# prob_px.update_xaxes(title_text="Time (s)")
# prob_px.update_yaxes(title_text="Pressure (Pa)")
# prob_px.update_layout(hovermode="x unified")
# prob_px.write_html("probe_pressure.html")


# probe_px = probe_df.plot(title="Pressure at Probes", template="simple_white",
#               labels=dict(index="Time (s)", value="Pressure", variable="Probe"))

# for time in times_lst:
#     probe_px.add_vline(x=time, line_width=3, line_dash="dash", line_color="green")


# probe_px.write_html("probe_pressure.html")

class Probe():

    def __init__(self, probe_path, nprobes, probes_text_path, dump2csv=True):
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

    def _calc_vel(self):
        bounds = []
        vel =[]
        
        self._read_probetxt()
        
        for i in range(len(self.t)-1):
            if self.t[i] == self.t[i+1]:
                bounds.append([self.t[i], self.t[i+1]])
                vel.append(self.v_z[i])
            else:
                pass
        
        def _vel_mask(df, bounds, vel):
            for i in range(len(bounds)):
                if df["Probe Time"] < bounds[i][0] and df["Probe Time"] > bounds[i][1]:
                    df["V_z"] = vel[i]
                else: pass
        
        self.probe_df = self.probe_df.apply(_vel_mask, args=(bounds, vel), axis=1)

    def plot_pressure_v(self):
        self._calc_vel()
        
        print(self.probe_df)

        plt.plot(x=self.probe_df["V_z"], xlabel="Time (s)", ylabel="Velocity (m/s)")
        plt.savefig("pressure_vel_plot.png")
    


if __name__ == "__main__":
    probe_path = '../../CFD/postProcessing/probes/0/p'
    probes_text_path = 'probes.txt'

    probe_cfdem = Probe(probe_path=probe_path, nprobes=5, probes_text_path=probes_text_path, dump2csv=True)
    # print(probe_cfdem.probe_df)
    probe_cfdem.plot_probe(plot_backend="plotly")
    probe_cfdem.plot_probe(plot_backend="matplotlib")
    probe_cfdem.plot_pressure_v()
