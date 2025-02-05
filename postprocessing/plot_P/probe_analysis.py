#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px


# probe_path = '../../CFD/postProcessing/probes/0/p'
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

# """
# Matplotlib plotting
# """
# # probe_df.plot(xlabel="Time (s)", ylabel="Pressure (Pa)", title="Pressure at Probes")
# # plt.savefig("probe_pressure.png")




"""
Plotly plotting
"""
pd.options.plotting.backend = "plotly"
times_lst = [0,0.5, 0.6, 1.6, 1.7, 2.7, 2.8, 3.8, 3.9, 4.9, 5, 6,
             6.1, 7.1, 7.2, 8.2, 8.3, 9.3, 9.4, 10.4, 10.5, 11.5,
             12.5, 12.6, 13.6, 13.7]

# prob_px = px.line(probe_df, title="Pressure at Probes")
# prob_px.update_xaxes(title_text="Time (s)")
# prob_px.update_yaxes(title_text="Pressure (Pa)")
# prob_px.update_layout(hovermode="x unified")
# prob_px.write_html("probe_pressure.html")


probe_px = probe_df.plot(title="Pressure at Probes", template="simple_white",
              labels=dict(index="Time (s)", value="Pressure", variable="Probe"))

for time in times_lst:
    probe_px.add_vline(x=time, line_width=3, line_dash="dash", line_color="green")


probe_px.write_html("probe_pressure.html")


# Look at delta pressure
delta_p1 = probe_df["Probe 1"] - probe_df["Probe 0"]
delta_p2 = probe_df["Probe 2"] - probe_df["Probe 0"]
delta_p3 = probe_df["Probe 3"] - probe_df["Probe 0"]
delta_p4 = probe_df["Probe 4"] - probe_df["Probe 0"]

df_deltaP = pd.concat([delta_p1, delta_p2, delta_p3, delta_p4], axis=1)

# prob_px = px.line(df_deltaP, title="Pressure Difference Probes")
# prob_px.update_xaxes(title_text="Time (s)")
# prob_px.update_yaxes(title_text="Pressure Difference(Pa)")
# prob_px.update_layout(hovermode="x unified")
# prob_px.write_html("probe_pressure_delta.html")

deltaP_px = df_deltaP.plot(title="Pressure Difference Probes", template="simple_white",
                labels=dict(index="Time (s)", value="Pressure Difference (Pa)", variable="Probe"))

for time in times_lst:
    deltaP_px.add_vline(x=time, line_width=3, line_dash="dash", line_color="green")


deltaP_px.write_html("probe_pressure_delta.html")
