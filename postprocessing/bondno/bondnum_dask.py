import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import dask.dataframe as dd
from dask_jobqueue import SLURMCluster

# class BondNumAnalysis:
    
#     def __init__(self, dump_dir, deltaT=5e-6):
#         self.dump_dir = dump_dir
#         self.deltaT = deltaT

#     def _read_dumpsteps(self):
#         files = os.listdir(self.dump_dir)
#         dumpsteps = np.zeros(len(files))
#         for i in range(len(files)):
#             dumpsteps[i] = int(
#                 files[i]
#                 .split('.')[0]
#                 .split('contactarea')[1]
#             )
#         self.timesteps = dumpsteps * self.deltaT

#     def _read_areas(self):
#         self._read_dumpsteps()
#         print("Making cluster")
#         cluster = SLURMCluster(
#             account="windowcr-astrazeneca-abhi",
#             cores=4, 
#             memory="50GB",
#             processes=1,
#             walltime="01:00:00"
#         )
#         cluster.scale(jobs=4)
#         print("Print Making Client")
#         client = cluster.get_client()

#         print("Starting analysis")
#         df = dd.read_csv(self.dump_dir+"*.liggghts_run", skiprows=8, include_path_column=True, sample_rows=100, low_memory=False, dtype={'ITEM: ENTRIES c_fc[1] ': 'object'})
#         col1 = df.columns[0]
#         df[col1] = dd.to_numeric(df[col1], errors="coerce")
#         area_df = (
#             df.groupby('path')
#             .mean()
#             .reset_index(drop=True)
#             .compute()
#         )
#         area_df.index = self.timesteps
#         area_df.to_csv("bondnum.csv")
        
#         client.close()
    


# if __name__ == "__main__":
#     DUMP_DIR = "../../DEM/post/bondnum/"

#     bna = BondNumAnalysis(dump_dir=DUMP_DIR)
#     bna._read_dumpsteps()

def calc_bondnum(df, k, r, rho=2700):
    """
    Calculate bond number from contact area dataframe

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing contact area values
    k : float
        Cohesion energy
    r : float
        Particle radius (m)
    rho : int
        Particle density (kg m^-3)
    """
    A = df.iloc[:,0]

    f_IMF = A * k
    f_weight = 4/3 * np.pi * r**3 * rho * 9.81

    bondnum = f_IMF / f_weight

    df["bondnum"] = bondnum


# !!! SPECIFY:
dump_dir = "../../DEM/post/bondnum/"
deltaT = 5e-6


# Get time from dump files
files = os.listdir(dump_dir)
dumpsteps = np.zeros(len(files))
for i in range(len(files)):
    dumpsteps[i] = int(
        files[i]
        .split('.')[0]
        .split('contactarea')[1]
    )
timesteps = dumpsteps * deltaT

# Start dask cluster
cluster = SLURMCluster(
    account="windowcr-astrazeneca-abhi",
    cores=4, 
    memory="50GB",
    processes=1,
    walltime="01:00:00"
)
cluster.scale(jobs=4)
client = cluster.get_client()

# Read area data in parallel
df = dd.read_csv(dump_dir+"*.liggghts_run", skiprows=8, include_path_column=True, sample_rows=100, low_memory=False, dtype={'ITEM: ENTRIES c_fc[1] ': 'object'})
col1 = df.columns[0]
df[col1] = dd.to_numeric(df[col1], errors="coerce")
# Aggregate data into average contact area per particle
area_df = (
    df.groupby('path')
    .mean()
    .reset_index(drop=True)
    .compute()
)
client.close()
area_df["timestep"] = timesteps

# Calculate bond number
calc_bondnum(area_df, 1e3, 75e-6)

# Save plot
fig = plt.figure([10,5])
plt.plot(area_df["timestep"], area_df["bondnum"])
plt.xlabel("Time (s)")
plt.ylabel("Average Bond Number")
plt.savefig("bondnum.png")