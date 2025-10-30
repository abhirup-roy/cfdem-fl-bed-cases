#!/usr/bin/env python3
# encoding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyevtk.hl import pointsToVTK


def _liggghtsdump2df(dump_filepath: str, col_names: str) -> pd.DataFrame:
    """
    Helper function to convert LIGGGHTS dump files to a pandas DataFrame.
    Args:
      dump_filepath:
        Path to the LIGGGHTS dump file.
      col_names:
        List of column names for the DataFrame.
    Returns:
      A pandas DataFrame containing the data from the LIGGGHTS dump file.
    """

    dump_df = pd.read_csv(
        dump_filepath, skiprows=9, sep=" ", header=None, engine="pyarrow"
    ).iloc[:, :-1]

    dump_df.columns = col_names

    return dump_df


def liggghts2vtk(
    timestep: float = 5e-6,
    vtk_dir: str = None,
    dump_every: int = None,
    liggghts_dump_dir: str = None,
    file_suffix: str = ".liggghts_run",
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
        raise Exception(f"LIGGGHTS dump directory {liggghts_dump_dir} does not exist")

    def is_liggghts_file(filename):
        return filename.endswith(file_suffix)

    liggghts_files = os.listdir(liggghts_dump_dir)
    liggghts_files = list(filter(is_liggghts_file, liggghts_files))

    if not liggghts_files:
        raise Exception(f"No LIGGGHTS dump files found in {liggghts_dump_dir}")

    def sort_key(filename):
        match = re.search(r"\d+", filename)
        return int(match.group(0)) if match else 0

    liggghts_files.sort(key=sort_key)

    if dump_every:
        liggghts_files = liggghts_files[::dump_every]

    with open(os.path.join(liggghts_dump_dir, liggghts_files[0]), "r") as f:
        header = f.read().split("\n")[8]
        header = header.split()[2:]

    dump_diff = sort_key(liggghts_files[1]) - sort_key(liggghts_files[0])
    step_diff = sort_key(liggghts_files[0]) - dump_diff

    for file in liggghts_files:
        vtk_name = file.split(file_suffix)[0]

        dump_df = _liggghtsdump2df(
            os.path.join(liggghts_dump_dir, file), col_names=header
        )

        x = np.array(dump_df["x"].values, dtype=np.float64)
        y = np.array(dump_df["y"].values, dtype=np.float64)
        z = np.array(dump_df["z"].values, dtype=np.float64)
        velocity = np.array(
            np.linalg.norm(dump_df[["vx", "vy", "vz"]].values, axis=1), dtype=np.float64
        )
        velocity_x = np.array(dump_df["vx"].values, dtype=np.float64)
        velocity_y = np.array(dump_df["vy"].values, dtype=np.float64)
        velocity_z = np.array(dump_df["vz"].values, dtype=np.float64)
        radius = np.array(dump_df["radius"].values, dtype=np.float64)

        time = (sort_key(file) - step_diff) * timestep

        pointsToVTK(
            os.path.join(vtk_dir, vtk_name),
            x,
            y,
            z,
            data=dict(
                time=np.full(len(x), time, dtype=np.float64),
                radius=radius,
                velocity=velocity,
                velocity_x=velocity_x,
                velocity_y=velocity_y,
                velocity_z=velocity_z,
            ),
        )


def msq_displ(
    time_rng: tuple[float, float] = None,
    dump_dir: str = "DEM/post",
    dump: bool = True,
    plot: bool = True,
    timestep: float = 5e-6,
    direction: str = None,
):
    """
    Calculate and plot the mean squared displacement of a particle or all particles in a LIGGGHTS simulation.
    Args:
        time_rng: A tuple specifying the start and end time for the analysis.
        dump_dir: Directory containing the LIGGGHTS dump files.
        dump: If True, saves the mean squared displacement data to a .npy file.
        plot: If True, generates and saves a histogram plot of the mean squared displacement.
        timestep: Time step of the simulation in seconds.
        direction: Direction for displacement calculation ('x', 'y', 'z'). If None, uses 'z' direction.
    Returns:
        A pandas Series containing the mean squared displacement for each particle.
    """
    if not os.path.isdir(dump_dir):
        raise FileNotFoundError(f"Dump directory {dump_dir} does not exist")

    dump_files = [f for f in os.listdir(dump_dir) if f.endswith(".liggghts_run")]
    if not dump_files:
        raise FileNotFoundError(f"No LIGGGHTS dump files found in {dump_dir}")

    # sort in order of time
    sort_key = lambda f: int(re.search(r"\d+", f).group(0))
    dump_files.sort(key=sort_key)

    # extract timesteps from filenames
    times = np.array([float(sort_key(f)) for f in dump_files])
    # normalise times to start from 0
    times -= times[0]
    # convert to seconds
    times *= timestep

    # get column names from the first dump file
    with open(os.path.join(dump_dir, dump_files[0]), "r") as f:
        col_names = f.read().split("\n")[8]
        col_names = col_names.split()[2:]

    if time_rng:
        if not (isinstance(time_rng, tuple) and len(time_rng) == 2):
            raise ValueError(
                "`time_rng` must be tuple with a start and end time"
                "Leave as None for all times"
            )
        # slice within time range
        start_t, stop_t = time_rng
        start_idx = np.searchsorted(times, start_t, side="left")
        stop_idx = np.searchsorted(times, stop_t, side="right")

        dump_files = dump_files[start_idx : stop_idx + 1]
        times = times[start_idx : stop_idx + 1]

    df_store = []
    for i, file in enumerate(dump_files):
        t_current = times[i]
        print(t_current)
        dump_df = _liggghtsdump2df(os.path.join(dump_dir, file), col_names=col_names)

        dump_df["time"] = t_current
        dump_df.set_index(["time", "id"], inplace=True)
        df_store.append(dump_df)

    msd_df = pd.concat(df_store, axis=0).sort_index()

    dz = msd_df.groupby("id")["z"].diff() ** 2
    msd_by_particle = dz.groupby(level="id").sum()

    if dump:
        ids = msd_by_particle.index.to_numpy()
        msd_vals = msd_by_particle.to_numpy()

        msd_2d = np.vstack((ids, msd_vals))
        np.save(os.path.join("pyoutputs", "msd.npy"), msd_2d)
    if plot:
        fig = plt.figure(figsize=(10, 10))
        hist, bins = np.histogram(msd_by_particle.values, bins=25, density=True)
        freq = hist / hist.sum()
        width = np.diff(bins)

        plt.bar(
            bins[1:],
            freq,
            width=width,
            align="edge",
            ec="k",
        )
        fig.tight_layout()
        plt.xlabel("Mean Squared Displacement (m$^2$)")
        plt.ylabel("Frequency")

        plt.savefig(os.path.join("pyoutputs", "msd_histogram.png"), bbox_inches="tight")

    return msd_by_particle
