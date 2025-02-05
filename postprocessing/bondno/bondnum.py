#!/usr/bin/env python3
import os

import numpy as np
import matplotlib.pyplot as plt


class BondNumAnalysis:

    def __init__(self, area_dir, run_script_path, init_script_path):
        self.area_dir = area_dir
        self.run_script_dir = run_script_path
        self.init_script_dir = init_script_path
        self._read_scripts()
        

    def _dump2area(self, file_path):

        with open(file_path, 'r') as file:
            # read file
            values = file.read().splitlines()[9:]
            value_arr = np.array(values, dtype=float)
        return value_arr

    def _read_dir(self):

        files = os.listdir(self.area_dir) 

        dump_nums = []
        areas = []

        for file in files:
            # extract dump number from file name
            timestep_num = file.split('.')[0].split('contactarea')[1]
            dump_nums.append(int(timestep_num))

            # add area values to list
            areas.append(self._dump2area(self.area_dir + file))
        
        return np.array(dump_nums), areas
    
    def _read_scripts(self):
        with open(self.run_script_dir, 'r') as file:
            script = file.read().splitlines()
            k_line = script[26]
            k = float(k_line.split()[-1])

            timestep_line = script[33]
            timestep = float(timestep_line.split()[-1]) 
                
            self.k = k
            self.timestep = timestep

        with open(self.init_script_dir, 'r') as file:
            script = file.read().splitlines()

            r_line = script[81]
            r = float(r_line.split()[-1])
            rho_line = script[48]
            rho = float(rho_line.split()[-1])
                
            self.r = r
            self.rho = rho
    
    def avg_bondnum(self, plot = True):
        dump_nums, areas = self._read_dir()

        avg_areas = np.array([area_arr.mean() for area_arr in areas])
        times = dump_nums * self.timestep

        time_idx = np.argsort(times)
        times = times[time_idx]
        avg_areas = avg_areas[time_idx]

        f_IMF = avg_areas * self.k
        f_weight = 4/3 * np.pi * self.r**3 * self.rho * 9.81

        avg_bond_nums = f_IMF / f_weight
        
        if plot:
            plt.plot(times, avg_bond_nums)
            plt.xlabel('Time (s)')
            plt.ylabel('Average Bond Number')
            plt.title('Average Bond Number vs Time')
            plt.savefig('bondnum.png')

        return times, avg_areas


if __name__ == '__main__':

    timestep = 0.000005
    area_dir = '../../DEM/post/bondnum/'
    run_script_path = '../../DEM/in.liggghts_run'
    init_script_path = '../../DEM/in.liggghts_init'

    analysis = BondNumAnalysis(
        area_dir=area_dir, 
        run_script_path=run_script_path,
        init_script_path=init_script_path
    )
    analysis.avg_bondnum(plot=True)

