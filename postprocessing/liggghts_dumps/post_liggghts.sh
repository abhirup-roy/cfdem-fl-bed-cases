#!/bin/bash
#SBATCH --job-name=postliggghts      # Job name to identify from queue
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=25                  # Number of processes
#SBATCH --qos bbdefault
#SBATCH --constraint=icelake
#SBATCH --time=01:00:00            # Time limit in HH:MM:SS

set -e

# Load the required modules
module purge; module load bluebear

export VENV_DIR="${HOME}/virtual-environments"
export VENV_PATH="${VENV_DIR}/my-virtual-env-${BB_CPU}"

# Create a master venv directory if necessary
mkdir -p ${VENV_DIR}

# Check if virtual environment exists and create it if not
if [[ ! -d ${VENV_PATH} ]]; then
    python3 -m venv --system-site-packages ${VENV_PATH}
fi

# Activate the virtual environment
echo "Activating virtual environment at ${VENV_PATH}!"      
source ${VENV_PATH}/bin/activate

echo "Virtual environment activated!"
# Store pip cache in /scratch directory, instead of the default home directory location
PIP_CACHE_DIR="/scratch/${USER}/pip"


# Perform any required pip installations.   
pip install --upgrade pip

echo "PIP install complete!"    

# Copy files to the post directory
cp -r ../../DEM/post/dump*.liggghts_run .


# Run the postprocessing script
: '
    usage: pizza [options] dump.example
    where dump.example is a filename or regular expression passing all relevant dump files to pizza.
    Important command line options:
        -o fname    : define output file name (default is liggghts + timestep number)
        --chunksize : sets the chunksize, default: 8
        --cpunum    : sets the number of processes to start, default (and maximum) is the amout of cpu cores avaliable at your system
        --help      : writes this help message and exits
        --no-overwrite: disables overwriting of already post-processed files.
        --timesteps: time steps to be converted, input as comma seperated list.
        --Nth: every Nth time step will be converted, cannot be combined with timesteps.
    For details, read README_GRANULAR.txt

 '

~/lpp dump*.liggghts_run --cpunum 12 

rm -r *_boundingBox.vtk
