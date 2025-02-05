#!/bin/bash
#SBATCH --job-name=calcBondNo          # Job name to identify from queue
#SBATCH --mail-type=FAIL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                 # Number of processes
#SBATCH --time=1:00:00                   # Time limit in days-hrs (can be hrs:min:sec)
#SBATCH --qos bbdefault
#SBATCH --constraint=icelake
#SBATCH --mem=32G

module purge; module load bluebear
module load bear-apps/2022b
module load Python/3.10.8-GCCcore-12.2.0
module load SciPy-bundle/2023.02-gfbf-2022b

export VENV_DIR="${HOME}/virtual-environments"
export VENV_PATH="${VENV_DIR}/my-virtual-env-${BB_CPU}"

source ${VENV_PATH}/bin/activate

export PIP_CACHE_DIR="/scratch/${USER}/pip"

python3 bondnum.py
