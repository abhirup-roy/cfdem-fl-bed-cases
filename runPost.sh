module purge; module load bluebear
module load bear-apps/2022b
module load Python/3.10.8-GCCcore-12.2.0
module load SciPy-bundle/2023.02-gfbf-2022b
export VENV_DIR="${HOME}/virtual-environments"
export VENV_PATH="${VENV_DIR}/my-virtual-env-${BB_CPU}"

source ${VENV_PATH}/bin/activate

export PIP_CACHE_DIR="/scratch/${USER}/pip"

current_dir=$PWD

echo "Analysing & plotting pressures"
cd postprocessing/plot_P

python3 post_scripts.py

echo "Probe analysis complete"
cd $current_dir

deactivate