module purge; module load bluebear
module load bear-apps/2022b
module load Python/3.10.8-GCCcore-12.2.0
module load SciPy-bundle/2023.02-gfbf-2022b
export VENV_DIR="${HOME}/virtual-environments"
export VENV_PATH="${VENV_DIR}/my-virtual-env-${BB_CPU}"

source ${VENV_PATH}/bin/activate

export PIP_CACHE_DIR="/scratch/${USER}/pip"

echo "Analysing & plotting pressures"

python3 plot_fluidn_curves.py

echo "Probe analysis complete"

deactivate