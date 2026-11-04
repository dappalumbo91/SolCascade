#!/bin/bash
# Run SolCascade inside the WSL conda environment that has OpenMC 0.16.
# Nuclear data is the official ENDF/B-VIII.0 HDF5 library.
set -euo pipefail
source /opt/miniforge3/etc/profile.d/conda.sh
conda activate solcascade
export PYTHONNOUSERSITE=1
export OPENMC_CROSS_SECTIONS=/opt/nuclear-data/endfb-viii.0-hdf5/cross_sections.xml
cd "/mnt/c/Users/damia/Desktop/SolCascade fsot"
exec python -m solcascade "$@"
