#!/bin/bash -l
set -euo pipefail

uenv run --view=modules prgenv-gnu -- bash -lc '
  module load boost cmake cray-mpich fftw gcc gsl hdf5 python cuda

  cd ~/pkdgrav3
  rm -rf build
  mkdir build
  cd build

  python -m venv ../.venv
  source ../.venv/bin/activate

  python -m pip install -r ../requirements.txt

  cmake ..
  make -j 40
'
