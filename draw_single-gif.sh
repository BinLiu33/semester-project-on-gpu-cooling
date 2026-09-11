#!/bin/bash -l

set -e

cd /users/binliu/sod_tube

uenv run --view=modules pkdgrav3/3.4 -- bash -lc '
    module load gcc pkdgrav3-python
    source .venv-plot/bin/activate
    export MPLBACKEND=Agg
    python plot_gif.py ./results sod_10_kwheat 1 100 10
'
