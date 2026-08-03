#!/bin/bash -l
#   ./sod_tube/draw_compare_gif.sh                         # default sod_16000 1 100 14400
#   ./sod_tube/draw_compare_gif.sh sod_10 1 100 9          # diy
#   parameters <file_keyword> <step_interval> <nsteps> <temperature>
set -euo pipefail

SOD_DIR=/users/binliu/sod_tube
VENV="$SOD_DIR/.venv-plot"
UENV_IMAGE=pkdgrav3/3.4

KEYWORD="${1:-sod_16000}"
STEP_INTERVAL="${2:-1}"
NSTEPS="${3:-100}"
TEMP="${4:-14400}"

do_draw() {
    module load gcc pkdgrav3-python >/dev/null 2>&1
    # shellcheck disable=SC1090
    source "$VENV/bin/activate"
    export MPLBACKEND=Agg
    cd "$SOD_DIR"
    echo "== plot_compare3.py ./results $KEYWORD $STEP_INTERVAL $NSTEPS $TEMP =="
    python plot_compare3.py ./results "$KEYWORD" "$STEP_INTERVAL" "$NSTEPS" "$TEMP"
    echo "== Done. Output: $SOD_DIR/results/${KEYWORD}_compare3_{density,pressure,temperature}.gif =="
}

if [ -n "${UENV_MOUNT_LIST:-}" ]; then
    do_draw
else
    export SOD_DIR VENV KEYWORD STEP_INTERVAL NSTEPS TEMP
    uenv run --view=modules "$UENV_IMAGE" -- bash -lc '
        set -euo pipefail
        module load gcc pkdgrav3-python >/dev/null 2>&1
        source "$VENV/bin/activate"
        export MPLBACKEND=Agg
        cd "$SOD_DIR"
        echo "== plot_compare3.py ./results $KEYWORD $STEP_INTERVAL $NSTEPS $TEMP =="
        python plot_compare3.py ./results "$KEYWORD" "$STEP_INTERVAL" "$NSTEPS" "$TEMP"
        echo "== Done. Output: $SOD_DIR/results/${KEYWORD}_compare3_{density,pressure,temperature}.gif =="
    '
fi
