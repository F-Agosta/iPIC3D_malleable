#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -A eupadmin
#SBATCH --mem=0
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --job-name=ipic3d-restart-repartition
#SBATCH -o restart-repartition-%j.out
#SBATCH -e restart-repartition-%j.err

module purge
module load OpenMPI/5.0.x

export UCX_IB_MLX5_DEVX=n
export LD_LIBRARY_PATH=/home/francesco.agosta/applications/hdf5/lib:${LD_LIBRARY_PATH:-}

set -euo pipefail

TEST_MODE="${TEST_MODE:-all}"
case "$TEST_MODE" in
    baseline|control|repartition|repartition-4-to-8|all)
        ;;
    *)
        echo "Invalid TEST_MODE: $TEST_MODE" >&2
        echo "Usage: TEST_MODE=baseline|control|repartition|repartition-4-to-8|all" >&2
        exit 2
        ;;
esac

ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
EXE="$ROOT/build-hdf5-ompi51/iPIC3D"
TEST_ROOT="$ROOT/restart-repartition-test"

CHECKPOINT_BASELINE="$TEST_ROOT/runs/restart/baseline-2"
CHECKPOINT_CONTROL="$TEST_ROOT/runs/restart/restart-same-2"
CHECKPOINT_2_TO_4="$TEST_ROOT/runs/restart/restart-2-to-4"
CHECKPOINT_4_TO_8="$TEST_ROOT/runs/restart/restart-4-to-8"
CHECKPOINT_8_TO_4="$TEST_ROOT/runs/restart/restart-4-to-8"


BASELINE_OUTPUT="$TEST_ROOT/runs/save/baseline-2"
CONTROL_OUTPUT="$TEST_ROOT/runs/save/restart-same-2"
REPARTITION_OUTPUT="$TEST_ROOT/runs/save/restart-2-to-4"
REPARTITION_4_TO_8_OUTPUT="$TEST_ROOT/runs/save/restart-4-to-8"
REPARTITION_8_TO_4_OUTPUT="$TEST_ROOT/runs/save/restart-8-to-4"

BASELINE_INPUT="$TEST_ROOT/inputs/baseline-2.inp"
CONTROL_INPUT="$TEST_ROOT/inputs/restart-same-2.inp"
REPARTITION_INPUT="$TEST_ROOT/inputs/restart-2-to-4.inp"
REPARTITION_4_TO_8_INPUT="$TEST_ROOT/inputs/restart-4-to-8.inp"
REPARTITION_8_TO_4_INPUT="$TEST_ROOT/inputs/restart-2-to-4.inp"

if [[ ! -x "$EXE" ]]; then
    echo "Executable not found: $EXE" >&2
    exit 1
fi

for input in "$BASELINE_INPUT" "$CONTROL_INPUT" "$REPARTITION_INPUT" "$REPARTITION_4_TO_8_INPUT"; do
    if [[ ! -f "$input" ]]; then
        echo "Input file not found: $input" >&2
        exit 1
    fi
done

verify_checkpoint_baseline() {
    for file in settings.hdf restart_x000_y000_z000.hdf restart_x001_y000_z000.hdf; do
        if [[ ! -f "$CHECKPOINT_BASELINE/$file" ]]; then
            echo "Expected baseline checkpoint file not found: $CHECKPOINT_BASELINE/$file" >&2
            exit 1
        fi
    done
}

verify_checkpoint_2_to_4() {
    for file in settings.hdf \
                restart_x000_y000_z000.hdf restart_x001_y000_z000.hdf \
                restart_x002_y000_z000.hdf restart_x003_y000_z000.hdf; do
        if [[ ! -f "$CHECKPOINT_2_TO_4/$file" ]]; then
            echo "Expected 2-to-4 checkpoint file not found: $CHECKPOINT_2_TO_4/$file" >&2
            exit 1
        fi
    done
}

verify_checkpoint_4_to_8() {
    for file in settings.hdf \
                restart_x000_y000_z000.hdf restart_x001_y000_z000.hdf \
                restart_x002_y000_z000.hdf restart_x003_y000_z000.hdf \
                restart_x004_y000_z000.hdf restart_x005_y000_z000.hdf \
                restart_x006_y000_z000.hdf restart_x007_y000_z000.hdf; do
        if [[ ! -f "$CHECKPOINT_4_TO_8/$file" ]]; then
            echo "Expected 4-to-8 checkpoint file not found: $CHECKPOINT_4_TO_8/$file" >&2
            exit 1
        fi
    done
}

verify_checkpoint_8_to_4() {
    for file in settings.hdf \
                restart_x000_y000_z000.hdf restart_x001_y000_z000.hdf \
                restart_x002_y000_z000.hdf restart_x003_y000_z000.hdf; do
        if [[ ! -f "$CHECKPOINT_8_TO_4/$file" ]]; then
            echo "Expected 8-to-4 checkpoint file not found: $CHECKPOINT_8_TO_4/$file" >&2
            exit 1
        fi
    done
}

run_baseline() {
    mkdir -p "$CHECKPOINT_BASELINE" "$BASELINE_OUTPUT"

    if [[ -n "$(find "$CHECKPOINT_BASELINE" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
        echo "Checkpoint directory is not empty: $CHECKPOINT_BASELINE" >&2
        echo "Move or clean the existing checkpoint before starting a new test." >&2
        exit 1
    fi

    echo "Generating the 2-process source checkpoint"
    (
        cd "$ROOT"
        mpirun -np 2 "$EXE" "$BASELINE_INPUT"
    )

    verify_checkpoint_baseline
}

run_control() {
    verify_checkpoint_baseline
    mkdir -p "$CONTROL_OUTPUT" "$CHECKPOINT_CONTROL"

    echo "Running the unchanged-topology 2-process control restart"
    (
        cd "$ROOT"
        mpirun -np 2 "$EXE" restart "$CONTROL_INPUT"
    )
}

run_repartition() {
    verify_checkpoint_baseline
    mkdir -p "$REPARTITION_OUTPUT" "$CHECKPOINT_2_TO_4"

    echo "Running the 2-to-4 restart"
    (
        cd "$ROOT"
        mpirun -np 4 "$EXE" restart "$REPARTITION_INPUT"
    )

    verify_checkpoint_2_to_4
}

run_repartition_4_to_8() {
    verify_checkpoint_2_to_4
    mkdir -p "$REPARTITION_4_TO_8_OUTPUT" "$CHECKPOINT_4_TO_8"

    echo "Running the 4-to-8 restart"
    (
        cd "$ROOT"
        mpirun -np 8 "$EXE" restart "$REPARTITION_4_TO_8_INPUT"
    )

    verify_checkpoint_4_to_8
}

run_repartition_8_to_4() {
    verify_checkpoint_4_to_8
    mkdir -p "$REPARTITION_8_TO_4_OUTPUT" "$CHECKPOINT_8_TO_4"

    echo "Running the 8-to-4 restart"
    (
        cd "$ROOT"
        mpirun -np 4 "$EXE" restart "$REPARTITION_INPUT"
    )

    verify_checkpoint_8_to_4
}

case "$TEST_MODE" in
    baseline)
        run_baseline
        ;;
    control)
        run_control
        ;;
    repartition)
        run_repartition
        ;;
    repartition-4-to-8)
        run_repartition_4_to_8
        ;;
    repartition-8-to-4)
        run_repartition_8_to_4
        ;;
    all)
        run_baseline
        run_control
        run_repartition
        run_repartition_4_to_8
        run_repartition_8_to_4
        ;;
esac

echo "Restart test completed in mode: $TEST_MODE"
echo "Run directory: $TEST_ROOT/runs"
