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

# Restart-repartition test with correctness verification.
#
# Usage:
#   sbatch test_restart_repartition.sh                                  # CASE=weibel, all stages + verify
#   TEST_MODE=baseline CASE=twostream sbatch test_restart_repartition.sh # single stage
#   TEST_MODE=verify   CASE=weibel   sbatch test_restart_repartition.sh # re-run only the checks
#
# CASE:      minimal | weibel | twostream   (master input in inputs/<case>.inp,
#            derived from the basic test input, inputfiles/Weibel256.inp and
#            inputfiles/two_stream.inp respectively)
# TEST_MODE: baseline | control | repartition | control-4 | repartition-4-to-8 |
#            control-8 | repartition-8-to-4 | verify | all
# TOL:       max relative error accepted by the verification (default 1e-8)
#
# Stages (all under runs/<case>/):
#   baseline-2        fresh run, 2 procs (XLEN=2), ncycles=BASE_CYCLES, checkpoint at cycle RESTART_CYCLE
#   restart-same-2    restart baseline-2 -> 2 procs  (control for the restart itself)
#   restart-2-to-4    restart baseline-2 -> 4 procs  (growing)
#   restart-same-4    restart restart-2-to-4 -> 4 procs (control)
#   restart-4-to-8    restart restart-2-to-4 -> 8 procs (growing)
#   restart-same-8    restart restart-4-to-8 -> 8 procs (control)
#   restart-8-to-4    restart restart-4-to-8 -> 4 procs (shrinking)
# Every stage runs STAGE_CYCLES cycles starting right after the source
# checkpoint, so each repartition stage and its same-topology control cover
# the SAME global cycles and can be compared row-by-row / field-by-field.
#
# verify compares, for each (control, repartition) pair and (baseline, control):
#   - ConservedQuantities.txt at common cycles  (scripts/plot_conserved.py --compare)
#   - reassembled fields at the last common cycle (scripts/read_hdf.py --compare)

module purge
module load OpenMPI/5.0.x

export UCX_IB_MLX5_DEVX=n
export LD_LIBRARY_PATH=/home/francesco.agosta/applications/hdf5/lib:${LD_LIBRARY_PATH:-}

set -euo pipefail

CASE="${CASE:-weibel}"
TEST_MODE="${TEST_MODE:-all}"
TOL="${TOL:-1e-8}"
PYTHON="${PYTHON:-python3}"

case "$CASE" in
    minimal|weibel|twostream) ;;
    *) echo "Invalid CASE: $CASE (use minimal|weibel|twostream)" >&2; exit 2 ;;
esac

case "$TEST_MODE" in
    baseline|control|repartition|control-4|repartition-4-to-8|control-8|repartition-8-to-4|verify|all) ;;
    *) echo "Invalid TEST_MODE: $TEST_MODE" >&2
       echo "Usage: TEST_MODE=baseline|control|repartition|control-4|repartition-4-to-8|control-8|repartition-8-to-4|verify|all" >&2
       exit 2 ;;
esac

# fields compared by read_hdf.py in verify mode (per case)
case "$CASE" in
    weibel)         FIELDS="Ez Bz" ;;   # 2D Weibel: Bz grows, Ez electrostatic
    minimal|twostream) FIELDS="Ez" ;;   # 1D: electrostatic field along x
esac

# stage lengths / output cadence (cycles)
BASE_CYCLES=40        # baseline run length
STAGE_CYCLES=20       # length of each restart stage
RESTART_CYCLE=20      # checkpoint cadence (baseline writes one checkpoint at cycle 20)
DIAG_CYCLE=5          # DiagnosticsOutputCycle  -> rows in ConservedQuantities.txt
FIELD_CYCLE=5         # FieldOutputCycle       -> /fields/*/cycle_NNNNNN in proc*.hdf
PART_CYCLE=0          # ParticlesOutputCycle   -> particle dump off (restart files carry particles)

ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
EXE="$ROOT/build-hdf5-ompi51/iPIC3D"
TEST_ROOT="$ROOT/restart-repartition-test"
MASTER_INPUT="$TEST_ROOT/inputs/$CASE.inp"

RUNS="$TEST_ROOT/runs/$CASE"
GEN_INPUT_DIR="$RUNS/inputs"
VERIFY_DIR="$RUNS/verify"

if [[ ! -x "$EXE" ]]; then
    echo "Executable not found: $EXE" >&2
    exit 1
fi
if [[ ! -f "$MASTER_INPUT" ]]; then
    echo "Master input not found: $MASTER_INPUT" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# set_key <file> <key> <value> : replace "key = ..." line or append it
set_key() {
    local f=$1 k=$2 v=$3
    if grep -qE "^[[:space:]]*${k}[[:space:]]*=" "$f"; then
        sed -i.bak -E "s|^[[:space:]]*${k}[[:space:]]*=.*|${k} = ${v}|" "$f"
        rm -f "$f.bak"
    else
        printf '%s = %s\n' "$k" "$v" >> "$f"
    fi
}

# make_stage_input <stage> <xlen> <source_stage|"">
# generates $GEN_INPUT_DIR/<stage>.inp from the case master
make_stage_input() {
    local stage=$1 xlen=$2 source=$3
    local out="$GEN_INPUT_DIR/$stage.inp"
    mkdir -p "$GEN_INPUT_DIR"
    cp "$MASTER_INPUT" "$out"

    set_key "$out" SaveDirName        "$RUNS/save/$stage"
    set_key "$out" RestartSaveDirName "$RUNS/restart/$stage"
    set_key "$out" XLEN "$xlen"
    set_key "$out" YLEN 1
    set_key "$out" ZLEN 1
    set_key "$out" WriteMethod shdf5
    set_key "$out" FieldOutputTag     "Eall + Ball + Jall + rhos + pressure"
    set_key "$out" MomentsOutputTag   "Eall + Ball + Jall + rhos + pressure"
    set_key "$out" ParticlesOutputTag "position + velocity + q"
    set_key "$out" FieldOutputCycle       "$FIELD_CYCLE"
    set_key "$out" ParticlesOutputCycle   "$PART_CYCLE"
    set_key "$out" DiagnosticsOutputCycle "$DIAG_CYCLE"

    if [[ -z "$source" ]]; then
        # baseline: fresh run, single mid-run checkpoint, no finalize checkpoint
        set_key "$out" RestartDirName     "$RUNS/restart/$stage"
        set_key "$out" ncycles            "$BASE_CYCLES"
        set_key "$out" RestartOutputCycle "$RESTART_CYCLE"
        set_key "$out" CallFinalize       0
    else
        # restart stage: read the source checkpoint, write our own at the last cycle
        set_key "$out" RestartDirName     "$RUNS/restart/$source"
        set_key "$out" ncycles            "$STAGE_CYCLES"
        set_key "$out" RestartOutputCycle "$STAGE_CYCLES"
        set_key "$out" CallFinalize       1
    fi
    echo "$out"
}

# expected checkpoint files for a stage with <nprocs> processes in x
checkpoint_files() {
    local nprocs=$1 i
    echo "settings.hdf"
    for ((i = 0; i < nprocs; i++)); do
        printf 'restart_x%03d_y000_z000.hdf\n' "$i"
    done
}

verify_checkpoint() {  # <stage> <nprocs>
    local stage=$1 nprocs=$2 f
    while IFS= read -r f; do
        if [[ ! -f "$RUNS/restart/$stage/$f" ]]; then
            echo "Expected checkpoint file not found: $RUNS/restart/$stage/$f" >&2
            exit 1
        fi
    done < <(checkpoint_files "$nprocs")
}

# run_stage <stage> <nprocs> <xlen> <source_stage|"">
run_stage() {
    local stage=$1 nprocs=$2 xlen=$3 source=$4
    local input
    input=$(make_stage_input "$stage" "$xlen" "$source")

    mkdir -p "$RUNS/save/$stage" "$RUNS/restart/$stage"

    if [[ -n "$source" ]]; then
        # the source checkpoint must exist before restarting from it
        if [[ ! -f "$RUNS/restart/$source/settings.hdf" ]] || \
           [[ ! -f "$RUNS/restart/$source/restart_x000_y000_z000.hdf" ]]; then
            echo "Source checkpoint incomplete: $RUNS/restart/$source (run the previous stage first)" >&2
            exit 1
        fi
    else
        if [[ -n "$(find "$RUNS/restart/$stage" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
            echo "Checkpoint directory is not empty: $RUNS/restart/$stage" >&2
            echo "Move or clean the existing checkpoint before starting a new test." >&2
            exit 1
        fi
    fi

    echo "=== stage $stage : np=$nprocs XLEN=$xlen source=${source:-<fresh>} ==="
    if [[ -n "$source" ]]; then
        (cd "$ROOT" && mpirun -np "$nprocs" "$EXE" restart "$input")
    else
        (cd "$ROOT" && mpirun -np "$nprocs" "$EXE" "$input")
    fi

    verify_checkpoint "$stage" "$nprocs"
}

# last cycle recorded in a run's ConservedQuantities.txt
last_cycle_of() {  # <stage>
    tail -n 1 "$RUNS/save/$1/ConservedQuantities.txt" | cut -f1
}

# ---------------------------------------------------------------------------
# verification
# ---------------------------------------------------------------------------

# compare a (reference, test) pair of run directories:
#   conserved quantities at common cycles + field maps at the last common cycle
verify_pair() {  # <ref_stage> <test_stage>
    local ref=$1 test=$2
    local ref_file="$RUNS/save/$ref/ConservedQuantities.txt"
    local test_file="$RUNS/save/$test/ConservedQuantities.txt"
    local cycle field pair_name

    pair_name="${ref}__vs__${test}"
    echo "=== verify pair: $ref (ref) vs $test (test) ==="

    if [[ ! -f "$ref_file" || ! -f "$test_file" ]]; then
        echo "Missing ConservedQuantities.txt for pair $pair_name" >&2
        return 1
    fi

    # conserved quantities: total energy, momentum, E, B, K at common cycles
    # (columns 2,3,4,5,6 in the whitespace-collapsed format of plot_conserved.py)
    "$PYTHON" "$ROOT/scripts/plot_conserved.py" \
        --compare "$RUNS/save/$ref" "$RUNS/save/$test" \
        --columns 2,3,4,5,6 --tol "$TOL" --yscale symlog \
        --out "$VERIFY_DIR/$pair_name.pdf" || return 1

    # field maps at the last cycle present in both runs
    cycle=$(last_cycle_of "$test")
    if [[ "$(last_cycle_of "$ref")" -lt "$cycle" ]]; then
        cycle=$(last_cycle_of "$ref")
    fi
    echo "--- field comparison at cycle $cycle ---"
    for field in $FIELDS; do
        "$PYTHON" "$ROOT/scripts/read_hdf.py" \
            --directory "$RUNS/save/$ref" \
            --compare "$RUNS/save/$test" \
            --cycle "$cycle" --group fields --comp "$field" \
            --tol "$TOL" --out-prefix "$VERIFY_DIR/${pair_name}_${field}" || return 1
    done
}

run_verify() {
    mkdir -p "$VERIFY_DIR"
    local rc=0
    # restart fidelity: uninterrupted baseline vs same-topology restart
    verify_pair baseline-2 restart-same-2 || rc=1
    # growing 2->4, control 2 procs
    verify_pair restart-same-2 restart-2-to-4 || rc=1
    # growing 4->8, control 4 procs
    verify_pair restart-same-4 restart-4-to-8 || rc=1
    # shrinking 8->4, control 8 procs
    verify_pair restart-same-8 restart-8-to-4 || rc=1

    if [[ $rc -ne 0 ]]; then
        echo "VERIFY: FAILED (tol=$TOL) - see $VERIFY_DIR" >&2
        exit 1
    fi
    echo "VERIFY: all pairs agree within tol=$TOL - plots in $VERIFY_DIR"
}

# ---------------------------------------------------------------------------
# stages
# ---------------------------------------------------------------------------

run_baseline()          { run_stage baseline-2      2 2 ""; }
run_control()           { run_stage restart-same-2  2 2 baseline-2; }
run_repartition()       { run_stage restart-2-to-4  4 4 baseline-2; }
run_control_4()         { run_stage restart-same-4  4 4 restart-2-to-4; }
run_repartition_4_to_8(){ run_stage restart-4-to-8  8 8 restart-2-to-4; }
run_control_8()         { run_stage restart-same-8  8 8 restart-4-to-8; }
run_repartition_8_to_4(){ run_stage restart-8-to-4  4 4 restart-4-to-8; }

case "$TEST_MODE" in
    baseline)              run_baseline ;;
    control)               run_control ;;
    repartition)           run_repartition ;;
    control-4)             run_control_4 ;;
    repartition-4-to-8)    run_repartition_4_to_8 ;;
    control-8)             run_control_8 ;;
    repartition-8-to-4)    run_repartition_8_to_4 ;;
    verify)                run_verify ;;
    all)
        run_baseline
        run_control
        run_repartition
        run_control_4
        run_repartition_4_to_8
        run_control_8
        run_repartition_8_to_4
        run_verify
        ;;
esac

echo "Restart test completed in mode: $TEST_MODE (case: $CASE)"
echo "Run directory: $RUNS"
