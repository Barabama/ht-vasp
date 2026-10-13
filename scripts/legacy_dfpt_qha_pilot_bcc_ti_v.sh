#!/usr/bin/env bash
#SBATCH --job-name=BCC-Ti-V-legacy-dfpt
#SBATCH --output=/nfs_hdd/2025/gaominliang/ht-vasp/logs/%x-%j.out
#SBATCH --error=/nfs_hdd/2025/gaominliang/ht-vasp/logs/%x-%j.err
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=partCPU
#SBATCH --time=100:00:00

# Separate B5070 legacy-reproduction and auto-kpoint comparison recipes.
# Default is read-only review. This file never calls sbatch.
#
# Fail-closed contract
# --------------------
# 1. A stage is only reused when `legacy_dfpt_qha_pilot_bcc_ti_v.py verify`
#    confirms *electronic convergence* for that stage. "OUTCAR has the VASP
#    completion banner" is not enough: VASP 53681 finished cleanly while its
#    ground-state SCF never reached EDIFF.
# 2. A DFPT directory whose VASP outputs exist but do not verify is renamed to
#    DFPT.failed-<jobid>-<timestamp>; the evidence is never deleted, never
#    overwritten, and never handed to phonopy.
# 3. FORCE_CONSTANTS and thermal properties are only produced downstream of a
#    verified DFPT vasprun.xml.
# 4. R3 and R3f-static outputs that fail verification stop the run in place
#    (no silent recompute of results that already look finished).
set -euo pipefail
umask 027

# Slurm executes a spooled copy of this batch file, so BASH_SOURCE[0] is not
# a reliable way to discover the original repository. Use the fixed checkout.
HTVASP_ROOT="${HTVASP_ROOT:-/nfs_hdd/2025/gaominliang/ht-vasp}"
SCRIPT_DIR="${HTVASP_ROOT}/scripts"
PYTHON="${HTVASP_ROOT}/.conda/bin/python"
PHONOPY="${HTVASP_ROOT}/.conda/bin/phonopy"
PHONOPY_QHA="${HTVASP_ROOT}/.conda/bin/phonopy-qha"
PREP="${SCRIPT_DIR}/legacy_dfpt_qha_pilot_bcc_ti_v.py"
PMG_VASP_PSP_DIR="${PMG_VASP_PSP_DIR:-/nfs_ssd/documents}"
I_MPI_FABRICS="${I_MPI_FABRICS:-ofi}"
PROFILE="${2:-legacy-fixed-grid}"
RUN_ID="${RUN_ID:-pilot-$(date -u +%Y%m%dT%H%M%SZ)}"
RESUME="${RESUME:-NO}"
PILOT_BASE="${PILOT_ROOT:-${HTVASP_ROOT}/data/legacy-dfpt-qha-pilot/BCC-Ti-V}"
RUN_ROOT="${PILOT_BASE}/${PROFILE}/runs/${RUN_ID}"
fail() { echo "ERROR: $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || fail "missing/empty file: $1"; }

# Volume grid, as lattice-constant scale factors (a/a0), where a0 is the
# R3-relaxed 2-site cell (a0 = 3.120007 A, V0 = 30.371528 A^3).
#
# Nine points, uniform in VOLUME (not in length): V/V0 = 0.95 .. 1.13, step
# 0.0225. Rationale, from /tmp/BCC and /tmp/FCC (B5070's reference drivers):
#
#  * Uniform in V is the convention of the legacy method: its `phonopy-qha --cu`
#    generates ten points at V/V0 = linspace(0.95, 1.10, 10). It also suits the
#    fit, which regresses F(V,T) onto a Vinet EOS with V as the abscissa.
#  * The top must reach V/V0 ~ 1.09-1.12 so the equilibrium volume stays inside
#    the grid up to TMAX=3000 K. With our Cv(T) and alpha_V = gamma*Cv/(B*V)
#    (gamma=2, B=140-160 GPa) the expansion at 3000 K is dV/V0 = 0.08-0.10.
#    The legacy 8-point grid stops at 1.0667 and only supports ~1500-2000 K;
#    its own --tmax=3000 run goes numerically unstable above ~2900 K.
#  * The bottom is held at 0.95. 0.885 (the previous tag 01) is ~+20 GPa of
#    hydrostatic compression, which V(T) never visits at ambient pressure, and
#    it is the volume with the largest imaginary-mode fraction. The legacy
#    convention and atomate2's own EOS grid both stop at 0.95.
#  * Nine points: phonopy-qha fits each temperature to a 4-parameter EOS, so
#    eight suffice; nine keeps margin without widening the range.
SCALES=(0.9830476 0.9907480 0.9983305 1.0057996 1.0131594 1.0204138 1.0275664 1.0346209 1.0415804)
TAGS=(01 02 03 04 05 06 07 08 09)
# Numeric tags to run in this submission. Default: all nine volume points.
# Use ONLY_TAGS=01 to smoke-test a single volume point before committing to the
# full QHA set. Stages outside this list are neither run nor archived.
ONLY_TAGS_RAW="${ONLY_TAGS:-${TAGS[*]}}"
read -r -a ONLY_TAGS <<< "$ONLY_TAGS_RAW"

# DFPT electronic-structure experiment. "legacy" reproduces B5070's
# INCAR-Phonopy exactly; the other variants change ALGO and/or the PAW tags.
# Each variant gets its OWN DFPT/charge-seed directories (suffix below) so
# several can run in parallel inside one RUN_ID without touching each other or
# the canonical, un-suffixed pair. Only the canonical "legacy" variant feeds
# phonopy-qha.
# "production" is ALGO=Normal + LASPH=True + LMAXMIX=4 and is the default. It
# owns the canonical (un-suffixed) DFPT/ and charge-seed/ directories.
DFPT_VARIANT="${DFPT_VARIANT:-production}"
case "$DFPT_VARIANT" in
    production)       VARIANT_SUFFIX="" ;;
    algo-normal)      VARIANT_SUFFIX="-algo-normal" ;;
    legacy)           VARIANT_SUFFIX="-legacy" ;;
    *) fail "unknown DFPT_VARIANT: $DFPT_VARIANT (production|algo-normal|legacy)" ;;
esac
DFPT_NAME="DFPT${VARIANT_SUFFIX}"
SEED_NAME="charge-seed${VARIANT_SUFFIX}"

tag_selected() {
    local needle="$1" candidate
    for candidate in "${ONLY_TAGS[@]}"; do
        [[ "$candidate" == "$needle" ]] && return 0
    done
    return 1
}

run_vasp() {
    local workdir="$1"
    echo "[VASP] $workdir"
    # Fail closed if the directory VASP is about to overwrite still holds the
    # outputs of an earlier run. Job 53702 walked into exactly that: a failed
    # archive left the 53697 seed in place and a fresh SCF then overwrote its
    # verified 32 MB CHGCAR. Nothing here is recoverable once VASP starts.
    if vasp_ran_here "$workdir"; then
        echo "ERROR: refusing to start VASP in a directory that already holds VASP output: $workdir" >&2
        echo "       archive or move the existing results first (they were not deleted)." >&2
        return 2
    fi
    (cd "$workdir" && srun --mpi=pmi2 vasp_std > vasp.log 2>&1)
}

# Write an auditable convergence report next to the outputs and return the
# verifier's status (0 = verified, 3 = failed checks, 2 = other). Callers decide
# whether to archive first; this never aborts the script on its own.
# See the fail-closed contract above.
verify_vasp() {
    local workdir="$1"
    local stage="$2"
    shift 2
    "$PYTHON" "$PREP" verify --stage "$stage" \
        --vasprun "${workdir}/vasprun.xml" --outcar "${workdir}/OUTCAR" \
        --report "${workdir}/vasp-verification.json" "$@"
}

# Same check, read-only: used by the resume plan so `review` writes nothing.
verify_quiet() {
    local workdir="$1"
    local stage="$2"
    shift 2
    [[ -s "${workdir}/vasprun.xml" && -s "${workdir}/OUTCAR" ]] || return 1
    "$PYTHON" "$PREP" verify --stage "$stage" \
        --vasprun "${workdir}/vasprun.xml" --outcar "${workdir}/OUTCAR" \
        "$@" > /dev/null 2>&1
}

vasp_complete() {
    local workdir="$1"
    [[ -s "${workdir}/vasprun.xml" && -s "${workdir}/OUTCAR" ]] \
        && grep -q "General timing and accounting informations for this job" "${workdir}/OUTCAR"
}

dfpt_has_vasp_artifacts() {
    local dfpt_root="$1"
    local file
    for file in OUTCAR vasprun.xml OSZICAR vasp.log CONTCAR DOSCAR EIGENVAL \
                IBZKPT WAVECAR CHGCAR CHG PCDAT REPORT XDATCAR; do
        [[ ! -e "$dfpt_root/$file" ]] || return 0
    done
    return 1
}

# Files that prove VASP actually ran. CHGCAR is excluded because the driver now
# places a charge density in the DFPT directory on purpose, and that must not
# be mistaken for leftover output of a previous attempt.
vasp_ran_here() {
    local root="$1"
    local file
    for file in OUTCAR vasprun.xml OSZICAR vasp.log CONTCAR DOSCAR EIGENVAL \
                IBZKPT WAVECAR CHG PCDAT REPORT XDATCAR DYNMAT; do
        [[ ! -e "$root/$file" ]] || return 0
    done
    return 1
}

# Preserve a whole unusable stage attempt (unconverged, killed, or partial) for
# audit. Never deletes or overwrites: the archive name must be new.
# Distinct, never-reused archive slot. A retry within the same second, or a
# second failure inside one Slurm job, would otherwise collide on the
# timestamp and hit the "archive already exists" guard.
free_archive_slot() {
    local volume_root="$1"
    local stage_name="$2"
    local attempt=0 archive
    while true; do
        archive="${volume_root}/${stage_name}.failed-${SLURM_JOB_ID:-manual}-$(date -u +%Y%m%dT%H%M%SZ)"
        [[ -e "$archive" ]] || break
        attempt=$((attempt + 1))
        [[ "$attempt" -gt 100 ]] && fail "cannot find a free ${stage_name} archive name under $volume_root"
        sleep 1
    done
    printf '%s\n' "$archive"
}

archive_stage_attempt() {
    local volume_root="$1"
    local stage_root="$2"
    local tag="$3"
    local stage_name="$4"
    [[ "$volume_root" == "$RUN_ROOT/volumes/QHA-${tag}" ]] \
        || fail "refusing unexpected QHA volume path: $volume_root"
    [[ "$stage_root" == "$volume_root/$stage_name" ]] \
        || fail "refusing unexpected ${stage_name} path: $stage_root"
    [[ -d "$stage_root" ]] || return 0

    local archive
    archive="$(free_archive_slot "$volume_root" "$stage_name")"
    mv -- "$stage_root" "$archive"
    printf '[ARCHIVE] unusable %s attempt preserved at %s\n' "$stage_name" "$archive"
    # `mv` discards the source's mtime when the target name does not exist, so
    # the archive would look "older" than the run that produced it. Restore it
    # from the newest file inside so the audit trail stays chronological.
    local newest
    newest="$(find "$archive" -type f -printf '%T@\n' 2>/dev/null | sort -n | tail -1)"
    if [[ -n "$newest" ]]; then
        touch -d "@${newest%.*}" -- "$archive" 2>/dev/null || true
    fi
}

archive_dfpt_attempt() {
    archive_stage_attempt "$1" "$2" "$3" "$DFPT_NAME"
}

# Only the helper-generated POSCAR residue may be removed; any other content is
# a hard failure rather than a recursive deletion.
cleanup_generated_inputs() {
    local dfpt_root="$1"
    local generated="${dfpt_root}/generated-inputs"
    [[ "$generated" == "$RUN_ROOT"/volumes/QHA-*/"$DFPT_NAME"/generated-inputs ]] \
        || fail "refusing to clean unexpected generated-inputs path: $generated"
    [[ -d "$generated" ]] || return 0

    # The helper writes POSCAR in addition to the four files consumed below.
    # Remove only this known, reproducible leftover; any other residue is a
    # hard failure rather than a recursive deletion.
    if [[ -e "$generated/POSCAR" ]]; then
        rm -f -- "$generated/POSCAR"
    fi
    if find "$generated" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
        echo "Unexpected files remain in $generated:" >&2
        find "$generated" -mindepth 1 -maxdepth 1 -printf "  %f\n" >&2
        fail "refusing to remove non-empty generated-inputs directory"
    fi
    rmdir -- "$generated"
}

# ---- read-only resume plan -------------------------------------------------
plan_stage() {
    local dir="$1"
    local stage="$2"
    local label="$3"
    local extra=()
    if [[ "$stage" != "dfpt" ]]; then
        extra=(--max-force-tolerance 0.01 --fail-on-force)
    fi
    if [[ -s "$dir/vasprun.xml" && -s "$dir/OUTCAR" ]]; then
        if verify_quiet "$dir" "$stage" "${extra[@]}"; then
            printf '  [SKIP] %-14s verified VASP outputs\n' "$label"
        else
            printf '  [STOP] %-14s outputs present but FAIL fail-closed convergence check\n' "$label"
        fi
    elif [[ -e "$dir/vasprun.xml" || -e "$dir/OUTCAR" ]]; then
        printf '  [STOP] %-14s partial VASP outputs, needs manual inspection\n' "$label"
    else
        printf '  [RUN ] %-14s no VASP outputs yet\n' "$label"
    fi
}

plan_dfpt() {
    local volume_root="$1"
    local tag="$2"
    local dfpt_root="$volume_root/$DFPT_NAME"
    local seed_root="$volume_root/$SEED_NAME"
    local label="QHA-${tag} ${DFPT_NAME}"
    local seed_label="QHA-${tag} seed${VARIANT_SUFFIX}"
    # The DFPT step needs a converged CHGCAR in its own directory, because its
    # INCAR uses ICHARG=1. Without it VASP starts from the atomic superposition
    # with an all-spin-up guess and the SCF does not converge (job 53681).
    if [[ -s "$seed_root/CHGCAR" ]] \
            && ! "$PYTHON" "$PREP" check-recipe --stage charge-seed \
                    --dfpt-variant "$DFPT_VARIANT" \
                    --manifest "$seed_root/input-manifest.json" > /dev/null 2>&1; then
        printf '  [NOTE] %-22s recipe check flags this CHGCAR as stale (advisory only)\n' "$seed_label"
    elif [[ -s "$seed_root/CHGCAR" ]]; then
        # A CHGCAR on disk is not proof the seed is reusable: the execute path
        # re-verifies with the fail-closed gate before trusting it, so the plan
        # must apply the same gate or it would advertise a stale seed as [SKIP].
        if verify_quiet "$seed_root" charge-seed --max-force-tolerance 0.01 --fail-on-force; then
            printf '  [SKIP] %-22s verified converged CHGCAR present\n' "$seed_label"
        else
            printf '  [STOP] %-22s CHGCAR present but FAILS fail-closed verification\n' "$seed_label"
        fi
    elif [[ -d "$seed_root" ]] && dfpt_has_vasp_artifacts "$seed_root"; then
        printf '  [STOP] %-22s charge-seed outputs present but no usable CHGCAR\n' "$seed_label"
    else
        printf '  [RUN ] %-22s charge seed SCF needed\n' "$seed_label"
    fi
    if [[ -d "$dfpt_root" ]] && dfpt_has_vasp_artifacts "$dfpt_root"; then
        if vasp_complete "$dfpt_root" && verify_quiet "$dfpt_root" dfpt; then
            if [[ -s "$dfpt_root/FORCE_CONSTANTS" && -s "$dfpt_root/thermal_properties-${tag}.yaml" ]]; then
                printf '  [SKIP] %-22s verified DFPT + thermal properties\n' "$label"
            elif [[ -s "$dfpt_root/FORCE_CONSTANTS" ]]; then
                printf '  [PART] %-22s verified DFPT + FORCE_CONSTANTS; thermal still to build\n' "$label"
            else
                printf '  [PART] %-22s verified DFPT; FORCE_CONSTANTS/thermal still to build\n' "$label"
            fi
        else
            printf '  [ARCH] %-22s invalid DFPT outputs -> archive, then re-run DFPT\n' "$label"
        fi
    else
        printf '  [RUN ] %-22s no VASP outputs yet\n' "$label"
    fi
}

plan_resume() {
    echo
    if [[ ! -d "$RUN_ROOT" ]]; then
        printf 'No existing run directory for RUN_ID=%s; nothing to resume.\n' "$RUN_ID"
        echo "Existing runs under ${PILOT_BASE}/${PROFILE}/runs:"
        ls -1 "${PILOT_BASE}/${PROFILE}/runs" 2>/dev/null | sed 's/^/  /' || true
        return 0
    fi
    printf 'RESUME PLAN for %s\n' "$RUN_ROOT"
    printf '  DFPT variant: %s\n' "$DFPT_VARIANT"
    printf '  tags selected by ONLY_TAGS: %s\n' "${ONLY_TAGS[*]}"
    plan_stage "$RUN_ROOT/00-R3" initial-relax "00-R3"
    for i in "${!SCALES[@]}"; do
        tag="${TAGS[$i]}"
        volume_root="$RUN_ROOT/volumes/QHA-${tag}"
        if ! tag_selected "$tag"; then
            printf '  [HOLD] %-14s not in ONLY_TAGS (left untouched)\n' "QHA-${tag}"
            continue
        fi
        plan_stage "$volume_root/R3f-static" volume-static "QHA-${tag} R3f"
        plan_dfpt "$volume_root" "$tag"
    done
    if [[ -s "$RUN_ROOT/Thermal/gibbs-temperature.dat" ]]; then
        printf '  [SKIP] %-14s phonopy-qha already produced\n' "Thermal"
    elif [[ "${#ONLY_TAGS[@]}" -lt "${#TAGS[@]}" ]]; then
        printf '  [HOLD] %-14s partial tag selection: phonopy-qha needs all 9 volumes\n' "Thermal"
    else
        printf '  [RUN ] %-14s phonopy-qha pending\n' "Thermal"
    fi
}

review() {
    [[ "$PROFILE" == "legacy-fixed-grid" || "$PROFILE" == "atomate2-auto-kpoints" ]] \
        || fail "unknown profile: $PROFILE"
    require_file "$PYTHON"
    require_file "$PREP"
    require_file "$PHONOPY"
    require_file "$PHONOPY_QHA"
    "$PYTHON" "$PREP" check-model --root "$HTVASP_ROOT" --potcar-root "$PMG_VASP_PSP_DIR"
    "$PYTHON" "$PREP" audit-kpoints --root "$HTVASP_ROOT" --potcar-root "$PMG_VASP_PSP_DIR"
    cat <<EOF

REVIEW ONLY — no sbatch, srun, VASP, phonopy, or QHA calculation was started.
Selected profile: $PROFILE
Output isolation: $RUN_ROOT
Legacy grids: R3/R3f 10x10x10 (2-site cell); DFPT 5x5x5 (16-site supercell).
Auto profile: generator is applied to the actual calculation cell; DFPT mesh is generated from SPOSCAR.
Volume points: R3f is static (IBRION=-1, NSW=0); no accidental ionic/cell relaxation.
Planned stages: 1 R3 relaxation + 9 static volume points + 9 charge seeds
                + 9 DFPT phonon points + QHA.
Slurm header: partCPU, 32 tasks, 64G, 100 hours, vasp-cpu.
Phonopy DFPT post-processing reads the 2-site POSCAR-unitcell with DIM 2 2 2,
so phonopy rebuilds the same 16-site supercell the DFPT forces came from.
DFPT INCAR starts from B5070's INCAR-Phonopy (ICHARG=1, NELM=120) with the
production override ALGO=Normal + LASPH=True + LMAXMIX=4; the template's
ALGO=Fast diverges on this cell (see the DFPT_VARIANTS comment in the helper).
The charge-seed stage exists to satisfy ICHARG=1: it runs a separate SCF on the
same 16-site SPOSCAR with MAGMOM=16*3.0 and ALGO=Normal (the recipe already
validated on this cell by data/diagnostics/BCC-Ti-V-phonon-20260921-v3) and
writes the CHGCAR that the legacy driver never produced. Jobs 53681 and 53688
failed precisely because that file was missing, leaving VASP on an all-spin-up
atomic guess that does not converge.
EOF
    printf 'DFPT variant: %s (DFPT dir "%s", seed dir "%s", run dir suffix "%s")\n' \
        "$DFPT_VARIANT" "$DFPT_NAME" "$SEED_NAME" "${VARIANT_SUFFIX:-<none>}"
    plan_resume
}

execute() {
    [[ "$PROFILE" == "legacy-fixed-grid" || "$PROFILE" == "atomate2-auto-kpoints" ]] \
        || fail "unknown profile: $PROFILE"
    [[ -n "${SLURM_JOB_ID:-}" ]] || fail "execute mode requires an existing Slurm allocation"
    [[ "${ALLOW_COMPUTE:-NO}" == "YES" ]] || fail "explicit ALLOW_COMPUTE=YES gate is required"
    [[ -x "$PYTHON" && -x "$PHONOPY" && -x "$PHONOPY_QHA" ]] || fail "ht-vasp Python/phonopy tools are missing"
    [[ -d "$PMG_VASP_PSP_DIR/potpaw_PBE" ]] || fail "PBE POTCAR library missing below $PMG_VASP_PSP_DIR"

    if [[ "$RESUME" == "YES" ]]; then
        [[ -d "$RUN_ROOT" ]] || fail "resume requested but run directory does not exist: $RUN_ROOT"
        require_file "$RUN_ROOT/run-manifest.json"
    else
        [[ ! -e "$RUN_ROOT" ]] || fail "run directory already exists; use RESUME=YES to continue: $RUN_ROOT"
    fi

    module purge
    module load vasp-cpu
    ulimit -s unlimited
    export PMG_VASP_PSP_DIR I_MPI_FABRICS
    cd "$HTVASP_ROOT"
    mkdir -p "$RUN_ROOT"
    if [[ "$RESUME" == "YES" ]]; then
        printf '{"resumed_at_utc":"%s","host":"%s","I_MPI_FABRICS":"%s","only_tags":"%s","dfpt_variant":"%s"}\n' \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(hostname)" "$I_MPI_FABRICS" "${ONLY_TAGS[*]}" "$DFPT_VARIANT" \
            > "$RUN_ROOT/resume-${DFPT_VARIANT}-$(date -u +%Y%m%dT%H%M%SZ).json"
    else
        printf '{"profile":"%s","run_id":"%s","created_at_utc":"%s","host":"%s","I_MPI_FABRICS":"%s","dfpt_variant":"%s"}\n' \
            "$PROFILE" "$RUN_ID" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(hostname)" "$I_MPI_FABRICS" "$DFPT_VARIANT" \
            > "$RUN_ROOT/run-manifest.json"
    fi

    # ---- 1. R3 initial relaxation ----------------------------------------
    r3_root="$RUN_ROOT/00-R3"
    if [[ -s "$r3_root/vasprun.xml" && -s "$r3_root/OUTCAR" ]]; then
        if ! verify_vasp "$r3_root" initial-relax --max-force-tolerance 0.01 --fail-on-force; then
            fail "existing R3 fails fail-closed verification; evidence left in place at $r3_root (see vasp-verification.json)"
        fi
        echo "[SKIP] verified R3: $r3_root"
    else
        [[ ! -e "$r3_root/vasprun.xml" && ! -e "$r3_root/OUTCAR" ]] \
            || fail "R3 has partial VASP outputs; inspect manually: $r3_root"
        "$PYTHON" "$PREP" prepare --root "$HTVASP_ROOT" --potcar-root "$PMG_VASP_PSP_DIR" \
            --profile "$PROFILE" --stage initial-relax --outdir "$r3_root"
        run_vasp "$r3_root"
        if ! verify_vasp "$r3_root" initial-relax --max-force-tolerance 0.01 --fail-on-force; then
            fail "R3 failed fail-closed verification; inspect $r3_root/vasp-verification.json"
        fi
    fi
    require_file "$r3_root/CONTCAR"

    for i in "${!SCALES[@]}"; do
        tag="${TAGS[$i]}"
        scale="${SCALES[$i]}"
        if ! tag_selected "$tag"; then
            echo "[HOLD] QHA-${tag} not selected by ONLY_TAGS=${ONLY_TAGS[*]}; leaving it untouched"
            continue
        fi
        volume_root="$RUN_ROOT/volumes/QHA-${tag}"
        static_root="$volume_root/R3f-static"
        dfpt_root="$volume_root/$DFPT_NAME"
        mkdir -p "$volume_root"

        # ---- 2. R3f static volume point ----------------------------------
        # Fixed scaled cell, static calculation matching B5070 INCAR-R3f.
        if [[ -s "$static_root/vasprun.xml" && -s "$static_root/OUTCAR" ]]; then
            if ! verify_vasp "$static_root" volume-static --max-force-tolerance 0.01 --fail-on-force; then
                fail "existing QHA-${tag} R3f fails fail-closed verification; evidence left in place at $static_root"
            fi
            echo "[SKIP] verified R3f-static QHA-${tag}: $static_root"
        else
            [[ ! -e "$static_root/vasprun.xml" && ! -e "$static_root/OUTCAR" ]] \
                || fail "QHA-${tag} R3f has partial VASP outputs; inspect manually: $static_root"
            "$PYTHON" "$PREP" prepare --root "$HTVASP_ROOT" --potcar-root "$PMG_VASP_PSP_DIR" \
                --profile "$PROFILE" --stage volume-static --scale "$scale" \
                --structure "$r3_root/CONTCAR" --outdir "$static_root"
            run_vasp "$static_root"
            if ! verify_vasp "$static_root" volume-static --max-force-tolerance 0.01 --fail-on-force; then
                fail "QHA-${tag} R3f failed fail-closed verification; inspect $static_root/vasp-verification.json"
            fi
        fi
        require_file "$static_root/CONTCAR"
        if [[ ! -s "$volume_root/eos-record.json" ]]; then
            "$PYTHON" "$PREP" record-e --scale "$scale" \
                --vasprun "$static_root/vasprun.xml" --output "$volume_root/eos-record.json"
        else
            echo "[SKIP] existing EOS record: $volume_root/eos-record.json"
        fi

        # ---- 3. DFPT preflight: keep verified results, archive the rest ---
        # An existing DFPT directory is reusable only when its VASP outputs
        # pass the fail-closed convergence check. Anything else — a run that
        # exited cleanly but never reached EDIFF, a killed job, partial output
        # — is archived whole so the failure stays auditable.
        dfpt_reusable=NO
        if [[ -d "$dfpt_root" ]] && dfpt_has_vasp_artifacts "$dfpt_root"; then
            # Verification is attempted only on a finished run: a killed job
            # (e.g. 53681) has no vasprun.xml at all and must go straight to
            # the archive instead of aborting the script on a missing file.
            if vasp_complete "$dfpt_root" && verify_vasp "$dfpt_root" dfpt; then
                dfpt_reusable=YES
                echo "[SKIP] verified DFPT VASP outputs QHA-${tag}: $dfpt_root"
            else
                archive_dfpt_attempt "$volume_root" "$dfpt_root" "$tag"
            fi
        fi
        mkdir -p "$dfpt_root"

        # Generate phonopy supercell once, then generate the actual DFPT VASP
        # inputs from SPOSCAR (not from the primitive unit cell).
        if [[ ! -s "$dfpt_root/SPOSCAR" ]]; then
            cp -- "$static_root/CONTCAR" "$dfpt_root/POSCAR-unitcell"
            cat > "$dfpt_root/mesh.conf" <<'PHONOPY_CONF'
ATOM_NAME = Ti V
DIM = 2 2 2
MP = 24 24 24
FORCE_CONSTANTS = READ
PHONOPY_CONF
            (
                cd "$dfpt_root"
                "$PHONOPY" -d --dim 2 2 2 -c POSCAR-unitcell > phonopy-generate.log 2>&1
                require_file SPOSCAR
            )
        else
            echo "[SKIP] existing phonopy supercell: $dfpt_root/SPOSCAR"
        fi
        require_file "$dfpt_root/SPOSCAR"
        require_file "$dfpt_root/mesh.conf"
        require_file "$dfpt_root/POSCAR-unitcell"

        # Prepare DFPT inputs only when the final input set is incomplete.
        if [[ ! -s "$dfpt_root/INCAR" || ! -s "$dfpt_root/KPOINTS" || \
              ! -s "$dfpt_root/POTCAR" || ! -s "$dfpt_root/input-manifest.json" ]]; then
            cleanup_generated_inputs "$dfpt_root"
            "$PYTHON" "$PREP" prepare --root "$HTVASP_ROOT" \
                --potcar-root "$PMG_VASP_PSP_DIR" --profile "$PROFILE" \
                --stage dfpt --dfpt-variant "$DFPT_VARIANT" \
                --structure "$dfpt_root/SPOSCAR" \
                --outdir "$dfpt_root/generated-inputs"
            for file in INCAR KPOINTS POTCAR input-manifest.json; do
                require_file "$dfpt_root/generated-inputs/$file"
                mv -- "$dfpt_root/generated-inputs/$file" "$dfpt_root/$file"
            done
        fi
        # Remove only the known helper-generated POSCAR residue; fail closed
        # if any other unexpected file remains in generated-inputs.
        cleanup_generated_inputs "$dfpt_root"
        require_file "$dfpt_root/INCAR"
        require_file "$dfpt_root/KPOINTS"
        require_file "$dfpt_root/POTCAR"
        cp -- "$dfpt_root/SPOSCAR" "$dfpt_root/POSCAR"

        # ---- 4. charge-seed SCF ------------------------------------------
        # DFPT runs ICHARG=1, i.e. "read the initial charge density". The legacy
        # driver never wrote one (R3 and R3f both use LCHARG=False), so VASP fell
        # back to the atomic superposition with an all-spin-up guess and job
        # 53681 failed to converge. Produce the charge density the DFPT template
        # always assumed was there, using MAGMOM + ALGO=Normal, keeping the DFPT
        # INCAR itself byte-for-byte equal to B5070's INCAR-Phonopy.
        seed_root="$volume_root/$SEED_NAME"
        if [[ -s "$seed_root/CHGCAR" ]]; then
            # A CHGCAR is only reusable if it was built by the CURRENT recipe.
            # After the 53690 diagnosis the seed recipe changed (its PAW tags
            # were dropped to match DFPT's Hamiltonian), so a seed from the old
            # recipe is stale evidence, not a shortcut.
            #
            # This check REPORTS, it does not delete. Job 53701 archived a seed
            # that a manual re-run of this very command scores as "current", so
            # a non-zero exit here is not by itself grounds for destroying a
            # 32 MB converged charge density. The fail-closed verification below
            # is what decides reuse, and check-stage-compat still gates the
            # seed/DFPT Hamiltonian match before any DFPT force is produced.
            recipe_rc=0
            "$PYTHON" "$PREP" check-recipe --stage charge-seed \
                --dfpt-variant "$DFPT_VARIANT" \
                --manifest "$seed_root/input-manifest.json" \
                > "$seed_root/check-recipe.json" 2>&1 || recipe_rc=$?
            if [[ "$recipe_rc" -ne 0 ]]; then
                echo "[WARN] QHA-${tag} charge seed does not match the current recipe (rc=$recipe_rc):"
                grep -E '"(status|changed_tags)"' "$seed_root/check-recipe.json" 2>/dev/null \
                    | sed 's/^/       /' || true
            fi
            # The seed's job is a converged charge density, not a force
            # criterion, so this gates on SCF convergence + CHGCAR only.
            if verify_vasp "$seed_root" charge-seed; then
                echo "[SKIP] verified charge seed QHA-${tag}: $seed_root"
            else
                fail "QHA-${tag} charge seed has a CHGCAR but fails verification; inspect $seed_root"
            fi
        else
            if dfpt_has_vasp_artifacts "$seed_root"; then
                archive_stage_attempt "$volume_root" "$seed_root" "$tag" "$SEED_NAME"
                # The archive must have actually taken effect *for this process*.
                # Jobs 53701 and 53702 both reached here, reported a successful
                # archive, and then ran VASP straight into a directory that in
                # fact still held the 53697 outputs — overwriting a verified
                # 32 MB CHGCAR with a fresh, unrelated SCF. A second observation
                # is cheap; destroying a converged charge density is not.
                if dfpt_has_vasp_artifacts "$seed_root"; then
                    fail "charge seed still holds VASP artifacts after archiving (QHA-${tag}); refusing to run VASP over existing results: $seed_root"
                fi
            fi
            "$PYTHON" "$PREP" prepare --root "$HTVASP_ROOT" --potcar-root "$PMG_VASP_PSP_DIR" \
                --profile "$PROFILE" --stage charge-seed --dfpt-variant "$DFPT_VARIANT" \
                --structure "$dfpt_root/SPOSCAR" --outdir "$seed_root"
            # The seed must run on precisely the cell DFPT will run on; the
            # helper rewrites POSCAR to its SPOSCAR source for that reason.
            if ! cmp -s "$seed_root/POSCAR" "$dfpt_root/SPOSCAR"; then
                fail "charge seed POSCAR differs from DFPT SPOSCAR; refusing to supply a mismatched charge density"
            fi
            run_vasp "$seed_root"
            if ! verify_vasp "$seed_root" charge-seed; then
                archive_stage_attempt "$volume_root" "$seed_root" "$tag" "$SEED_NAME"
                fail "charge seed QHA-${tag} failed fail-closed verification and was archived; review it before resubmitting"
            fi
        fi
        require_file "$seed_root/CHGCAR"

        # The supplied density is only the converged density of the Hamiltonian
        # DFPT is about to solve if both stages agree on the tags that define
        # it. Job 53690 failed exactly here.
        if ! "$PYTHON" "$PREP" check-stage-compat \
                --seed-root "$seed_root" --dfpt-root "$dfpt_root" > "$dfpt_root/seed-compat.json"; then
            cat "$dfpt_root/seed-compat.json" >&2
            fail "QHA-${tag} seed/DFPT Hamiltonian mismatch; see $dfpt_root/seed-compat.json"
        fi

        # Hand the converged charge density to DFPT, then make DFPT own it: the
        # shell driver re-copies it on every pass so a resumed run cannot start
        # from a stale or absent file.
        cp -- "$seed_root/CHGCAR" "$dfpt_root/CHGCAR"
        require_file "$dfpt_root/CHGCAR"

        # ---- 5. DFPT VASP ------------------------------------------------
        # CHGCAR alone is not evidence that VASP ran here: the driver wrote it
        # on purpose just above. Use the real output set for that judgement.
        if [[ "$dfpt_reusable" == "YES" ]]; then
            echo "[SKIP] reusing verified DFPT VASP outputs QHA-${tag}"
        else
            # The preflight above removed every artifact-bearing directory, so
            # a fresh run may only start from a clean state (CHGCAR excepted).
            vasp_ran_here "$dfpt_root" \
                && fail "DFPT has VASP output after the archive preflight: $dfpt_root"
            run_vasp "$dfpt_root"
            if ! verify_vasp "$dfpt_root" dfpt; then
                # Archive the attempt (with its verification report and its
                # supplied CHGCAR) instead of letting phonopy consume an
                # unconverged vasprun.xml.
                archive_dfpt_attempt "$volume_root" "$dfpt_root" "$tag"
                fail "DFPT QHA-${tag} failed fail-closed verification and was archived; review it before resubmitting"
            fi
        fi

        # ---- 6. FORCE_CONSTANTS + thermal properties ---------------------
        # Reached only with a verified DFPT vasprun.xml in place.
        (
            cd "$dfpt_root"
            # Build the force constants FIRST. This audit reads FORCE_CONSTANTS
            # as one of its inputs, so running it before the file exists fails
            # the run for the wrong reason.
            if [[ ! -s FORCE_CONSTANTS ]]; then
                "$PHONOPY" --fc vasprun.xml > phonopy-force-constants.log 2>&1
            else
                echo "[SKIP] existing FORCE_CONSTANTS"
            fi
            require_file FORCE_CONSTANTS
            # Now verify the whole data flow: unit cell -> DIM -> supercell vs
            # SPOSCAR vs the 16x16 force constants vs the DFPT vasprun.
            if ! "$PYTHON" "$PREP" check-phonopy-inputs --dfpt-root "$dfpt_root" \
                    --vasprun "$dfpt_root/vasprun.xml" > phonopy-input-audit.json; then
                cat phonopy-input-audit.json >&2
                fail "phonopy DFPT inputs are inconsistent for QHA-${tag}"
            fi
            thermal_file="thermal_properties-${tag}.yaml"
            if [[ ! -s "$thermal_file" ]]; then
                # -c POSCAR-unitcell is required: POSCAR in this directory is
                # phonopy's 16-site SPOSCAR, and applying DIM=2 2 2 to it again
                # would build a 128-site supercell that cannot be mapped onto
                # the 16x16 FORCE_CONSTANTS.
                "$PHONOPY" -t --mesh 24 24 24 --tmin 0 --tmax 3000 --tstep 50 \
                    mesh.conf -c POSCAR-unitcell > phonopy-thermal.log 2>&1
                require_file thermal_properties.yaml
                mv -- thermal_properties.yaml "$thermal_file"
            else
                echo "[SKIP] existing thermal properties: $thermal_file"
            fi
        )
    done

    # ---- 7. phonopy-qha ---------------------------------------------------
    # phonopy-qha needs all nine (v, E) points and nine thermal-property files.
    # A partial ONLY_TAGS run therefore stops after the selected volume points.
    thermal_root="$RUN_ROOT/Thermal"
    if [[ "$DFPT_VARIANT" != "production" ]]; then
        echo "[HOLD] DFPT_VARIANT=$DFPT_VARIANT is an experiment, not the production recipe;"
        echo "       phonopy-qha is produced only from the production variant."
        echo "Pilot stage run finished for: $RUN_ROOT"
        return 0
    fi
    if [[ "${#ONLY_TAGS[@]}" -lt "${#TAGS[@]}" ]]; then
        echo "[HOLD] ONLY_TAGS=${ONLY_TAGS[*]} is a partial set; phonopy-qha is deferred"
        echo "       Re-submit with ONLY_TAGS unset (all nine) once the smoke test passes."
        echo "Pilot stage run finished for: $RUN_ROOT"
        return 0
    fi
    if [[ -s "$thermal_root/gibbs-temperature.dat" ]]; then
        echo "[SKIP] completed phonopy-qha: $thermal_root"
    else
        "$PYTHON" "$PREP" write-ve --records "$RUN_ROOT"/volumes/QHA-*/eos-record.json \
            --output "$RUN_ROOT/v-e.dat"
        mkdir -p "$thermal_root"
        cp -- "$RUN_ROOT/v-e.dat" "$thermal_root/v-e.dat"
        for tag in "${TAGS[@]}"; do
            require_file "$RUN_ROOT/volumes/QHA-${tag}/$DFPT_NAME/thermal_properties-${tag}.yaml"
            cp -- "$RUN_ROOT/volumes/QHA-${tag}/$DFPT_NAME/thermal_properties-${tag}.yaml" \
                "$thermal_root/thermal_properties-${tag}.yaml"
        done
        (
            cd "$thermal_root"
            "$PHONOPY_QHA" v-e.dat thermal_properties-*.yaml --tmax 3000 > phonopy-qha.log 2>&1
            require_file gibbs-temperature.dat
        )
    fi
    echo "Pilot completed/resumed: $RUN_ROOT"
}
case "${1:-review}" in
    review) review ;;
    execute) execute ;;
    *) echo "Usage: $0 [review|execute] [legacy-fixed-grid|atomate2-auto-kpoints]" >&2; exit 2 ;;
esac
