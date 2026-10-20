# DFPT-QHA endpoint manifest

Status is derived from each endpoint's `*-dfpt-qha.json` on the cluster that
ran it; the ledger and that file remain the ground truth.  Regenerate with the
snippet at the bottom of this file.

| endpoint | class | cell | supercell | magnetic | batch | status |
|---|---|---|---|---|---|---|
| BCC-Al-Al-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | queued A (429pro) |
| BCC-Al-Cr-2 | BCC_B2 | 2 | 16 | confirmed_branch | A (429pro) | done 429pro: ok=true, harmonic_QHA_applicable |
| BCC-Al-Hf-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | done 429pro: ok=true, finite_q_soft_mode |
| BCC-Al-Nb-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | done 429pro: ok=true, finite_q_soft_mode |
| BCC-Al-Ta-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | done 429pro: ok=true, finite_q_soft_mode |
| BCC-Al-Ti-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | done 429pro: ok=true, finite_q_soft_mode |
| BCC-Al-Zr-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | done 429pro: ok=true, finite_q_soft_mode |
| BCC-Cr-Cr-2 | BCC_B2 | 2 | 16 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| BCC-Cr-Hf-2 | BCC_B2 | 2 | 16 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| BCC-Cr-Nb-2 | BCC_B2 | 2 | 16 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| BCC-Cr-Ta-2 | BCC_B2 | 2 | 16 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| BCC-Cr-Ti-2 | BCC_B2 | 2 | 16 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| BCC-Cr-Zr-2 | BCC_B2 | 2 | 16 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| BCC-Hf-Hf-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | done 429pro: ok=true, finite_q_soft_mode |
| BCC-Hf-Nb-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | queued A (429pro) |
| BCC-Hf-Ta-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | done 429pro: ok=true, finite_q_soft_mode |
| BCC-Hf-Ti-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | queued A (429pro) |
| BCC-Hf-Zr-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | done 429pro: ok=true, finite_q_soft_mode |
| BCC-Nb-Nb-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | done 429pro: ok=false (old schema) |
| BCC-Nb-Ta-2 | BCC_B2 | 2 | 16 | nonmagnetic | B (B5070) | done B5070: ok=true, harmonic_QHA_applicable |
| BCC-Nb-Ti-2 | BCC_B2 | 2 | 16 | nonmagnetic | B (B5070) | done B5070: ok=true, numerical_unresolved |
| BCC-Nb-Zr-2 | BCC_B2 | 2 | 16 | nonmagnetic | B (B5070) | done B5070: ok=true, finite_q_soft_mode |
| BCC-Ta-Ta-2 | BCC_B2 | 2 | 16 | nonmagnetic | A (429pro) | done 429pro: ok=true (old schema) |
| BCC-Ta-Ti-2 | BCC_B2 | 2 | 16 | nonmagnetic | B (B5070) | done B5070: ok=true, finite_q_soft_mode |
| BCC-Ta-Zr-2 | BCC_B2 | 2 | 16 | nonmagnetic | B (B5070) | done B5070: ok=true, finite_q_soft_mode |
| BCC-Ti-Ti-2 | BCC_B2 | 2 | 16 | nonmagnetic | B (B5070) | done B5070: ok=true, finite_q_soft_mode |
| BCC-Ti-Zr-2 | BCC_B2 | 2 | 16 | nonmagnetic | B (B5070) | done B5070: ok=true, finite_q_soft_mode |
| BCC-Zr-Zr-2 | BCC_B2 | 2 | 16 | nonmagnetic | B (B5070) | done B5070: ok=true, finite_q_soft_mode |
| FCC-Al-Al-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | done 429pro: ok=true (old schema) · done B5070: ok=true, harmonic_QHA_applicable |
| FCC-Al-Cr-4 | FCC_L12 | 4 | 32 | confirmed_branch | A (429pro) | queued A (429pro) |
| FCC-Al-Hf-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Al-Nb-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Al-Ta-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Al-Ti-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Al-Zr-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Cr-Al-4 | FCC_L12 | 4 | 32 | confirmed_branch | A (429pro) | queued A (429pro) |
| FCC-Cr-Cr-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Cr-Hf-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Cr-Nb-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Cr-Ta-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Cr-Ti-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Cr-Zr-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Hf-Al-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Hf-Cr-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Hf-Hf-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Hf-Nb-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Hf-Ta-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Hf-Ti-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Hf-Zr-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Nb-Al-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Nb-Cr-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Nb-Hf-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Nb-Nb-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Nb-Ta-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Nb-Ti-4 | FCC_L12 | 4 | 32 | nonmagnetic | A (429pro) | queued A (429pro) |
| FCC-Nb-Zr-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ta-Al-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ta-Cr-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Ta-Hf-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ta-Nb-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ta-Ta-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ta-Ti-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ta-Zr-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ti-Al-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ti-Cr-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Ti-Hf-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ti-Nb-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ti-Ta-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ti-Ti-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Ti-Zr-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Zr-Al-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Zr-Cr-4 | FCC_L12 | 4 | 32 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| FCC-Zr-Hf-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Zr-Nb-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Zr-Ta-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Zr-Ti-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| FCC-Zr-Zr-4 | FCC_L12 | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| SER-Al-1 | SER_FCC | 4 | 32 | nonmagnetic | B (B5070) | queued B (B5070) |
| SER-Cr-1 | SER_BCC | 2 | 16 | unresolved | - | blocked: Cr magnetic_policy=unresolved |
| SER-Hf-1 | SER_HCP | 2 | 16 | nonmagnetic | B (B5070) | queued B (B5070) |
| SER-Nb-1 | SER_BCC | 2 | 16 | nonmagnetic | B (B5070) | queued B (B5070) |
| SER-Ta-1 | SER_BCC | 2 | 16 | nonmagnetic | B (B5070) | queued B (B5070) |
| SER-Ti-1 | SER_HCP | 2 | 16 | nonmagnetic | B (B5070) | queued B (B5070) |
| SER-Zr-1 | SER_HCP | 2 | 16 | nonmagnetic | B (B5070) | queued B (B5070) |

Totals: 84 endpoints = 33 batch A + 33 batch B + 18 blocked (Cr).

## Cr spin screening (2026-10-20)

Free relaxations, ISIF=3, nm/fm(/afm). `gain` is the best admissible state
below the non-magnetic one; below 5 meV/atom the endpoint is called
nonmagnetic.  AFM states that changed the cell by more than 2% are excluded
(the L1_2 cell collapses to tetragonal under a staggered start).

| endpoint | decision | gain (meV/atom) | note |
|---|---|---|---|
| BCC-Cr-Cr-2 | confirmed_branch | 14.5 |  |
| BCC-Cr-Hf-2 | confirmed_branch | 77.5 |  |
| BCC-Cr-Nb-2 | nonmagnetic | -0.0 |  |
| BCC-Cr-Ta-2 | nonmagnetic | -7.0 |  |
| BCC-Cr-Ti-2 | confirmed_branch | 21.4 |  |
| BCC-Cr-Zr-2 | confirmed_branch | 91.5 |  |
| FCC-Cr-Cr-4 | nonmagnetic | 0.0 |  |
| FCC-Cr-Hf-4 | confirmed_branch | 102.9 |  |
| FCC-Cr-Nb-4 | confirmed_branch | 63.8 |  |
| FCC-Cr-Ta-4 | nonmagnetic | 0.0 |  |
| FCC-Cr-Ti-4 | nonmagnetic | -0.0 |  |
| FCC-Cr-Zr-4 | confirmed_branch | 128.4 |  |
| FCC-Hf-Cr-4 | confirmed_branch | 15.7 | afm excluded (cell collapse) |
| FCC-Nb-Cr-4 | confirmed_branch | 6.5 | afm excluded (cell collapse); constrained check: fm |
| FCC-Ta-Cr-4 | confirmed_branch | 5.8 | afm excluded (cell collapse); constrained check: fm |
| FCC-Ti-Cr-4 | nonmagnetic | 2.9 | afm excluded (cell collapse) |
| FCC-Zr-Cr-4 | confirmed_branch | 41.1 | afm excluded (cell collapse) |
| SER-Cr-1 | confirmed_branch | 14.5 |  |

## Constrained AFM check (ISIF=2 at the fm geometry)

| endpoint | nm vs fm | afm vs fm | decision |
|---|---|---|---|
| FCC-Ta-Cr-4 | +0.2 meV/atom | +1.0 meV/atom | fm (the nm start itself converged with |m|=1.73) |
| FCC-Nb-Cr-4 | +0.3 meV/atom | did not converge (3 custodian corrections) | fm (the nm start converged with 3×|m|=0.62) |

## Notes

- BCC-Nb-Ti-2 came back `numerical_unresolved` on volumes 3, 4 and 6: the raw
  DFPT Hessian has an acoustic-sum-rule residual of 13–16%, above the 10%
  above which the matrix is not repaired.  The endpoint is not usable for QHA.
- `BCC-Al-Ta-2` has two run directories: the old failed run and the current
  one (`finite_q_soft_mode`, V0 = 33.59 A^3).  Runs are immutable; the newest
  is the record.

## Submission log

| date | batch | cluster | Slurm jobs | note |
|---|---|---|---|---|
| 2026-10-20 | A | 429pro | 56962–56994 (33) | 29 new + 4 regenerated under the current schema |
| 2026-10-20 | B | B5070 | 2291–2323 (33) | 33 new |

## Refresh

```bash
python - <<'EOF'
import glob, json
from pathlib import Path
pool = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")
for name in sorted(json.loads(pool.read_text())["endpoints"]):
    hits = glob.glob(f"data/endpoints/{name}/dfpt-qhaflow/runs/*/{name}-dfpt-qha.json")
    print(f"{name:16s}", "done" if hits else "-")
EOF
```
