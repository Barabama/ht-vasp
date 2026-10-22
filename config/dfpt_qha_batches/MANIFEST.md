# DFPT-QHA endpoint manifest

Derived from each endpoint's `*-dfpt-qha.json` on the cluster that produced it.
Regenerate: `scripts/gather_status.py <cluster> <out>` on each cluster, merge,
then the table builder at the bottom of this file.

| endpoint    | class   | cell | supercell | magnetic         | batch  | status                                                          | qha_applicability       |
|-------------|---------|------|-----------|------------------|--------|-----------------------------------------------------------------|-------------------------|
| BCC-Al-Al-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | queued B5070 (job 2325, NCORE=16); 429pro attempt failed        |                         |
| BCC-Al-Cr-2 | BCC_B2  | 2    | 16        | confirmed_branch | 429pro | done 429pro                                                     | harmonic_QHA_applicable |
| BCC-Al-Hf-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| BCC-Al-Nb-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| BCC-Al-Ta-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| BCC-Al-Ti-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| BCC-Al-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| BCC-Cr-Cr-2 | BCC_B2  | 2    | 16        | confirmed_branch | -      | not submitted                                                   |                         |
| BCC-Cr-Hf-2 | BCC_B2  | 2    | 16        | confirmed_branch | -      | not submitted                                                   |                         |
| BCC-Cr-Nb-2 | BCC_B2  | 2    | 16        | nonmagnetic      | -      | not submitted                                                   |                         |
| BCC-Cr-Ta-2 | BCC_B2  | 2    | 16        | nonmagnetic      | -      | not submitted                                                   |                         |
| BCC-Cr-Ti-2 | BCC_B2  | 2    | 16        | confirmed_branch | -      | not submitted                                                   |                         |
| BCC-Cr-Zr-2 | BCC_B2  | 2    | 16        | confirmed_branch | -      | not submitted                                                   |                         |
| BCC-Hf-Hf-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| BCC-Hf-Nb-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| BCC-Hf-Ta-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| BCC-Hf-Ti-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| BCC-Hf-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| BCC-Nb-Nb-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | harmonic_QHA_applicable |
| BCC-Nb-Ta-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | harmonic_QHA_applicable |
| BCC-Nb-Ti-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | numerical_unresolved    |
| BCC-Nb-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | finite_q_soft_mode      |
| BCC-Ta-Ta-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                                     | unknown                 |
| BCC-Ta-Ti-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | finite_q_soft_mode      |
| BCC-Ta-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | finite_q_soft_mode      |
| BCC-Ti-Ti-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | finite_q_soft_mode      |
| BCC-Ti-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | finite_q_soft_mode      |
| BCC-Zr-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | finite_q_soft_mode      |
| FCC-Al-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                                      | harmonic_QHA_applicable |
| FCC-Al-Cr-4 | FCC_L12 | 4    | 32        | confirmed_branch | 429pro | done 429pro                                                     | gamma_soft_mode         |
| FCC-Al-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                                     | harmonic_QHA_applicable |
| FCC-Al-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Al-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Al-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| FCC-Al-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Cr-Al-4 | FCC_L12 | 4    | 32        | confirmed_branch | 429pro | done 429pro                                                     | harmonic_QHA_applicable |
| FCC-Cr-Cr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | -      | not submitted                                                   |                         |
| FCC-Cr-Hf-4 | FCC_L12 | 4    | 32        | confirmed_branch | -      | not submitted                                                   |                         |
| FCC-Cr-Nb-4 | FCC_L12 | 4    | 32        | confirmed_branch | -      | not submitted                                                   |                         |
| FCC-Cr-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | -      | not submitted                                                   |                         |
| FCC-Cr-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | -      | not submitted                                                   |                         |
| FCC-Cr-Zr-4 | FCC_L12 | 4    | 32        | confirmed_branch | -      | not submitted                                                   |                         |
| FCC-Hf-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                                     | harmonic_QHA_applicable |
| FCC-Hf-Cr-4 | FCC_L12 | 4    | 32        | confirmed_branch | -      | not submitted                                                   |                         |
| FCC-Hf-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| FCC-Hf-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Hf-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                                     | finite_q_soft_mode      |
| FCC-Hf-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Hf-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Nb-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Nb-Cr-4 | FCC_L12 | 4    | 32        | confirmed_branch | -      | not submitted                                                   |                         |
| FCC-Nb-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Nb-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Nb-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Nb-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                                   |                         |
| FCC-Nb-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | partial 3/7 volumes, timed out 10-21 18:03, resubmitted as 2326 |                         |
| FCC-Ta-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                                      | harmonic_QHA_applicable |
| FCC-Ta-Cr-4 | FCC_L12 | 4    | 32        | confirmed_branch | -      | not submitted                                                   |                         |
| FCC-Ta-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | partial 4/7 volumes, timed out 10-21 23:07, needs --resume-from |                         |
| FCC-Ta-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | partial 4/7 volumes, timed out 10-22 18:04, needs --resume-from |                         |
| FCC-Ta-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                                      | finite_q_soft_mode      |
| FCC-Ta-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Ta-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Ti-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Ti-Cr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | -      | not submitted                                                   |                         |
| FCC-Ti-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Ti-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Ti-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Ti-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Ti-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Zr-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Zr-Cr-4 | FCC_L12 | 4    | 32        | confirmed_branch | -      | not submitted                                                   |                         |
| FCC-Zr-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Zr-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Zr-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Zr-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| FCC-Zr-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| SER-Al-1    | SER_FCC | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                                    |                         |
| SER-Cr-1    | SER_BCC | 2    | 16        | confirmed_branch | -      | not submitted                                                   |                         |
| SER-Hf-1    | SER_HCP | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | harmonic_QHA_applicable |
| SER-Nb-1    | SER_BCC | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | harmonic_QHA_applicable |
| SER-Ta-1    | SER_BCC | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | harmonic_QHA_applicable |
| SER-Ti-1    | SER_HCP | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | harmonic_QHA_applicable |
| SER-Zr-1    | SER_HCP | 2    | 16        | nonmagnetic      | B5070  | done B5070                                                      | harmonic_QHA_applicable |

**35 / 84 done.** 3 more have partial scratch state that
`--resume-from` can continue. The remaining are queued or not submitted.


## How to regenerate

```bash
# on each cluster, from the repo root
./.conda/bin/python scripts/gather_status.py 429pro /tmp/status_429.txt
./.conda/bin/python scripts/gather_status.py B5070  /tmp/status_507.txt
```

Merge the two lists (later wins for an endpoint that ran on both), then rebuild
the table from `pool_recipes.json` for the class / cell / magnetic columns.  The
`status` column is the ground truth: it comes from the artifact, not from a
booking sheet.

## Resuming a run that hit the time limit

A job killed at its Slurm limit gets SIGTERM, so `Worker.run_flow` never reaches
its `finally` block and `store_dir` is never created.  The completed volumes sit
in `/nfs_ssd/tmp/<name>-<uuid>` instead.  Continue from there without copying:

```bash
./.conda/bin/python main_dfpt_qha.py -n <endpoint> --slurm --resume \
    --resume-from /nfs_ssd/tmp/<name>-<uuid>
```

`--resume-from` is threaded into the submitted command, so the job itself sees
it.  Without `--resume` it is rejected rather than silently ignored.
