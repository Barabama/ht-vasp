# DFPT-QHA endpoint manifest

Derived from each endpoint's `*-dfpt-qha.json` on the cluster that produced
it -- not from a booking sheet.  Regenerate with `scripts/gather_status.py`
on each cluster, then `scripts/build_manifest.py`.

| endpoint    | class   | cell | supercell | magnetic         | batch  | status                                          | qha_applicability       |
|-------------|---------|------|-----------|------------------|--------|-------------------------------------------------|-------------------------|
| BCC-Al-Al-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| BCC-Al-Cr-2 | BCC_B2  | 2    | 16        | confirmed_branch | 429pro | done 429pro                                     | harmonic_QHA_applicable |
| BCC-Al-Hf-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| BCC-Al-Nb-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| BCC-Al-Ta-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| BCC-Al-Ti-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| BCC-Al-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| BCC-Cr-Cr-2 | BCC_B2  | 2    | 16        | confirmed_branch | B5070  | done B5070                                      | harmonic_QHA_applicable |
| BCC-Cr-Hf-2 | BCC_B2  | 2    | 16        | confirmed_branch | B5070  | done B5070                                      | harmonic_QHA_applicable |
| BCC-Cr-Nb-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| BCC-Cr-Ta-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| BCC-Cr-Ti-2 | BCC_B2  | 2    | 16        | confirmed_branch | B5070  | done B5070                                      | harmonic_QHA_applicable |
| BCC-Cr-Zr-2 | BCC_B2  | 2    | 16        | confirmed_branch | B5070  | done B5070                                      | harmonic_QHA_applicable |
| BCC-Hf-Hf-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| BCC-Hf-Nb-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | 429pro partial, 6/7 volumes done, missing 5     |                         |
| BCC-Hf-Ta-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| BCC-Hf-Ti-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| BCC-Hf-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| BCC-Nb-Nb-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | harmonic_QHA_applicable |
| BCC-Nb-Ta-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| BCC-Nb-Ti-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | numerical_unresolved    |
| BCC-Nb-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| BCC-Ta-Ta-2 | BCC_B2  | 2    | 16        | nonmagnetic      | 429pro | done 429pro                                     | harmonic_QHA_applicable |
| BCC-Ta-Ti-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| BCC-Ta-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| BCC-Ti-Ti-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| BCC-Ti-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| BCC-Zr-Zr-2 | BCC_B2  | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| FCC-Al-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| FCC-Al-Cr-4 | FCC_L12 | 4    | 32        | confirmed_branch | 429pro | done 429pro                                     | gamma_soft_mode         |
| FCC-Al-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | harmonic_QHA_applicable |
| FCC-Al-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | 429pro partial, 4/7 volumes done, missing 4,6,7 |                         |
| FCC-Al-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | 429pro partial, 5/7 volumes done, missing 6,7   |                         |
| FCC-Al-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| FCC-Al-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | 429pro partial, 4/7 volumes done, missing 5,6,7 |                         |
| FCC-Cr-Al-4 | FCC_L12 | 4    | 32        | confirmed_branch | 429pro | done 429pro                                     | harmonic_QHA_applicable |
| FCC-Cr-Cr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| FCC-Cr-Hf-4 | FCC_L12 | 4    | 32        | confirmed_branch | B5070  | done B5070                                      | finite_q_soft_mode      |
| FCC-Cr-Nb-4 | FCC_L12 | 4    | 32        | confirmed_branch | 429pro | done 429pro                                     | numerical_unresolved    |
| FCC-Cr-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| FCC-Cr-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| FCC-Cr-Zr-4 | FCC_L12 | 4    | 32        | confirmed_branch | B5070  | queued B5070                                    |                         |
| FCC-Hf-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | harmonic_QHA_applicable |
| FCC-Hf-Cr-4 | FCC_L12 | 4    | 32        | confirmed_branch | B5070  | done B5070                                      | gamma_soft_mode         |
| FCC-Hf-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| FCC-Hf-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | 429pro partial, 5/7 volumes done, missing 1,2   |                         |
| FCC-Hf-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| FCC-Hf-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                   |                         |
| FCC-Hf-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | harmonic_QHA_applicable |
| FCC-Nb-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| FCC-Nb-Cr-4 | FCC_L12 | 4    | 32        | confirmed_branch | B5070  | done B5070                                      | numerical_unresolved    |
| FCC-Nb-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                   |                         |
| FCC-Nb-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| FCC-Nb-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | numerical_unresolved    |
| FCC-Nb-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | queued 429pro                                   |                         |
| FCC-Nb-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                    |                         |
| FCC-Ta-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| FCC-Ta-Cr-4 | FCC_L12 | 4    | 32        | confirmed_branch | B5070  | done B5070                                      | numerical_unresolved    |
| FCC-Ta-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                    |                         |
| FCC-Ta-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| FCC-Ta-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | finite_q_soft_mode      |
| FCC-Ta-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                    |                         |
| FCC-Ta-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                    |                         |
| FCC-Ti-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| FCC-Ti-Cr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | numerical_unresolved    |
| FCC-Ti-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                    |                         |
| FCC-Ti-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                    |                         |
| FCC-Ti-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| FCC-Ti-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                    |                         |
| FCC-Ti-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | harmonic_QHA_applicable |
| FCC-Zr-Al-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| FCC-Zr-Cr-4 | FCC_L12 | 4    | 32        | confirmed_branch | B5070  | done B5070                                      | gamma_soft_mode         |
| FCC-Zr-Hf-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                    |                         |
| FCC-Zr-Nb-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| FCC-Zr-Ta-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | finite_q_soft_mode      |
| FCC-Zr-Ti-4 | FCC_L12 | 4    | 32        | nonmagnetic      | B5070  | queued B5070                                    |                         |
| FCC-Zr-Zr-4 | FCC_L12 | 4    | 32        | nonmagnetic      | 429pro | done 429pro                                     | harmonic_QHA_applicable |
| SER-Al-1    | SER_FCC | 4    | 32        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| SER-Cr-1    | SER_BCC | 2    | 16        | confirmed_branch | B5070  | done B5070                                      | harmonic_QHA_applicable |
| SER-Hf-1    | SER_HCP | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| SER-Nb-1    | SER_BCC | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| SER-Ta-1    | SER_BCC | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| SER-Ti-1    | SER_HCP | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |
| SER-Zr-1    | SER_HCP | 2    | 16        | nonmagnetic      | B5070  | done B5070                                      | harmonic_QHA_applicable |

**66 / 84 done.** 5 more have partial scratch state that `--resume-from` continues.

## Resuming a run that hit the time limit

A job killed at its Slurm limit gets SIGTERM, so `Worker.run_flow` never
reaches its `finally` block and `store_dir` is never created.  The completed
volumes sit in `/nfs_ssd/tmp/<name>-<uuid>` instead.  Continue from there:

```bash
./.conda/bin/python main_dfpt_qha.py -n <endpoint> --slurm --resume \
    --resume-from /nfs_ssd/tmp/<name>-<uuid>
```

`--resume-from` is threaded into the submitted command, so the job itself
sees it.  Without `--resume` it is rejected rather than silently ignored.
