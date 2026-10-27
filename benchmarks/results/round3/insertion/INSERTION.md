# In-silico insertion on real replicates

Commit `5a5be1371e`; uncommitted src files 4; regional.py `15a257546ec0aa00`; src diff `2995aee9d62eb521`; 2026-10-26T21:36:14 → 2026-10-26T21:41:48.


## yeast

| method | null: share called / FWER | amplitude-only: share called | shift20: power / FDP | shift50: power / FDP | shift100: power / FDP | broaden1.3: power / FDP | broaden1.5: power / FDP |
|---|---|---|---|---|---|---|---|
| shape_kernel | 0.0167 / 0.17 | 0.0160 | 0.10 / 0.125 | 0.77 / 0.110 | 0.98 / 0.111 | 0.98 / 0.119 | 0.96 / 0.123 |
| shift | 0.0000 / 0.00 | 0.0000 | 1.00 / 0.051 | 0.88 / 0.051 | 0.34 / 0.054 | 0.20 / 0.046 | 0.21 / 0.045 |
| shape_any | 0.0080 / 0.11 | 0.0050 | 1.00 / 0.099 | 0.93 / 0.099 | 0.97 / 0.101 | 0.98 / 0.106 | 0.96 / 0.111 |
| v2 centroid (moderated t) | 0.0001 / 0.06 | 0.0001 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.056 | 0.00 / 0.056 | 0.00 / 0.056 |
| Welch t on centroid | 0.0000 / 0.00 | 0.0000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 |
| global Wilcoxon over regions | 0.78 of runs significant | 0.78 | 0.78 (share of runs significant) | 0.75 (share of runs significant) | 0.72 (share of runs significant) | 0.72 (share of runs significant) | 0.67 (share of runs significant) |

## proseq

| method | null: share called / FWER | amplitude-only: share called | shift20: power / FDP | shift50: power / FDP | shift100: power / FDP | broaden1.3: power / FDP | broaden1.5: power / FDP |
|---|---|---|---|---|---|---|---|
| shape_kernel | 0.0000 / 0.00 | 0.0000 | 0.36 / 0.010 | 0.91 / 0.013 | 0.99 / 0.014 | 0.12 / 0.005 | 0.40 / 0.007 |
| shift | 0.0000 / 0.00 | 0.0000 | 0.96 / 0.017 | 1.00 / 0.016 | 0.35 / 0.015 | 0.21 / 0.012 | 0.44 / 0.015 |
| shape_any | 0.0000 / 0.00 | 0.0000 | 0.94 / 0.015 | 1.00 / 0.014 | 0.99 / 0.015 | 0.23 / 0.008 | 0.52 / 0.014 |
| v2 centroid (moderated t) | 0.0000 / 0.25 | 0.0000 | 0.00 / 0.250 | 0.00 / 0.250 | 0.43 / 0.008 | 0.00 / 0.125 | 0.00 / 0.104 |
| Welch t on centroid | 0.0000 / 0.00 | 0.0000 | 0.00 / 0.000 | 0.02 / 0.000 | 0.32 / 0.007 | 0.00 / 0.000 | 0.00 / 0.000 |
| global Wilcoxon over regions | 1.00 of runs significant | 1.00 | 1.00 (share of runs significant) | 1.00 (share of runs significant) | 1.00 (share of runs significant) | 0.75 (share of runs significant) | 1.00 (share of runs significant) |
| pausing index (Welch) | 0.0000 / 0.00 | 0.0000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 |
| pausing index (moderated t) | 0.0001 / 0.75 | 0.0005 | 0.00 / 0.625 | 0.04 / 0.074 | 0.27 / 0.043 | 0.00 / 0.438 | 0.00 / 0.277 |

