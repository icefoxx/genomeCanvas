# Replicate structure versus calls in yeast real-null splits

ratio = mean across-pair / mean within-pair distance of per-sample meta-profiles (label-free). Null calls: share of regions called with no insertion; FDP: mean over shift/broadening insertions.

| genotype | split | pairing ratio | shape_kernel null calls | shape_any null calls | shift null calls | FDP kernel / any / shift |
|---|---|---|---|---|---|---|
| WT | 1+2 vs rest | 0.96 | 0.00 % | 0.00 % | 0.00 % | 0.081 / 0.059 / 0.050 |
| WT | 1+3 vs rest | 0.93 | 0.00 % | 0.00 % | 0.00 % | 0.026 / 0.046 / 0.055 |
| WT | 1+4 vs rest | 1.13 | 0.00 % | 0.00 % | 0.00 % | 0.047 / 0.033 / 0.025 |
| dKO | 1+2 vs rest | 1.44 | 8.68 % | 5.54 % | 0.00 % | 0.521 / 0.397 / 0.118 |
| dKO | 1+3 vs rest | 0.86 | 0.00 % | 0.00 % | 0.00 % | 0.000 / 0.015 / 0.040 |
| dKO | 1+4 vs rest | 0.83 | 0.00 % | 0.00 % | 0.00 % | 0.000 / 0.012 / 0.010 |
| rtf1-8-9A | 1+2 vs rest | 0.96 | 0.00 % | 0.00 % | 0.00 % | 0.033 / 0.060 / 0.063 |
| rtf1-8-9A | 1+3 vs rest | 1.01 | 0.00 % | 0.00 % | 0.00 % | 0.053 / 0.039 / 0.031 |
| rtf1-8-9A | 1+4 vs rest | 1.03 | 0.00 % | 0.00 % | 0.00 % | 0.031 / 0.033 / 0.027 |
| rtf1-11A | 1+2 vs rest | 1.56 | 0.00 % | 0.00 % | 0.00 % | 0.224 / 0.209 / 0.158 |
| rtf1-11A | 1+3 vs rest | 0.82 | 0.00 % | 0.00 % | 0.00 % | 0.008 / 0.019 / 0.030 |
| rtf1-11A | 1+4 vs rest | 0.82 | 0.00 % | 0.00 % | 0.00 % | 0.004 / 0.008 / 0.011 |
| rtf1 | 1+2 vs rest | 1.66 | 3.14 % | 0.00 % | 0.00 % | 0.488 / 0.363 / 0.100 |
| rtf1 | 1+3 vs rest | 0.81 | 0.00 % | 0.00 % | 0.00 % | 0.003 / 0.018 / 0.022 |
| rtf1 | 1+4 vs rest | 0.79 | 0.00 % | 0.00 % | 0.00 % | 0.000 / 0.006 / 0.014 |
| isw1rtf1 | 1+2 vs rest | 2.28 | 18.28 % | 8.80 % | 0.00 % | 0.602 / 0.515 / 0.097 |
| isw1rtf1 | 1+3 vs rest | 0.72 | 0.00 % | 0.00 % | 0.00 % | 0.000 / 0.014 / 0.018 |
| isw1rtf1 | 1+4 vs rest | 0.72 | 0.00 % | 0.00 % | 0.00 % | 0.000 / 0.014 / 0.018 |

Splits with pairing ratio > 1.2: 4 of 18; their mean `shape_kernel` FDP 0.459, null calls in 3 of them. Other splits: mean FDP 0.020, null calls in 0 of 14. Spearman (ratio, FDP) = 0.91.

