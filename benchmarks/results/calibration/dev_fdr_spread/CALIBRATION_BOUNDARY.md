# Calibration boundaries of the per-region shape tests (nominal FDR 0.05)

Generated from `calibration_boundary.json` of the same run. Commit `5a5be1371e` (profileseg-v3); uncommitted changes under src/: 4; regional.py sha256 `15a257546ec0aa00`; src diff sha256 `2995aee9d62eb521`; options {"reps": 100, "null": "strata", "fdr_strata": "spread", "seed_base": 5000000, "factors": "mixed,mixed_large"}; run 2026-10-26T20:58:23 → 2026-10-26T21:01:46.

FDP: mean [95 % Monte Carlo CI] / power; at 0 % changed the FWER [exact binomial CI]. Verdict: EXCEEDS = CI above 0.05; WATCH = estimate above 0.075; ok = no clear exceedance (not proof of FDR ≤ 0.05). p01 = share of unchanged regions with p < 0.01 (nominal 0.01).

| factor | level | design | class | shape_kernel: FDP / power / p01 | shift: FDP / power / p01 | shape_any: FDP / power / p01 | verdict |
|---|---|---|---|---|---|---|---|
| mixed | 3 classes | 2v2 | all | 0.045 [0.037, 0.053] / 0.34 / 0.0116 | 0.032 [0.023, 0.041] / 0.16 / 0.0102 | 0.035 [0.028, 0.042] / 0.31 / 0.0113 | ok |
| mixed | 3 classes | 2v2 | class0 | 0.039 [0.031, 0.046] / 0.83 / 0.0225 | 0.002 [0.000, 0.004] / 0.40 / 0.0023 | 0.020 [0.014, 0.025] / 0.76 / 0.0131 | ok |
| mixed | 3 classes | 2v2 | class1 | 0.000 [0.000, 0.000] / 0.00 / 0.0000 | 0.160 [0.097, 0.223] / 0.01 / 0.0109 | 0.162 [0.096, 0.227] / 0.01 / 0.0054 | shift:EXCEEDS; shape_any:EXCEEDS |
| mixed | 3 classes | 2v2 | class2 | 0.075 [0.055, 0.096] / 0.18 / 0.0123 | 0.087 [0.058, 0.116] / 0.07 / 0.0174 | 0.066 [0.048, 0.084] / 0.16 / 0.0152 | shape_kernel:EXCEEDS; shift:EXCEEDS |
| mixed | 3 classes | 3v3 | all | 0.047 [0.041, 0.054] / 0.42 / 0.0122 | 0.033 [0.026, 0.040] / 0.22 / 0.0097 | 0.041 [0.035, 0.047] / 0.41 / 0.0110 | ok |
| mixed | 3 classes | 3v3 | class0 | 0.033 [0.028, 0.039] / 0.94 / 0.0181 | 0.007 [0.003, 0.011] / 0.49 / 0.0060 | 0.022 [0.017, 0.027] / 0.91 / 0.0124 | ok |
| mixed | 3 classes | 3v3 | class1 | 0.000 [0.000, 0.000] / 0.01 / 0.0000 | 0.167 [0.107, 0.227] / 0.03 / 0.0117 | 0.163 [0.105, 0.222] / 0.04 / 0.0054 | shift:EXCEEDS; shape_any:EXCEEDS |
| mixed | 3 classes | 3v3 | class2 | 0.085 [0.068, 0.102] / 0.31 / 0.0186 | 0.077 [0.052, 0.102] / 0.13 / 0.0115 | 0.073 [0.058, 0.089] / 0.29 / 0.0152 | shape_kernel:EXCEEDS; shift:EXCEEDS; shape_any:EXCEEDS |
| mixed_large | 3 classes × 1000 | 2v2 | all | 0.059 [0.055, 0.062] / 0.34 / 0.0116 | 0.032 [0.029, 0.036] / 0.16 / 0.0106 | 0.048 [0.044, 0.051] / 0.32 / 0.0111 | shape_kernel:EXCEEDS |
| mixed_large | 3 classes × 1000 | 2v2 | class0 | 0.048 [0.045, 0.051] / 0.85 / 0.0145 | 0.018 [0.015, 0.021] / 0.44 / 0.0077 | 0.036 [0.032, 0.039] / 0.81 / 0.0115 | ok |
| mixed_large | 3 classes × 1000 | 2v2 | class1 | 0.008 [0.000, 0.020] / 0.00 / 0.0089 | 0.000 [0.000, 0.000] / 0.00 / 0.0109 | 0.000 [0.000, 0.000] / 0.00 / 0.0093 | ok |
| mixed_large | 3 classes × 1000 | 2v2 | class2 | 0.108 [0.096, 0.120] / 0.17 / 0.0115 | 0.143 [0.122, 0.164] / 0.05 / 0.0132 | 0.110 [0.100, 0.120] / 0.14 / 0.0126 | shape_kernel:EXCEEDS; shift:EXCEEDS; shape_any:EXCEEDS |
| mixed_large | 3 classes × 1000 | 3v3 | all | 0.044 [0.041, 0.046] / 0.43 / 0.0113 | 0.033 [0.029, 0.036] / 0.19 / 0.0108 | 0.041 [0.038, 0.044] / 0.40 / 0.0114 | ok |
| mixed_large | 3 classes × 1000 | 3v3 | class0 | 0.039 [0.036, 0.041] / 0.96 / 0.0125 | 0.023 [0.019, 0.026] / 0.49 / 0.0097 | 0.034 [0.031, 0.036] / 0.93 / 0.0120 | ok |
| mixed_large | 3 classes × 1000 | 3v3 | class1 | 0.014 [0.000, 0.029] / 0.01 / 0.0096 | 0.000 [0.000, 0.000] / 0.00 / 0.0110 | 0.000 [0.000, 0.000] / 0.00 / 0.0103 | ok |
| mixed_large | 3 classes × 1000 | 3v3 | class2 | 0.059 [0.054, 0.064] / 0.31 / 0.0117 | 0.085 [0.072, 0.098] / 0.09 / 0.0117 | 0.067 [0.060, 0.074] / 0.26 / 0.0120 | shape_kernel:EXCEEDS; shift:EXCEEDS; shape_any:EXCEEDS |
