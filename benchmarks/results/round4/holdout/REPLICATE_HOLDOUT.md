# Held-out-region check: pairing ratio on half A vs null calls on half B (yeast)

Commit `5a5be1371e`; regional.py `28584ac7dd480341`; src diff `ede132e32bdbd65b`; seed 9600000; 1796 / 1796 regions; 2026-10-26T22:51:54 → 2026-10-26T22:52:03.

| split | pairing ratio (A) | null calls on B: shape_kernel / shape_any / shift |
|---|---|---|
| isw1rtf1:isw1rtf1_1+isw1rtf1_2 | 1.90 | 23.06 % / 12.12 % / 0.00 % |
| rtf1:rtf1_1+rtf1_2 | 1.51 | 3.44 % / 3.18 % / 0.00 % |
| dKO:dKO_1+dKO_2 | 1.41 | 9.23 % / 6.15 % / 0.00 % |
| rtf1-11A:rtf1-11A_1+rtf1-11A_2 | 1.35 | 0.00 % / 0.00 % / 0.00 % |
| WT:WT_1+WT_4 | 1.09 | 0.00 % / 0.00 % / 0.00 % |
| rtf1-8-9A:rtf1-8-9A_1+rtf1-8-9A_3 | 1.00 | 0.00 % / 0.00 % / 0.00 % |
| rtf1-8-9A:rtf1-8-9A_1+rtf1-8-9A_2 | 1.00 | 0.00 % / 0.00 % / 0.00 % |
| rtf1-8-9A:rtf1-8-9A_1+rtf1-8-9A_4 | 1.00 | 0.00 % / 0.00 % / 0.00 % |
| WT:WT_1+WT_3 | 0.96 | 0.00 % / 0.00 % / 0.00 % |
| WT:WT_1+WT_2 | 0.96 | 0.00 % / 0.00 % / 0.00 % |
| rtf1-11A:rtf1-11A_1+rtf1-11A_4 | 0.89 | 0.00 % / 0.00 % / 0.00 % |
| dKO:dKO_1+dKO_3 | 0.87 | 0.00 % / 0.00 % / 0.00 % |
| rtf1-11A:rtf1-11A_1+rtf1-11A_3 | 0.85 | 0.00 % / 0.00 % / 0.00 % |
| rtf1:rtf1_1+rtf1_3 | 0.85 | 0.00 % / 0.00 % / 0.00 % |
| dKO:dKO_1+dKO_4 | 0.84 | 0.00 % / 0.00 % / 0.00 % |
| rtf1:rtf1_1+rtf1_4 | 0.81 | 0.00 % / 0.00 % / 0.00 % |
| isw1rtf1:isw1rtf1_1+isw1rtf1_3 | 0.77 | 0.00 % / 0.00 % / 0.00 % |
| isw1rtf1:isw1rtf1_1+isw1rtf1_4 | 0.76 | 0.00 % / 0.00 % / 0.00 % |

Pre-specified predictions: (1) higher null-call share on B for ratio(A) > 1.2: **met** (8.93 % vs 0.00 %; splits with calls 3 of 4 vs 0 of 14); (2) Spearman > 0.5: **met** (0.51).

