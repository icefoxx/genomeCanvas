# Mouse 2-cell reference demo (GRCm38 / mm10)

Only the small cytoband file is stored in git. The large inputs (about 104 MB) are downloaded and
checked against their SHA-256 with:

```bash
python tests/fixtures/reference_demo/download.py
```

| File | Source |
|---|---|
| `GSM5746910_MS_2cell_Input_rep1.bigWig` | GEO GSM5746910 supplementary file |
| `GSM5746911_MS_2cell_Input_rep2.bigWig` | GEO GSM5746911 supplementary file |
| `GSM5746912_MS_2cell_IP_rep1.bigWig` | GEO GSM5746912 supplementary file |
| `GSM5746913_MS_2cell_IP_rep2.bigWig` | GEO GSM5746913 supplementary file |
| `Mus_musculus.GRCm38.102.gtf.gz` | Ensembl release 102, `pub/release-102/gtf/mus_musculus/` |
| `mm10_cytoBandIdeo.txt.gz` | UCSC mm10 `cytoBandIdeo` (in git) |

The test suite needs only `mm10_cytoBandIdeo.txt.gz`. The downloaded files are used by
`examples/reference_comparison/` and `examples/deeptools_realdata/`.
