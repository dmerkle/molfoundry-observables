# A/B composition measurements

Native Linux x86_64 on `algochem-compute`, Intel Xeon w7-2595X.
MolFoundry 0.1.9 / Python 3.12.3, measured 14 October 2026;
FRAGGER `0d9595c` / Scala 2.13.2 / Java 11 with `-Xmx4g`, measured 9 October 2026.

`summary.csv` contains medians of five repetitions after two per-case warm-ups. Each workload makes twelve ordered compositions for six terminal-extension rules. Startup and compilation are excluded.

The comparison covers even n=2–28. The FRAGGER samples are retained unchanged;
all twelve ordered counts agree with the current MolFoundry results at each size.

[Protocol and code](../../../../scripts/comparison/ab_scaling/README.md).
