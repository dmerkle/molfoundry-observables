# Read retained A/B benchmark measurements. This does not rerun either benchmark or establish equality of all composite expressions.

import csv
from pathlib import Path

package = Path(__file__).resolve().parents[1]
summary = package / "notebooks/visual/evidence/ab-scaling/summary.csv"
with summary.open(encoding="utf-8") as stream:
    measurements = list(csv.DictReader(stream))

for row in measurements:
    print(row["word"], "MolFoundry seconds:", row["molfoundry_median_seconds"],
          "FRAGGER seconds:", row["fragger_median_seconds"])

result = {"measurements": measurements}
