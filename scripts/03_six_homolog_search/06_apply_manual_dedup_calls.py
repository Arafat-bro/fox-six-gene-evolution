"""
Apply manual resolutions for the two ambiguous scaffold-duplication cases
(see docs/decisions/manual_duplication_calls.md) on top of the automated
dedup output.

SIX3 / GCA_003615085.1: copy2 (MRCU01000019.1) excluded as a duplicate
    scaffold representation of copy1 (MRCU01000016.1) -- sequences
    confirmed byte-identical by direct extraction.
SIX7 / GCA_013423245.1: no change -- both copies already retained by the
    automated script (correctly classified as independent_copy).
"""
import csv
from pathlib import Path

IN_FILE = Path("results/six_homologs/filtered/genuine_homologs_deduped.tsv")
OUT_FILE = Path("results/six_homologs/filtered/genuine_homologs_deduped_final.tsv")

MANUALLY_EXCLUDED_COPY_IDS = {
    "SIX3_GCA_003615085.1_copy2",
}

with open(IN_FILE) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)
    fieldnames = reader.fieldnames

kept_rows = [r for r in rows if r["copy_id"] not in MANUALLY_EXCLUDED_COPY_IDS]

with open(OUT_FILE, "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(kept_rows)

print(f"Input rows: {len(rows)}")
print(f"Manually excluded: {len(rows) - len(kept_rows)}")
print(f"Final rows written to {OUT_FILE}: {len(kept_rows)}")
