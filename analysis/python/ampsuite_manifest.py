#!/usr/bin/env python3
import argparse, csv, pathlib, sys

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tumor",  required=True, help="space-separated BAMs")
    ap.add_argument("--normal", required=True, help="space-separated BAMs")
    ap.add_argument("--output", required=True, help="output TSV path")
    return ap.parse_args()

def rows(bam_list: list[str], sample_type: str):
    for bam in bam_list:
        if not pathlib.Path(bam).exists():
            sys.stderr.write(f"[AmpSuite]  skip  {bam}  (missing)\n")
            continue
        aliquot = pathlib.Path(bam).name.split(".")[0]
        yield [aliquot, pathlib.Path(bam).resolve(), sample_type]

def main():
    args = parse_args()
    tumors  = args.tumor.split()
    normals = args.normal.split()

    out = pathlib.Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerows(rows(tumors,  "tumor"))
        writer.writerows(rows(normals, "normal"))

if __name__ == "__main__":
    main()
