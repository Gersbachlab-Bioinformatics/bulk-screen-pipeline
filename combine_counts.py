#!/usr/bin/env python3
"""
Combine per-sample Bowtie2/samtools counts files into a single counts table
and generate a DESeq2-compatible metadata file.

Sample sheet CSV format (one row per sample):
    sample,sorted_bin,replicate,counts_file
    D1_Top,High,1,/path/to/D1_Top.counts.txt
    D1_Bot,Low,1,/path/to/D1_Bot.counts.txt
    ...

The sample names must be valid R identifiers (no spaces or special characters).
"""

import argparse
import csv
import sys


def read_counts(filepath):
    counts = {}
    with open(filepath) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                counts[parts[0]] = parts[1]
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sample_sheet",
                        help="CSV with columns: sample, sorted_bin, replicate, counts_file")
    parser.add_argument("--counts_out", required=True,
                        help="Output path for combined counts CSV")
    parser.add_argument("--metadata_out", required=True,
                        help="Output path for DESeq2 metadata CSV")
    args = parser.parse_args()

    # Read sample sheet
    samples = []
    with open(args.sample_sheet) as fh:
        reader = csv.DictReader(fh)
        required = {"sample", "sorted_bin", "replicate", "counts_file"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"Sample sheet is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            samples.append(row)

    if not samples:
        sys.exit("Sample sheet contains no rows.")

    # Read counts files, preserving barcode order from the first sample
    sample_counts = {}
    barcode_order = []
    seen = set()
    for s in samples:
        counts = read_counts(s["counts_file"])
        sample_counts[s["sample"]] = counts
        for barcode in counts:
            if barcode not in seen:
                barcode_order.append(barcode)
                seen.add(barcode)

    # Write combined counts CSV
    sample_names = [s["sample"] for s in samples]
    with open(args.counts_out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["barcode"] + sample_names)
        for barcode in barcode_order:
            row = [barcode] + [sample_counts[name].get(barcode, "0") for name in sample_names]
            writer.writerow(row)

    # Write metadata CSV (row names must match counts table column names)
    with open(args.metadata_out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Sample", "sorted_bin", "replicate"])
        for s in samples:
            writer.writerow([s["sample"], s["sorted_bin"], s["replicate"]])

    print(f"Wrote {len(barcode_order)} barcodes × {len(samples)} samples → {args.counts_out}")
    print(f"Wrote {len(samples)}-sample metadata → {args.metadata_out}")


if __name__ == "__main__":
    main()
