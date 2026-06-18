#!/usr/bin/env bash
# End-to-end pipeline test using downsampled fixtures.
# Usage: bash test/run_test.sh
# Requires: bowtie2, samtools, Rscript, and python3 on PATH (activate the
# custom-alignment-pipeline conda environment first).
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$REPO/test/data"
OUT="$(mktemp -d)"
trap 'rm -rf "$OUT"' EXIT

echo "==> Output dir: $OUT"

INDEX="$OUT/index/barcodes"
BAM_DIR="$OUT/bam"
mkdir -p "$(dirname "$INDEX")" "$BAM_DIR"

BARCODES="$DATA/tf.orf.lib.fa.gz"
ALIGN="$REPO/alignment-barcode-SE-no-reverse-comp.sh"
COMBINE="$REPO/combine_counts.py"
DESEQ2="$REPO/run_deseq2.R"

# Sample definitions: "sample_id  fastq_basename  sorted_bin  replicate"
SAMPLES=(
    "D1_Top  D1-acute-Top-10_S1_L001_R1_001.small.fastq.gz  High  1"
    "D1_Bot  D1-acute-Bot-10_S2_L001_R1_001.small.fastq.gz  Low   1"
    "D2_Top  D2-acute-Top-10_S3_L001_R1_001.small.fastq.gz  High  2"
    "D2_Bot  D2-acute-Bot-10_S4_L001_R1_001.small.fastq.gz  Low   2"
    "D3_Top  D3-acute-Top-10_S5_L001_R1_001.small.fastq.gz  High  3"
    "D3_Bot  D3-acute-Bot-10_S6_L001_R1_001.small.fastq.gz  Low   3"
)

# ── Step 1: Alignment ─────────────────────────────────────────────────────────
echo ""
echo "==> Step 1: Alignment"
first=1
for entry in "${SAMPLES[@]}"; do
    read -r sample fastq condition replicate <<< "$entry"
    echo "    Aligning $sample..."

    if [[ $first -eq 1 ]]; then
        bash "$ALIGN" "$INDEX" "$DATA" "$fastq" "$BAM_DIR/$sample" 1 "$BARCODES" \
            2>&1 | grep -E "alignment rate|Error" || true
        first=0
    else
        bash "$ALIGN" "$INDEX" "$DATA" "$fastq" "$BAM_DIR/$sample" 1 \
            2>&1 | grep -E "alignment rate|Error" || true
    fi

    [[ -s "$BAM_DIR/$sample.counts.txt" ]] \
        || { echo "ERROR: $BAM_DIR/$sample.counts.txt missing or empty"; exit 1; }
done

# ── Step 2: Combine counts + metadata ────────────────────────────────────────
echo ""
echo "==> Step 2: Combining counts and generating metadata"
SHEET="$OUT/sample_sheet.csv"
printf 'sample,sorted_bin,replicate,counts_file\n' > "$SHEET"
for entry in "${SAMPLES[@]}"; do
    read -r sample fastq condition replicate <<< "$entry"
    printf '%s,%s,%s,%s\n' "$sample" "$condition" "$replicate" "$BAM_DIR/$sample.counts.txt" >> "$SHEET"
done

python3 "$COMBINE" "$SHEET" \
    --counts_out  "$OUT/counts_table.csv" \
    --metadata_out "$OUT/metadata.csv"

# ── Step 3: DESeq2 ───────────────────────────────────────────────────────────
echo ""
echo "==> Step 3: DESeq2"
Rscript "$DESEQ2" \
    "$OUT/counts_table.csv" \
    "$OUT/metadata.csv" \
    "$OUT/deseq2_results.tsv"

# ── Assertions ───────────────────────────────────────────────────────────────
echo ""
echo "==> Validating outputs"

N_BARCODES=$(awk 'NR>1' "$OUT/counts_table.csv" | wc -l | tr -d ' ')
N_RESULTS=$(awk 'NR>1' "$OUT/deseq2_results.tsv" | wc -l | tr -d ' ')

[[ "$N_BARCODES" -gt 0 ]] \
    || { echo "ERROR: counts table is empty"; exit 1; }

[[ "$N_RESULTS" -eq "$N_BARCODES" ]] \
    || { echo "ERROR: DESeq2 results rows ($N_RESULTS) != barcode count ($N_BARCODES)"; exit 1; }

echo ""
echo "==> All tests passed"
echo "    Barcodes quantified : $N_BARCODES"
echo "    DESeq2 results rows : $N_RESULTS"
