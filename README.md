# Custom alignment pipeline

[![Pipeline test](https://github.com/Gersbachlab-Bioinformatics/bulk-screen-pipeline/actions/workflows/test.yml/badge.svg)](https://github.com/Gersbachlab-Bioinformatics/bulk-screen-pipeline/actions/workflows/test.yml)

Scripts to align reads to a custom reference (e.g. barcode sequences, protospacer sequences, etc.) and detect differential abundance. The workflow has three steps:

1. **`alignment-barcode-SE-no-reverse-comp.sh`** — align reads and count barcodes per sample with `Bowtie2`/`samtools`
2. **`combine_counts.py`** — merge per-sample count files into a single counts table and generate a `DESeq2`-compatible metadata file
3. **`run_deseq2.R`** — run `DESeq2` on the combined counts table to detect differential barcode abundance

## Dependencies

Use the provided `environment.yml` to create a conda environment with all required tools:

```bash
conda env create -f environment.yml
conda activate custom-alignment-pipeline
```

## alignment-barcode-SE-no-reverse-comp.sh

Alignment shell script to run `Bowtie2` and compute per-barcode counts with
`samtools`. Reads are mapped single-end in `--end-to-end --very-sensitive`
mode against a custom barcode reference, without reverse-complement mapping
(`--norc`).

Because the sequenced reads are longer than the barcodes, the number of bases
to trim from the 3' end is exposed as a parameter (`read_length - barcode_length`),
so only the barcode portion of each read is aligned. For example, 21 bp reads
with 20 nt barcodes would use a trim of `1`.

### Parameters

```
$1: Reference genome location (Bowtie2 index basename)
$2: Folder with FASTQ files
$3: Sample Read 1 FASTQ file
$4: Output basename
$5: Number of bp to trim from the 3' end (read_length - barcode_length)
$6: Barcodes FASTA file [Optional, only required the first time to build the index]
```

### Usage

The first time, pass the barcodes FASTA file so the Bowtie2 index is built:

```bash
./alignment-barcode-SE-no-reverse-comp.sh \
    /path/to/barcode_index \
    /path/to/fastq_dir \
    sample_R1.fastq.gz \
    /path/to/output_basename \
    1 \
    /path/to/barcodes.fasta
```

Once the index exists, the FASTA argument can be omitted:

```bash
./alignment-barcode-SE-no-reverse-comp.sh \
    /path/to/barcode_index \
    /path/to/fastq_dir \
    sample_R1.fastq.gz \
    /path/to/output_basename \
    1
```

This produces a sorted, indexed BAM file, a shuffled BAM for library
complexity estimation, and a `output_basename.counts.txt` table of barcode
counts (two tab-separated columns: barcode name and read count).

## combine_counts.py

Merges per-sample `counts.txt` files (produced by the alignment script) into a
single counts table and generates a `DESeq2`-compatible metadata file.

### Sample sheet format

Provide a CSV with one row per sample:

```
sample,sorted_bin,replicate,counts_file
D1_Top,High,1,/path/to/D1_Top.counts.txt
D1_Bot,Low,1,/path/to/D1_Bot.counts.txt
D2_Top,High,2,/path/to/D2_Top.counts.txt
D2_Bot,Low,2,/path/to/D2_Bot.counts.txt
```

Sample names must be valid R identifiers (no spaces or special characters).

### Usage

```bash
python combine_counts.py sample_sheet.csv \
    --counts_out counts.csv \
    --metadata_out metadata.csv
```

This produces:

- `counts.csv` — barcodes × samples counts table, ready for `DESeq2`
- `metadata.csv` — sample metadata matched to the counts table columns

## run_deseq2.R

DESeq2 analysis to detect significant differences in barcode abundance between
two conditions (sorted bins). Expects the counts table and metadata produced by
`combine_counts.py`.

`counts.csv` format:

```
barcode,D1_Top,D1_Bot,D2_Top,D2_Bot
barcode1,1723,819,4081,2242
barcode2,1456,851,3703,1552
barcode3,902,1226,2840,1475
...
```

`metadata.csv` format:

```
Sample,sorted_bin,replicate
D1_Top,High,1
D1_Bot,Low,1
D2_Top,High,2
D2_Bot,Low,2
```

### Usage

```bash
Rscript run_deseq2.R \
    counts.csv \
    metadata.csv \
    results.tsv
```

The `results.tsv` file will contain the DESeq2 results (log2 fold change,
p-value, adjusted p-value, etc.) for each barcode.
