# Custom alignment pipeline

Scripts to align reads to a custom reference (e.g. barcode sequences, protospacer sequences, etc.) and detect differential abundance. The workflow has two steps: align reads and
count barcodes with `Bowtie2`/`samtools`, then run `DESeq2` on the combined
counts table.

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

```
./alignment-barcode-SE-no-reverse-comp.sh \
    /path/to/barcode_index \
    /path/to/fastq_dir \
    sample_R1.fastq.gz \
    /path/to/output_basename \
    1 \
    /path/to/barcodes.fasta
```

Once the index exists, the FASTA argument can be omitted:

```
./alignment-barcode-SE-no-reverse-comp.sh \
    /path/to/barcode_index \
    /path/to/fastq_dir \
    sample_R1.fastq.gz \
    /path/to/output_basename \
    1
```

This produces a sorted, indexed BAM file, a shuffled BAM for library
complexity estimation, and a `output_basename.counts.txt` table of barcode
counts.

## run_deseq2.R

DESeq2 analysis to detect significant differences in barcode abundance between
two conditions. A combined table of raw counts is expected.

Assuming a `path/to/counts_table.txt` counts table file like this:

```
barcode,barcode_sequence,Low_R1,Low_R2,Low_R3,High_R1,High_R2,High_R3
barcode1,CCTTGTTTCAAATGGATTTT,1723,819,4081,2242,2216,2158
barcode2,TCGAGAAAATCCATTTGAAA,1456,851,3703,1552,2262,2197
barcode3,CGAGAAAATCCATTTGAAAC,902,1226,2840,1475,1859,2178
...
```

And a `path/to/metadata.txt` metadata file like this:

```
Sample,id,sorted_bin,replicate
1,Low R1,Low,1
2,Low R2,Low,2
3,Low R3,Low,3
4,High R1,High,1
5,High R2,High,2
6,High R3,High,3
```

You can run DESeq2 like this:

```
Rscript run_deseq2.R \
    path/to/counts_table.txt \
    path/to/metadata.txt \
    path/to/output.tsv \
    --skip_columns barcode_sequence
```

The `output.tsv` file will contain the DESeq2 results.
