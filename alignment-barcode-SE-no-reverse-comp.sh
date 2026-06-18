#! /bin/bash

## command line parameters
# $1: Reference genome location (Bowtie2 index basename)
# $2: Folder with FASTQ files
# $3: Sample Read 1 FASTQ file
# $4: Output basename
# $5: Number of bp to trim from the 3' end of each read
#     (i.e. read_length - barcode_length, so only the barcode is aligned)
#     Defaults to 1 if not provided. Use 0 when reads are exactly barcode length,
#     but note that with --end-to-end bowtie2 will penalise any overhanging bases,
#     so alignment rates will drop significantly if reads are longer than barcodes.
# $6: Barcodes FASTA file [Optional, only required the first time to build the index]

## dependencies:
# - samtools >=1.3.1
# - bowtie2 >=2.3.5.1

TRIM=${5:-1}

# Build Bowtie2 custom reference index if it doesn't exist
if [[ ! -e $1.1.bt2 ]];
then
	bowtie2-build $6 $1
fi


perform_alignment ()
{

# Map single-end reads in end-to-end --very-sensitive mode.
# The barcodes are shorter than the reads, so we trim the extra bases from the
# 3' end ($5) and align the remaining barcode-length sequence end-to-end.
bowtie2 \
	--very-sensitive \
	--end-to-end \
	--threads 32 `# Using 32 CPUs, adjust as needed` \
	--norc `# Do not try to map reverse complements` \
	-3 $5 `# Trim $5 bp from the 3' end so only the barcode is aligned` \
	-I 0 `# no minimum length alignment` \
	-X 200 `# to reveal chimeric alignments, but faster than the 500bp default` \
	-x $1 `# Bowtie2 index basename` \
	-U $2/$3 `# FASTQ file with reads to align` \
| samtools view -bS - > $4.bam `# Save BAM file` \
&& samtools sort \
	-@ 32 `# Using 32 CPUs, adjust as needed` \
	-m 2G `# Memory per thread, adjust as needed` \
	$4.bam \
	-o $4.sorted.bam  \
&& samtools index $4.sorted.bam \
&& samtools idxstats `# Extract count table for each barcode in the reference` \
	$4.sorted.bam \
	| awk '{print $1"\t"$3}' \
	| grep -v "^\*" \
> $4.counts.txt

# Create shuffled BAM file for library complexity estimation
samtools bamshuf $4.sorted.bam $4.shuffle

}

perform_alignment $1 $2 $3 $4 $TRIM
