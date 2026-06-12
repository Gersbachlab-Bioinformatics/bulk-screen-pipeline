#!/usr/bin/env Rscript

## Load appropriate libraries in R without startup messages
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(argparse))

## Function to read data with specified separator
read_data <- function(file_path, sep) {
  if (sep == ",") {
    return(read.csv(file_path, row.names = 1))
  } else {
    return(read.table(file_path, row.names = 1, sep = sep))
  }
}

## Function to remove specified columns from a data frame
remove_columns <- function(data, columns) {
  if (!is.null(columns)) {
    columns <- unlist(strsplit(columns, ","))
    data <- data[, !(colnames(data) %in% columns)]
  }
  return(data)
}

## Create a parser object
parser <- ArgumentParser(description = "Run DESeq2 analysis")

## Add arguments
parser$add_argument("counts_table_path", help = "Path to the counts table file")
parser$add_argument("metadata_path", help = "Path to the metadata file")
parser$add_argument("output_path", help = "Path to the output file")
parser$add_argument("--sep", default = ",", help = "File separator delimiter (default is comma)")
parser$add_argument("--skip_columns", default = NULL, help = "Comma-separated list of columns to skip in the counts table")

## Parse the command line arguments
args <- parser$parse_args()

## Load the data
counts_table <- read_data(args$counts_table_path, args$sep)
metadata <- read_data(args$metadata_path, args$sep)

## Remove specified columns from the counts table
counts_table <- remove_columns(counts_table, args$skip_columns)

## Generate a DeSeq object by inputting the count table, metadata, and specifying the design for the differential analysis
## Here the design is to compare barcode abundance across sorted bins while keeping replicates paired

metadata$replicate <- factor(metadata$replicate)
dds <- DESeqDataSetFromMatrix(countData = counts_table, colData = metadata, design = ~sorted_bin + replicate, tidy = FALSE)
dds <- DESeq(dds)

## Specify the groups you are contrasting - here we are comparing the counts of barcodes in low to high

res_hl <- results(dds, contrast = c("sorted_bin", "High", "Low"))
High_low <- as.data.frame(res_hl)

## Save the results to a tab-delimited file without quotes
write.table(High_low, file = args$output_path, sep = "\t", quote = FALSE, row.names = TRUE)
