# Description: This script generates the samples.csv file needed per each cellranger run in multiome datasets
# Author: Coral Fustero-Torre
# email: coral.fustero@bih-charite.de

# -- Read input parameters -- #
path <- snakemake@params[["input_dir"]]
sample_ID <- snakemake@params[["sample_id"]]
out <- snakemake@params[["out"]]
samples <- snakemake@params[["samples"]]
input_format <- snakemake@params[["input_format"]]
sample_type <- snakemake@params[["sample_type"]]

# -- Read samples file -- #
message("Reading samples.tsv")
samples <- read.table(samples, sep = "\t", header = T,  check.names=FALSE)
