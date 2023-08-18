# Description: This script generates the samples.csv file needed per each cellranger run in multiome datasets
# Author: Coral Fustero-Torre
# email: coral.fustero@bih-charite.de

# -- Read input parameters -- #
out <- snakemake@params[["outdir"]]
sample_ID <- snakemake@params[["sample"]]


# -- Function -- #
per_sample<- function(x, outdir, sample_ID){
  GEX <- data.frame(fastqs = paste0(x[["fqs"]],"/GEX"),
                    sample = paste0(x[["sample"]],"_GEX"),
                    library_type = "Gene Expression") 
  ATAC <- data.frame(fastqs = paste0(x[["fqs"]],"/ATAC"),
                sample = paste0(x[["sample"]],"_ATAC"),
                library_type = "Chromatine Accesibility") 
  out <- rbind(GEX, ATAC)
  message(paste0("Saving individual sample.csv file for sample:", sample_ID))
  write.csv(out, file.path(outdir, sample_ID, "cellranger_count", "samples.csv"))
}

# -- Read samples file -- #
message("Reading units.tsv file")
if(file.exists("config/units.tsv")){
    units_file <- read.table("config/units.tsv", sep = "\t", header = T,  check.names=FALSE)  
} else {
    message("The units.tsv couldn't be found. Please check it is in the right location.")
}

# -- Select row -- #
message(paste0("Performing sample matching in units.tsv for sample ", sample_ID))
position <- match(sample_ID, units_file$sample)
# -- Generate samples.csv -- #
message("Generating per sample samples.csv file")
per_sample(x = units_file[position,], outdir = out, sample_ID = sample_ID)

                                                                                   
