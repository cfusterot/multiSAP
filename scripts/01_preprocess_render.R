# Snakemake params 
outs = snakemake@params$outs
print(outs)
directory = snakemake@params$directory
print(directory)
sample_id = unlist(sapply(strsplit(unlist(sapply(strsplit(outs, directory),`[`, 2, simplify=FALSE)), "/cellranger_count/outs"), `[`, 1, simplify=FALSE))
print(sample_id)
reference = snakemake@params$reference
print(reference)
ncount_atac_max = snakemake@params$ncount_atac_max
print(ncount_atac_max)
ncount_rna_max = snakemake@params$ncount_rna_max
print(ncount_rna_max)
ncount_atac_min = snakemake@params$ncount_atac_min
print(ncount_atac_min)
ncount_rna_min = snakemake@params$ncount_rna_min
print(ncount_rna_min)
nuclesome_signal = snakemake@params$nuclesome_signal
print(nuclesome_signal)
tss_enrichment = snakemake@params$tss_enrichment
print(tss_enrichment)
percent_mt = snakemake@params$percent_mt
print(percent_mt)
min_depth = snakemake@params$min_depth
print(min_depth)
n_cells_conf_detected = snakemake@params$n_cells_conf_detected
print(n_cells_conf_detected)
strand_correlation = snakemake@params$strand_correlation
print(strand_correlation)
min_cell_var = snakemake@params$min_cell_var
print(min_cell_var)
haplogroup = snakemake@params$haplogroup
print(haplogroup)
message(paste0("Analysing sample ", sample_id))
message(paste0("File path: ", file.path(directory, sample_id,"signac", paste0('01_preprocessing_', sample_id,'.html'))))

message("Rendering analysis report: 01 sample pre-processing")
#render_report = function(sample) {
  rmarkdown::render('scripts/01_pre-process.Rmd',
  params = list(sample = sample_id, 
                outs = outs, 
                directory = directory,
                reference = reference, 
                ncount_atac_max = ncount_atac_max, 
                ncount_rna_max = ncount_rna_max, 
                ncount_atac_min = ncount_atac_min, 
                ncount_rna_min = ncount_rna_min, 
                nuclesome_signal = nuclesome_signal, 
                tss_enrichment = tss_enrichment, 
                percent_mt = percent_mt, 
                min_depth = min_depth, 
                n_cells_conf_detected = n_cells_conf_detected, 
                strand_correlation = strand_correlation, 
                min_cell_var = min_cell_var,
                haplogroup = haplogroup),
  output_file = file.path(directory, sample_id,"signac", paste0('01_preprocessing_', sample_id,'.html')))
#}


