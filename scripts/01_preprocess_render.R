# Snakemake params 
outs = snakemake@input$outs
print(outs)
directory = snakemake@params$directory
print(directory)
sample_id = unlist(sapply(strsplit(unlist(sapply(strsplit(outs, directory),`[`, 2, simplify=FALSE)), "/cellranger_count/outs"), `[`, 1, simplify=FALSE))
reference = snakemake@params$reference
ncount_atac_max = snakemake@params$ncount_atac_max
ncount_rna_max = snakemake@params$ncount_rna_max
ncount_atac_min = snakemake@params$ncount_atac_min
ncount_rna_min = snakemake@params$ncount_rna_min
nuclesome_signal = snakemake@params$nuclesome_signal
tss_enrichment = snakemake@params$tss_enrichment
percent_mt = snakemake@params$percent_mt
min_depth = snakemake@params$min_depth
n_cells_conf_detected = snakemake@params$n_cells_conf_detected
strand_correlation = snakemake@params$strand_correlation
min_cell_var = snakemake@params$min_cell_var

message(paste0("Analysing sample ", sample_id))

message("Rendering analysis report: 01 sample pre-processing")
render_report = function(sample) {
  rmarkdown::render('01_pre-process.Rmd',
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
                min_cell_var = min_cell_var),
  output_file = file.path(directory, "signac", paste0('01_preprocessing_', sample_id,'.html'))
)}


