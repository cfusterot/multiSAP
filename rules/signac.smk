import glob

if config["signac"]["enable"]:
    rule step1_preprocess:
        input:
            outs="{}/{{sample}}/cellranger_count/outs/".format(OUTDIR), 
            mgatk="{}/{{sample}}/mgatk/final".format(OUTDIR), 
            amulet="{}/{{sample}}/amulet/MultipletBarcodes_01.txt".format(OUTDIR)
        output:
            directory=expand("{OUTDIR}/signac", OUTDIR = OUTDIR),
            qc_plots=expand("{OUTDIR}/signac/plots/vlnplot_qc_beforefiltering_ATAC.pdf", OUTDIR = OUTDIR)
        conda:
            "../envs/signac.yaml"
        params: 
            min_peak_width=config['signac']['qc']['min_peak_width'],
            max_peak_width=config['signac']['qc']['max_peak_width'],
            min_counts=config['signac']['qc']['min_counts'],
            integration=config['signac']['integrate_samples'],
            nb_cores=config['signac']['cores'],
            granges_ens=config['signac']['annotation']['granges_ens'],
            genome=config['signac']['annotation']['genome']
        threads: get_resource("signac", "threads")
        resources:
            mem_mb=get_resource("signac", "mem_mb"),
            walltime=get_resource("signac", "walltime")
        log:
            err=expand("{LOGDIR}/signac/step1_signac.err", LOGDIR = LOGDIR),
            out=expand("{LOGDIR}/signac/step1_signac.out", LOGDIR = LOGDIR)
        # Add here something 
        script:
            "../scripts/01_pre-process.Rmd"
