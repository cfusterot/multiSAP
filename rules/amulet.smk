import glob

rule atac_bc:
    input:
        per_barcode_metrics="{}/{{sample}}/cellranger_count/outs/per_barcode_metrics.csv".format(OUTDIR)
    output:
        atac_bc="{}/{{sample}}/cellranger_count/outs/atac_bc.csv".format(OUTDIR)
    threads: get_resource("default", "threads")
    log:
        err="{}/{{sample}}/atac_bc.err".format(LOGDIR),
        out="{}/{{sample}}/atac_bc.out".format(LOGDIR),
    resources:
        mem_mb=get_resource("default", "mem_mb"),
        walltime=get_resource("default", "walltime")
    shell:
        """ 
        awk -F"\\"*,\\"*" '{{print $1sep$4 ;sep=","}}' {input.per_barcode_metrics} | sed -e "1s/barcodeis_cell/barcode,is__cell_barcode/" > {output.atac_bc}  
        """


rule amulet:
    input:
        atac_bc="{}/{{sample}}/cellranger_count/outs/atac_bc.csv".format(OUTDIR)
    output:
        multiplets="{}/{{sample}}/amulet/MultipletSummary.txt".format(OUTDIR)
    conda:
        "../envs/amulet.yaml"
    params:
        fragments="{}/{{sample}}/cellranger_count/outs/atac_fragments.tsv.gz".format(OUTDIR),
        amulet="resources/AMULET",
        autosomes=config['amulet']['autosomes'],
        blacklist=config['amulet']['blacklist']
    threads: get_resource("amulet", "threads")
    resources:
        mem_mb=get_resource("amulet", "mem_mb"),
        walltime=get_resource("amulet", "walltime")
    log:
        err="{}/{{sample}}/amulet.err".format(LOGDIR),
        out="{}/{{sample}}/amulet.out".format(LOGDIR)
    shell:
        """
        chmod +x {params.amulet}/AMULET.sh &&
        {params.amulet}/AMULET.sh {params.fragments} {input.atac_bc} \
        {params.autosomes} {params.blacklist} \
        {OUTDIR}/{wildcards.sample}/amulet {params.amulet}
        """
