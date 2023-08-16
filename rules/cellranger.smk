import glob
import os
from pathlib import Path

def get_sample_option(wc):
    option_str=""
    prefix=units["sample"][wc.sample] 
    option_str += f"--sample={prefix}"
    return option_str

def get_fq_path(wc):
    filepath = units["fq1"][wc.sample][0]
    dirname, basename = os.path.split(filepath)
    return dirname

rule samples_csv:
    output:
        csv = "{}/{{sample}}/cellranger_count/samples.csv".format(OUTDIR)
    params:
        outdir = OUTDIR,
        sample = lambda wc: units["fqs"][wc.sample]
    threads:
        get_resource("default", "threads")
    resources:
        mem_mb=get_resource("default", "mem_mb"),
        walltime=get_resource("default", "walltime")
    log:
        err="{}/{{sample}}/samples_csv.err".format(LOGDIR),
        out="{}/{{sample}}/samples_csv.out".format(LOGDIR),
        time="{}/{{sample}}/samples_csv.time".format(LOGDIR)
    shell:
        "../scripts/create_samples_csv.R"

rule cellranger_count:
    input:
        fq=lambda wc: units["fqs"][wc.sample],
        csv = "{}/{{sample}}/cellranger_count/samples.csv".format(OUTDIR)
    output:
        finish="{}/{{sample}}/cellranger_count/cellranger.finish".format(OUTDIR)
    params:    
        reference=config['cellranger']['reference']
    envmodules:
        config['envmodules']['cellranger']
    threads: get_resource("cellranger", "threads")
    resources:
        mem_mb=get_resource("cellranger", "mem_mb"),
        walltime=get_resource("cellranger", "walltime")
    log:
        err="{}/{{sample}}/cellranger.err".format(LOGDIR),
        out="{}/{{sample}}/cellranger.out".format(LOGDIR),
        time="{}/{{sample}}/cellranger.time".format(LOGDIR)
    shell:
        """
        {DATETIME} > {log.time} &&
        cellranger-arc count --id={wildcards.sample} \
        --reference={params.reference} \
        --libraries={samples.csv} 2> {log.err} > {log.out}
        mv {wildcards.sample}/* {OUTDIR}/{wildcards.sample}/cellranger_count
        touch {OUTDIR}/{wildcards.sample}/cellranger_count/cellranger.finish
        rm -r {wildcards.sample}
        {DATETIME} >> {log.time} 
        """