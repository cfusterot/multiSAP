import glob

rule cellbender:
    input:
        h5="{}/{{sample}}/cellranger_count/outs/raw_feature_bc_matrix.h5".format(OUTDIR)
    output:
        cb="{}/{{sample}}/cellbender/cellbender_output_file.h5".format(OUTDIR)
    conda:
        "../envs/cellbender.yaml"
    threads: get_resource("cellbender", "threads")
    resources:
        time="12:00:00",
        mem=60000,
        nodes=1,
        nvidia_gpu=1,
        partition="gpu",
        slurm="gres=gpu:tesla:1",
        walltime=get_resource("cellbender", "walltime")
    log:
        err="{}/{{sample}}/cellbender.err".format(LOGDIR),
        out="{}/{{sample}}/cellbender.out".format(LOGDIR)
    shell:
        """
        cellbender remove-background --cuda --input {input.h5} --output {output.cb} --expected-cells 15000 --total-droplets-included 50000 --fpr 0.01 --epochs 150
        """
