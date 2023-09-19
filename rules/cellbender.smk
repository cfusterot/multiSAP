import glob
if config["cellbender"]["enable"]:
    rule cellbender:
        input:
            h5="{}/{{sample}}/cellranger_count/outs/raw_feature_bc_matrix.h5".format(OUTDIR)
        output:
            cb="{}/{{sample}}/cellbender/cellbender_output_file.h5".format(OUTDIR)
        conda:
            "../envs/cellbender.yaml"
        threads: get_resource("cellbender", "threads")
        params:
            fpr = config['cellbender']['params']['fpr'],
            epochs = config['cellbender']['params']['epochs'],
            expected_cells = config['cellbender']['params']['expected_cells'],
            extra = config['cellbender']['params']['extra']             
        resources:
            time="4:00:00",
            partition='gpu',
            nvidia_gpu=1,
            gpu='gpu:tesla:1'
        log:
            err="{}/{{sample}}/cellbender.err".format(LOGDIR),
            out="{}/{{sample}}/cellbender.out".format(LOGDIR)
        shell:
            """
            cellbender remove-background --input {input.h5} --output {output.cb} --fpr {params.fpr} --cuda --epochs {params.epochs} --expected-cells {params.expected_cells} {params.extra}
            mv ckpt.tar.gz {OUTDIR}/{wildcards.sample}/cellbender
            """
