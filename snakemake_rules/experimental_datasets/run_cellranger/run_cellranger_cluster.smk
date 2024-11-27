from os.path import join
import yaml

rdata_dir = join(out_cellranger_dir, "rdata") ## out_chromap in r
img_dir = join(out_cellranger_dir, "images")
out_rdata = join(rdata_dir, "seurat_ob_sub.rdata")

rule all_cellranger_clusters:
    input:
        expand(out_rdata, data = ["10k_pbmc_ATACv2_nextgem_Chromium_Controller_fastqs", 
        "8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_fastqs"])

rule run_cellranger_rscript:
    input:
        map_bed = frag_file, ## run_chromap_macs2.smk
        macs2_out = rules.run_cellranger_macs2.output, ## run_chromap_macs2.smk
        tab_out = map_tb ## run_chromap_macs2.smk
    output:
        out_rdata
    params:
        data = lambda wildcards: wildcards.data,
        Rscript = join(config["R_path"], "Rscript"),
        out_dir = out_cellranger_dir,
        org = lambda wildcards: data_dict[wildcards.data]["org"],
        time = join(out_cellranger_dir, "time_rscript.out") 
    shell:
        """
            /usr/bin/time -o {params.time} {params.Rscript} ../bash_scripts/proc_peaks.R \
                {input.map_bed} \
                {input.macs2_out} \
                {params.data} \
                cell_ranger \
                {params.org} \
                {params.out_dir}
        """