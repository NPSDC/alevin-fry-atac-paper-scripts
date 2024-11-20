from os.path import join
import yaml

# map_bed = join(out_data, "map.bed")
# map_sorted_bed_gz = f"{map_bed}.gz"
# macs2_pref = join(out_data, "macs2")
# macs2_out = f"{macs2_pref}_peaks.narrowPeak"
rdata_dir = join(out_chromap_data, "rdata") ## out_chromap in r
img_dir = join(out_chromap_data, "images")
out_rdata = join(rdata_dir, "seurat_ob_sub.rdata")

rule all_chromap_clusters:
    input:
        expand(out_rdata, data = data_names)

rule run_chromap_rscript:
    input:
        map_bed = rules.run_chromap_tabix.input, ## run_chromap_macs2.smk
        macs2_out = rules.run_chromap_macs2.output, ## run_chromap_macs2.smk
        tab_out = rules.run_chromap_tabix.output ## run_chromap_macs2.smk
    output:
        out_rdata
    params:
        data = lambda wildcards: wildcards.data,
        Rscript = join(config["R_path"], "Rscript"),
        out_dir = out_chromap_data,
        org = lambda wildcards: data_dict[wildcards.data]["org"]
    shell:
        """
            {params.Rscript} ../bash_scripts/proc_peaks.R \
                {input.map_bed} \
                {input.macs2_out} \
                {params.data} \
                chromap \
                {params.org} \
                {params.out_dir}
        """