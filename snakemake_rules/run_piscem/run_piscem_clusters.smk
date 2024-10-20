from os.path import join
import yaml

# map_sorted_bed = join(out_dir_k_m_thr, "map_sorted.bed")
# map_sorted_bed_gz = f"{map_sorted_bed}.gz"
# macs2_pref = join(out_dir_k_m_thr, "macs2")
# macs2_out = f"{macs2_pref}_peaks.narrowPeak"
rdata_dir = join(out_dir_k_m_rem, "rdata")
img_dir = join(out_dir_k_m_rem, "images")
out_rdata = join(rdata_dir, "seurat_ob_sub.rdata")

rule all_piscem_clusters:
    input:
        expand(out_rdata, thr = [0.7], data = data_names, m = config["m"],
         k = [25], bin_size = [1000], orp = ["false"])

rule run_piscem_rscript:
    input:
        map_bed = rules.run_piscem_tabix.input, ## run_piscem_macs2.smk
        macs2_out = rules.run_piscem_macs2.output, ## run_piscem_macs2.smk
        tab_out = rules.run_piscem_tabix.output ## run_piscem_macs2.smk
    output:
        out_rdata
    params:
        data = lambda wildcards: wildcards.data,
        Rscript = join(config["R_path"], "Rscript"),
        org = lambda wildcards: data_dict[wildcards.data]["org"],
        out_dir = out_dir_k_m_rem
    shell:
        """
            {params.Rscript} ../bash_scripts/proc_peaks.R \
                {input.map_bed} \
                {input.macs2_out} \
                {params.data} \
                piscem \
                {params.org} \
                {params.out_dir}
        """