from os.path import join
import yaml

map_sorted_bed = join(out_dir_k_m_thr, "map_sorted.bed")
map_sorted_bed_gz = f"{map_sorted_bed}.gz"
macs2_pref = join(out_dir_k_m_thr, "macs2")
macs2_out = f"{macs2_pref}_peaks.narrowPeak"
rdata_dir = join(out_dir_k_m_thr, "rdata")
img_dir = join(out_dir_k_m_thr, "images")
out_rdata = join(rdata_dir, "seurat_ob_sub.rdata")

rule all_piscem_clusters:
    input:
        expand(out_rdata, thr = thr, data = data_names)

rule run_piscem_rscript:
    input:
        map_bed = map_sorted_bed_gz,
        macs2_out = macs2_out
    output:
        out_rdata
    params:
        Rscript = join(config["R_path"], "Rscript"),
        data = lambda wildcards: wildcards.data,
        out_dir = out_dir_k_m_thr
    shell:
        """
            {params.Rscript} ../bash_scripts/proc_peaks.R \
                {input.map_bed} \
                {input.macs2_out} \
                {params.data} \
                piscem \
                hg38 \
                human \
                {params.out_dir}
        """