from os.path import join
import yaml

map_bed = join(out_data, "map.bed")
map_sorted_bed_gz = f"{map_bed}.gz"
macs2_pref = join(out_data, "macs2")
macs2_out = f"{macs2_pref}_peaks.narrowPeak"
rdata_dir = join(out_data, "rdata")
img_dir = join(out_data, "images")
out_rdata = join(rdata_dir, "seurat_ob_sub.rdata")

rule all_chromap_clusters:
    input:
        expand(out_rdata, data = data_names)

rule run_chromap_rscript:
    input:
        map_bed = map_sorted_bed_gz,
        macs2_out = macs2_out
    output:
        out_rdata
    params:
        Rscript = join(config["R_path"], "Rscript"),
        data = lambda wildcards: wildcards.data,
        out_dir = out_data
    shell:
        """
            {params.Rscript} ../bash_scripts/proc_peaks.R \
                {input.map_bed} \
                {input.macs2_out} \
                {params.data} \
                chromap \
                hg38 \
                human \
                {params.out_dir}
        """