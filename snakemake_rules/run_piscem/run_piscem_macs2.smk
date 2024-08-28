from os.path import join
import yaml

map_sorted_bed = join(out_dir_k_m_thr, "map_sorted.bed")
map_sorted_bed_gz = f"{map_sorted_bed}.gz"
map_tb = f"{map_sorted_bed}.gz.tbi"
macs2_pref = join(out_dir_k_m_thr, "macs2")
macs2_out = f"{macs2_pref}_peaks.narrowPeak"

rule all_piscem_macs2:
    input:
        expand(macs2_out, data = data_names, thr = thr),
        expand(map_tb, data = data_names, thr = thr)

rule run_piscem_macs2:
    input:
        map_sorted_bed
    output:
        macs2_out
    params:
        macs2_path = join(config["macs2_path"], "macs2"),
        macs2_pref = macs2_pref,
        q = config["q_macs2"],
        g = config["g_macs2"]
    shell:
        """
            {params.macs2_path} callpeak \
                -f BEDPE \
                -g {params.g} \
                --nomodel \
                --extsize 50 \
                --keep-dup all \
                -q {params.q} \
                -t {input} \
                -n {params.macs2_pref}
        """

rule run_piscem_tabix:
    input:
        map_sorted_bed_gz
    output:
        map_tb
    params:
        tabix = join(config["htslib_path"], "tabix")
    shell:
        """
            {params.tabix} -0 -p bed {input}
        """

rule run_piscem_bgzip:
    input:
        map_sorted_bed
    output:
        map_sorted_bed_gz
    params:
        bgzip = join(config["htslib_path"], "bgzip")
    shell:
        """
            {params.bgzip} -k {input}
        """

rule run_piscem_sort:
    input:
        out_bed
    output:
        map_sorted_bed
    shell:
        """
            sort -k1,1 -k2,2n -k3,3n {input} -o {output}
        """