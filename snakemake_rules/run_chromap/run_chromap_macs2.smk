from os.path import join
import yaml

map_output_path = join(chromap_output_path, "map_output")
out_data = join(map_output_path, "{data}")
map_bed = join(out_data, "map.bed")
map_sorted_bed_gz = f"{out_file}.gz"
map_tb = f"{map_sorted_bed_gz}.tbi"
macs2_pref = join(out_data, "macs2")
macs2_out = f"{macs2_pref}_peaks.narrowPeak"

rule all_chromap_macs2:
    input:
        expand(macs2_out, data = data_names)

rule run_chromap_macs2:
    input:
        map_bed
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

rule run_chromap_tabix:
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

rule run_chromap_bgzip:
    input:
        map_bed
    output:
        map_sorted_bed_gz
    params:
        bgzip = join(config["htslib_path"], "bgzip")
    shell:
        """
            {params.bgzip} -k {input}
        """