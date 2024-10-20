from os.path import join
import yaml

map_output_path = join(chromap_output_path, "map_output")
out_chromap_data = join(map_output_path, "{data}")
map_bed = join(out_chromap_data, "map_wc.bed")
map_sorted_bed_gz = f"{map_bed}.gz"
map_tb = f"{map_sorted_bed_gz}.tbi"
macs2_pref = join(out_chromap_data, "macs2")
macs2_out = f"{macs2_pref}_peaks.narrowPeak"
time_macs2 = join(out_chromap_data, "time_macs2.out")
time_bgzip = join(out_chromap_data, "time_bgzip.out")
time_tabix = join(out_chromap_data, "time_tabix.out")

rule all_chromap_macs2:
    input:
        expand(macs2_out, data = data_names),
        expand(map_tb, data = data_names)

rule run_chromap_macs2:
    input:
        rules.run_chromap_map_wc.output.out_file
    output:
        macs2_out
    params:
        macs2_path = join(config["macs2_path"], "macs2"),
        macs2_pref = macs2_pref,
        time_macs2 = time_macs2,
        q = config["q_macs2"],
        g = lambda wc:config["g_macs2"][data_dict[wc.data]["org"]]
    shell:
        """
            /usr/bin/time -o {params.time_macs2} {params.macs2_path} callpeak \
                -f BEDPE \
                -g {params.g} \
                --nomodel \
                --extsize 50 \
                --keep-dup all \
                -q {params.q} \
                -t {input} \
                -n {params.macs2_pref}
        """

rule run_chromap_bgzip:
    input:
        rules.run_chromap_map_wc.output.out_file
    output:
        map_sorted_bed_gz
    params:
        bgzip = join(config["htslib_path"], "bgzip"),
        time_bgzip = time_bgzip
    shell:
        """
            /usr/bin/time -o {params.time_bgzip} {params.bgzip} -k {input}
        """

rule run_chromap_tabix:
    input:
        rules.run_chromap_bgzip.output
    output:
        map_tb
    params:
        tabix = join(config["htslib_path"], "tabix"),
        time_tabix = time_tabix
    shell:
        """
            /usr/bin/time -o {params.time_tabix} {params.tabix} -0 -p bed {input}
        """