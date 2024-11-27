from os.path import join
import yaml

map_sorted_bed_gz = f"{map_bed}.gz"
map_tb = f"{map_sorted_bed_gz}.tbi"
macs2_pref = join(out_dir_k_m_rem, "macs2")
macs2_out = f"{macs2_pref}_peaks.narrowPeak"

rule all_piscem_macs2:
    input:
        # expand(macs2_out, data = data_names, thr = [0.7], m = config["m"],
        #     k = [25], bin_size = [1000], orp = ["false"]),
        # expand(map_tb, data = data_names, thr = [0.7], m = config["m"], 
        #     k = [25], bin_size = [1000], orp = ["false"])
        expand(macs2_out, data = data_names, thr = [0.7], m = config["m"],
            k = [23], bin_size = [1000], orp = ["false"]),
        expand(map_tb, data = data_names, thr = [0.7], m = config["m"], 
            k = [23], bin_size = [1000], orp = ["false"])


rule run_piscem_macs2:
    input:
        rules.run_piscem_dedup.output
    output:
        macs2_out
    params:
        macs2_path = join(config["macs2_path"], "macs2"),
        macs2_pref = macs2_pref,
        q = config["q_macs2"],
        g = lambda wc:config["g_macs2"][data_dict[wc.data]["org"]],
        time_macs2 = join(out_dir_k_m_rem, "time_macs2.out")
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

rule run_piscem_bgzip:
    input:
        rules.run_piscem_dedup.output
    output:
        map_sorted_bed_gz
    params:
        bgzip = join(config["htslib_path"], "bgzip"),
        time_bgzip = join(out_dir_k_m_rem, "time_bgzip.out")
    shell:
        """
            /usr/bin/time -o {params.time_bgzip} {params.bgzip} -k {input}
        """

rule run_piscem_tabix:
    input:
        rules.run_piscem_bgzip.output
    output:
        map_tb
    params:
        tabix = join(config["htslib_path"], "tabix"),
        time_tabix = join(out_dir_k_m_rem, "time_tabix.out")
    shell:
        """
            /usr/bin/time -o {params.time_tabix} {params.tabix} -0 -p bed {input}
        """