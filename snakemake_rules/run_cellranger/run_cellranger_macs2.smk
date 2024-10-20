from os.path import join
import yaml

out_cellranger_dir = join(config["cellranger_output_path"], "{data}")
frag_file = join(out_cellranger_dir, "fragments.tsv.gz")
map_tb = join(out_cellranger_dir, "fragments.tsv.gz")
macs2_pref = join(out_cellranger_dir, "macs2")
macs2_out = f"{macs2_pref}_peaks.narrowPeak"

rule all_piscem_macs2:
    input:
        expand(macs2_out, data = data_names)

rule run_cellrange_map:
    input:
        frag_file
    output:
        macs2_out
    params:
        macs2_path = join(config["macs2_path"], "macs2"),
        macs2_pref = macs2_pref,
        q = config["q_macs2"],
        g = lambda wc:config["g_macs2"][data_dict[wc.data]["org"]],
        time_macs2 = join(out_cellranger_dir, "time_macs2.out")
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