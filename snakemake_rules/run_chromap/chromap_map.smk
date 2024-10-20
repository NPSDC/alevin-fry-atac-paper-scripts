from os.path import join
import yaml
if not workflow.overwrite_configfiles:
    configfile: "../config.yml"

def get_chrom_ind(wildcards):
    ref_gen = config['data_dict'][wildcards.data]["org"]
    return rules.run_chromap_index.output.ind_ref.format(org=ref_gen)

### map path
map_output_path = join(chromap_output_path, "map_output")
out_data = join(map_output_path, "{data}")
out_wc_file = join(out_data, "map_wc.bed")
out_wc_log = join(out_data, "map_wc.log")
time_wc_out = join(out_data, "time_wc_align.out")

out_only_file = join(out_data, "map_only.bed")
out_only_log = join(out_data, "map_only_wc.log")
time_only_out = join(out_data, "time_only_align.out")


rule all_chromap_map:
    input:
        expand(time_only_out, data = data_names),
        expand(out_only_file, data = data_names),
        expand(out_wc_file, data = data_names)

rule run_chromap_map_wc:
    input:
        ind = get_chrom_ind,
        read1 = lambda wildcards:get_fastq(wildcards.data, "read1"),
        read2 = lambda wildcards:get_fastq(wildcards.data, "read2"),
        barcode = lambda wildcards:get_fastq(wildcards.data, "barcode")
    output:
        time_out = time_wc_out,
        out_file = out_wc_file
    params:
        threads = get_qos("run_chromap_map")["cpus_per_task"],
        chromap_soft = config["chromap_path"],
        read1 = lambda wildcards,input: ",".join(input.read1),
        read2 = lambda wildcards,input: ",".join(input.read2),
        barcode = lambda wildcards,input: ",".join(input.barcode),
        whitelist_file = lambda wildcards: whl_map[data_dict[wildcards.data]['whl_type']],
        rc = lambda wildcards: "--read-format bc:0:-1:-" if data_dict[wildcards.data]['rc'] else "",
        ref_file = lambda wc:input_ref_dict[data_dict[wc.data]["org"]],
        out_log = out_wc_log
    shell:
        """
            /usr/bin/time -o {output.time_out} {params.chromap_soft} \
                --preset atac \
                -t {params.threads} \
                -x {input.ind} \
                -r {params.ref_file} \
                -1 {params.read1} \
                -2 {params.read2} \
                -b {params.barcode} \
                -o {output.out_file} \
                --barcode-whitelist {params.whitelist_file} \
                {params.rc} >{params.out_log} 2>&1 
        """

rule run_chromap_only_map:
    input:
        ind = get_chrom_ind,
        read1 = lambda wildcards:get_fastq(wildcards.data, "read1"),
        read2 = lambda wildcards:get_fastq(wildcards.data, "read2"),
        barcode = lambda wildcards:get_fastq(wildcards.data, "barcode")
    output:
        time_out = time_only_out,
        out_file = out_only_file
    params:
        threads = get_qos("run_chromap_map")["cpus_per_task"],
        chromap_soft = config["chromap_path"],
        read1 = lambda wildcards,input: ",".join(input.read1),
        read2 = lambda wildcards,input: ",".join(input.read2),
        barcode = lambda wildcards,input: ",".join(input.barcode),
        rc = lambda wildcards: "--read-format bc:0:-1:-" if data_dict[wildcards.data]['rc'] else "",
        ref_file = lambda wc:input_ref_dict[data_dict[wc.data]["org"]],
        out_log = out_only_log
    shell:
        """
            /usr/bin/time -o {output.time_out} {params.chromap_soft} \
                -t {params.threads} \
                -x {input.ind} \
                -r {params.ref_file} \
                -1 {params.read1} \
                -2 {params.read2} \
                -b {params.barcode} \
                -o {output.out_file} \
                {params.rc} >{params.out_log} 2>&1
        """