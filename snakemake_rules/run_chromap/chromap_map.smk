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
out_file = join(out_data, "map.bed")
time_out = join(out_data, "time_align.out")

rule all_chromap_map:
    input:
        expand(time_out, data = data_names),
        expand(out_file, data = data_names),

rule run_chromap_map:
    input:
        ind = get_chrom_ind,
        read1 = lambda wildcards:get_fastq(wildcards.data, "read1"),
        read2 = lambda wildcards:get_fastq(wildcards.data, "read2"),
        barcode = lambda wildcards:get_fastq(wildcards.data, "barcode")
    output:
        time_out = time_out,
        out_file = out_file
    params:
        threads = get_qos("run_chromap_map")["cpus_per_task"],
        chromap_soft = config["chromap_path"],
        read1 = lambda wildcards,input: ",".join(input.read1),
        read2 = lambda wildcards,input: ",".join(input.read2),
        barcode = lambda wildcards,input: ",".join(input.barcode),
        whitelist = config["whitelist_file"],
        ref_file = lambda wc:input_ref_dict[config['data_dict'][wc.data]["org"]]
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
                --barcode-whitelist {params.whitelist} \
                --read-format bc:0:-1:-
        """