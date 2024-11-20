from os.path import join
import yaml
if not workflow.overwrite_configfiles:
    configfile: "../config.yml"

def get_chrom_ind(wildcards):
    ref_gen = config['data_dict'][wildcards.data]["org"]
    return rules.run_chromap_index.output.ind_ref.format(org=ref_gen)

### map path
map_output_path = join(chromap_output_path, "map_output")

#### Whitelist
out_wh_data = join(map_output_path, "{data}", "map_wh")
chromap_wh_log1 = join(out_wh_data, "log_threads={threads}_1.out")
chromap_wh_log2 = join(out_wh_data, "log_threads={threads}_2.out")
chromap_wh_log3 = join(out_wh_data, "log_threads={threads}_3.out")
out_wh_file = join(out_wh_data, "map.bed")
time_wh_out_1 = join(out_wh_data, "time_align_threads={threads}_1.out")
time_wh_out_2 = join(out_wh_data, "time_align_threads={threads}_2.out")
time_wh_out_3 = join(out_wh_data, "time_align_threads={threads}_3.out")

#### map only
out_monly_data = join(map_output_path, "{data}", "map_only")
chromap_monly_log1 = join(out_monly_data, "log_threads={threads}_1.out")
chromap_monly_log2 = join(out_monly_data, "log_threads={threads}_2.out")
chromap_monly_log3 = join(out_monly_data, "log_threads={threads}_3.out")
out_monly_file = join(out_monly_data, "map.bed")
time_monly_out_1 = join(out_monly_data, "time_align_threads={threads}_1.out")
time_monly_out_2 = join(out_monly_data, "time_align_threads={threads}_2.out")
time_monly_out_3 = join(out_monly_data, "time_align_threads={threads}_3.out")

rule all_chromap_map:
    input:
        expand(time_wh_out_1, data = ["10k_pbmc_ATACv2_nextgem_Chromium_Controller_fastqs", "8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_fastqs"], threads = threads),
        expand(time_monly_out_1, data = ["10k_pbmc_ATACv2_nextgem_Chromium_Controller_fastqs", "8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_fastqs"], threads = threads),
#        expand(out_wh_file, data = data_names),
#        expand(out_monly_file, data = data_names)

rule run_chromap_wh:
    input:
        ind = get_chrom_ind,
        read1 = lambda wildcards:get_fastq(wildcards.data, "read1"),
        read2 = lambda wildcards:get_fastq(wildcards.data, "read2"),
        barcode = lambda wildcards:get_fastq(wildcards.data, "barcode")
    output:
        time_out_1 = time_wh_out_1,
        time_out_2 = time_wh_out_2,
        time_out_3 = time_wh_out_3,
    params:
        out_log1 = chromap_wh_log1,
        out_log2 = chromap_wh_log2,
        out_log3 = chromap_wh_log3,
        out_file = out_wh_file,
        threads = lambda wc:wc.threads,
        chromap_soft = config["chromap_path"],
        read1 = lambda wildcards,input: ",".join(input.read1),
        read2 = lambda wildcards,input: ",".join(input.read2),
        barcode = lambda wildcards,input: ",".join(input.barcode),
        whitelist_file = lambda wildcards: whl_map[data_dict[wildcards.data]['whl_type']],
        rc = lambda wildcards: "--read-format bc:0:-1:-" if data_dict[wildcards.data]['rc'] else "",
        ref_file = lambda wc:input_ref_dict[data_dict[wc.data]["org"]]
    shell:
        """
            /usr/bin/time -o {output.time_out_1} {params.chromap_soft} \
                --preset atac \
                -t {params.threads} \
                -x {input.ind} \
                -r {params.ref_file} \
                -1 {params.read1} \
                -2 {params.read2} \
                -b {params.barcode} \
                -o {params.out_file} \
                --barcode-whitelist {params.whitelist_file} \
                {params.rc} >{params.out_log1} 2>&1
	   rm {params.out_file}
           /usr/bin/time -o {output.time_out_2} {params.chromap_soft} \
                --preset atac \
                -t {params.threads} \
                -x {input.ind} \
                -r {params.ref_file} \
                -1 {params.read1} \
                -2 {params.read2} \
                -b {params.barcode} \
                -o {params.out_file} \
                --barcode-whitelist {params.whitelist_file} \
                {params.rc} >{params.out_log2} 2>&1
	   rm {params.out_file}
           /usr/bin/time -o {output.time_out_3} {params.chromap_soft} \
                --preset atac \
                -t {params.threads} \
                -x {input.ind} \
                -r {params.ref_file} \
                -1 {params.read1} \
                -2 {params.read2} \
                -b {params.barcode} \
                -o {params.out_file} \
                --barcode-whitelist {params.whitelist_file} \
                {params.rc} >{params.out_log3} 2>&1
        """

rule run_chromap_map:
    input:
        ind = get_chrom_ind,
        read1 = lambda wildcards:get_fastq(wildcards.data, "read1"),
        read2 = lambda wildcards:get_fastq(wildcards.data, "read2"),
        barcode = lambda wildcards:get_fastq(wildcards.data, "barcode")
    output:
        time_out_1 = time_monly_out_1,
        time_out_2 = time_monly_out_2,
        time_out_3 = time_monly_out_3,
    params:
        threads = lambda wc:wc.threads,
        chromap_soft = config["chromap_path"],
        read1 = lambda wildcards,input: ",".join(input.read1),
        read2 = lambda wildcards,input: ",".join(input.read2),
        barcode = lambda wildcards,input: ",".join(input.barcode),
        whitelist_file = lambda wildcards: whl_map[data_dict[wildcards.data]['whl_type']],
        rc = lambda wildcards: "--read-format bc:0:-1:-" if data_dict[wildcards.data]['rc'] else "",
        ref_file = lambda wc:input_ref_dict[data_dict[wc.data]["org"]],
        out_file = out_monly_file,
        out_log1 = chromap_monly_log1,
        out_log2 = chromap_monly_log2,
        out_log3 = chromap_monly_log3
    shell:
        """
            /usr/bin/time -o {output.time_out_1} {params.chromap_soft} \
                -t {params.threads} \
                -x {input.ind} \
                -r {params.ref_file} \
                -1 {params.read1} \
                -2 {params.read2} \
                -b {params.barcode} \
                -o {params.out_file} \
                {params.rc} >{params.out_log1} 2>&1
	   rm {params.out_file}
           /usr/bin/time -o {output.time_out_2} {params.chromap_soft} \
                -t {params.threads} \
                -x {input.ind} \
                -r {params.ref_file} \
                -1 {params.read1} \
                -2 {params.read2} \
                -b {params.barcode} \
                -o {params.out_file} \
                {params.rc} >{params.out_log2} 2>&1
	   rm {params.out_file}
           /usr/bin/time -o {output.time_out_3} {params.chromap_soft} \
                -t {params.threads} \
                -x {input.ind} \
                -r {params.ref_file} \
                -1 {params.read1} \
                -2 {params.read2} \
                -b {params.barcode} \
                -o {params.out_file} \
                {params.rc} >{params.out_log3} 2>&1
        """
