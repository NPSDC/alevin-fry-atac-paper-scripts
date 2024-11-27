chrom_params = config["chrom_params"]
chrom_k = list(chrom_params["k_w"].keys())
chrom_w = [chrom_params["k_w"][k] for k in chrom_k]


chromap_output_path = join(config["out_path"], "chromap")
chromap_ind_path = join(chromap_output_path, "index")
chrom_ind_ref_genome_path = join(chromap_ind_path, f"{{ref_genome}}_k={{k}}_w={{w}}")
chrom_time_ind_path = join(chrom_ind_ref_genome_path, f"time_ind_{{ref_genome}}.out")
chrom_ind_ref = join(chrom_ind_ref_genome_path, f"{{ref_genome}}.ind")

chromap_map_output_path = join(chromap_output_path, "map", f"{{ref_genome}}_k={{k}}_w={{w}}")
map_out_sam_file = join(chromap_map_output_path, f"map_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.sam")
map_out_sam_mf_file = join(chromap_map_output_path, f"map_mf_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.sam")

time_out_sam_file = join(chromap_map_output_path, f"time_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}_sam.out")

rule run_chromap_index:
    input:
        ref_file = lambda wc:refs[wc.ref_genome]
    output:
        ind_ref = chrom_ind_ref,
        time_out = chrom_time_ind_path

    params:
        threads = get_qos("run_chromap_index")["cpus_per_task"],
        chromap_soft = config["chromap_path"],
        k = lambda wc:wc['k'],
        w = lambda wc:wc['w']
    shell:
        """
            /usr/bin/time -o {output.time_out} {params.chromap_soft} -i -r {input.ref_file} -o {output.ind_ref} -t {params.threads} \
            -k {params.k} -w {params.w}
        """

rule run_chromap_map_sam:
    input:
        r1 = rules.run_mason.output.read1,
        r2 = rules.run_mason.output.read2,
        ref = lambda wc:refs[wc.ref_genome],
        index = rules.run_chromap_index.output.ind_ref
    params:
        threads = get_qos("run_chromap_map")["cpus_per_task"],
        chromap_soft = config["chromap_path"]
    output:
        out_sam = map_out_sam_file,
        time_out = time_out_sam_file
    shell:
        """
            /usr/bin/time -o {output.time_out} {params.chromap_soft} \
                -t {params.threads} \
                -x {input.index} \
                -r {input.ref} \
                -1 {input.r1} \
                -2 {input.r2} \
                -q 0 \
                --SAM \
                -o {output.out_sam}
        """

rule all_chromap:
    input:
        expand(expand(chrom_ind_ref, zip, k = chrom_k, w = chrom_w, allow_missing=True), ref_genome = r_gen),
        expand(expand(rules.run_chromap_map_sam.output, zip, k = chrom_k, w = chrom_w, allow_missing=True), 
            ref_genome = r_gen, nfrag = config["num_frags"], length = config["read_length"], 
            error_rate = config["error_rate"])