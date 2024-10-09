eval_piscem_sam = join(piscem_map_out_path_dir, "eval.out")
eval_chromap_sam = join(chromap_map_output_path, f"eval_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.out")
eval_bowtie2_sam = join(bowtie2_map_output_path, f"eval_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.out")

# eval_paf_piscem_sam = join(piscem_map_out_path_dir, "eval_paf.out")
# eval_paf_chromap_sam = join(chromap_map_output_path, f"eval_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}_paf.out")
# eval_paf_bowtie2_sam = join(bowtie2_map_output_path, f"eval_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}_paf.out")

paftools = f"{config["k8_path"]} {config["paftools_path"]}"

rule run_eval_piscem:
    input:
        piscem = rules.run_piscem_map.output.out_piscem,
        # piscem_mf = rules.run_piscem_map_mf.output.out_piscem
    output:
        piscem = eval_piscem_sam,
        # paf = eval_paf_piscem_sam
    params:
        eval_path = join(config["eval_path"], "get_accuracy.py"),
        sam_truth = rules.run_mason.output.rsam,
        paftools = paftools
    shell:
        """
            python {params.eval_path} --truth {params.sam_truth} --predicted_sam {input.piscem} > {output.piscem}
            
        """
        #{params.paftools} mapeval {input.piscem_mf} > {output.paf} 
#{params.paftools} mapeval {input.piscem} > {output.paf} 
rule run_eval_chromap:
    input:
        chromap = rules.run_chromap_map_sam.output.out_sam,
        # chromap2 = rules.run_chromap_map_mf_sam.output.out_sam_mf
    output:
        chromap = eval_chromap_sam,
        # paf = eval_paf_chromap_sam
    params:
        eval_path = join(config["eval_path"], "get_accuracy.py"),
        sam_truth = rules.run_mason.output.rsam,
        paftools = paftools
    shell:
        """
            python {params.eval_path} --truth {params.sam_truth} --predicted_sam {input.chromap} > {output.chromap}
        """
     #{params.paftools} mapeval {input.chromap2} > {output.paf}

rule run_eval_bowtie2:
    input:
        bowtie2 = rules.run_bowtie2_map.output,
        # bowtie2_mf = rules.run_bowtie2_mf_map.output
    output:
        bowtie2 = eval_bowtie2_sam,
        # paf = eval_paf_bowtie2_sam
    params:
        eval_path = join(config["eval_path"], "get_accuracy.py"),
        sam_truth = rules.run_mason.output.rsam,
        paftools = paftools
    shell:
        """
            python {params.eval_path} --truth {params.sam_truth} --predicted_sam {input.bowtie2} > {output.bowtie2}
        """
        #        {params.paftools} mapeval {input.bowtie2_mf} > {output.paf}     

rule all_eval:
    input:
        expand(eval_piscem_sam, k = config["k"], m = config["m"], ref_genome = r_gen, nfrag = config["num_frags"],
            length = config["read_length"], error_rate = config["error_rate"], thr = config['thr'], 
            bin_size = config["bin_size"], orp = config["kmers_orphans"]),
        # expand(eval_paf_piscem_sam, k = config["k"], m = config["m"], ref_genome = r_gen, nfrag = config["num_frags"],
        #     length = config["read_length"], error_rate = config["error_rate"], thr = config['thr'],
        #     bin_size = config["bin_size"]),
        expand(expand(eval_chromap_sam, zip, k = chrom_k, w = chrom_w, allow_missing=True), ref_genome = r_gen, nfrag = config["num_frags"],
            length = config["read_length"], error_rate = config["error_rate"]),
        # expand(eval_paf_chromap_sam, ref_genome = r_gen, nfrag = config["num_frags"],
        #     length = config["read_length"], error_rate = config["error_rate"]),
        expand(eval_bowtie2_sam, ref_genome = r_gen, nfrag = config["num_frags"],
            length = config["read_length"], error_rate = config["error_rate"]),
        # expand(eval_paf_bowtie2_sam, ref_genome = r_gen, nfrag = config["num_frags"],
        #     length = config["read_length"], error_rate = config["error_rate"])