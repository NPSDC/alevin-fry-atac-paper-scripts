eval_piscem_sam = join(piscem_map_out_path_dir, "eval.out")
eval_chromap_sam = join(chromap_map_output_path, f"eval_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.out")
eval_bowtie2_sam = join(bowtie2_map_output_path, f"eval_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.out")

rule run_eval_piscem:
    input:
        piscem = rules.run_piscem_map.output.out_piscem
    output:
        piscem = eval_piscem_sam
    params:
        eval_path = join(config["eval_path"], "get_accuracy.py"),
        sam_truth = rules.run_mason.output.rsam
    shell:
        """
            python {params.eval_path} --truth {params.sam_truth} --predicted_sam {input.piscem} > {output.piscem}
        """

rule run_eval_chromap:
    input:
        chromap = rules.run_chromap_map_sam.output.out_sam
    output:
        chromap = eval_chromap_sam
    params:
        eval_path = join(config["eval_path"], "get_accuracy.py"),
        sam_truth = rules.run_mason.output.rsam
    shell:
        """
            python {params.eval_path} --truth {params.sam_truth} --predicted_sam {input.chromap} > {output.chromap}
        """

rule run_eval_bowtie2:
    input:
        bowtie2 = rules.run_bowtie2_map.output
    output:
        bowtie2 = eval_bowtie2_sam
    params:
        eval_path = join(config["eval_path"], "get_accuracy.py"),
        sam_truth = rules.run_mason.output.rsam
    shell:
        """
            python {params.eval_path} --truth {params.sam_truth} --predicted_sam {input.bowtie2} > {output.bowtie2}
        """
rule all_eval:
    input:
        expand(eval_piscem_sam, k = config["k"], m = config["m"], ref_genome = r_gen, nfrag = config["num_frags"],
            length = config["read_length"], error_rate = config["error_rate"], thr = config['thr']),
        expand(eval_chromap_sam, ref_genome = r_gen, nfrag = config["num_frags"],
            length = config["read_length"], error_rate = config["error_rate"]),
        expand(eval_bowtie2_sam, ref_genome = r_gen, nfrag = config["num_frags"],
            length = config["read_length"], error_rate = config["error_rate"])