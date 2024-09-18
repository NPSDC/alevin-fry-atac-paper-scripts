piscem_output_path = join(config["out_path"], "piscem")
piscem_ind_path = join(piscem_output_path, "index")
piscem_ind_ref_genome_path = join(piscem_ind_path, "{ref_genome}", f"k={{k}}_m={{m}}")
piscem_time_ind_path = join(piscem_ind_ref_genome_path, f"time_ind.out")
piscem_ind_pref = join(piscem_ind_ref_genome_path, "k={k}_m={m}")

piscem_map_out_path = join(piscem_output_path, "map", "{ref_genome}")
piscem_map_out_path_dir = join(piscem_map_out_path, f"map_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}_k{{k}}_m{{m}}_thr{{thr}}")
piscem_map_out_sam_file = join(piscem_map_out_path_dir, "map.sam")
time_out_piscem_map = join(piscem_map_out_path_dir, f"time_map.out")

rule run_pisc_index:
    input:
        lambda wc:refs[wc.ref_genome]
    output:
        [f"{piscem_ind_pref}{x}" for x in [".sshash", ".ectab", ".ctab", ".refinfo", "_cfish.json"]]
    params:
        ind_pref = piscem_ind_pref,
        time_ind = piscem_time_ind_path,
        threads = get_qos("run_chromap_index")["cpus_per_task"],
        k = lambda wildcards:wildcards.k,
        m = lambda wildcards:wildcards.m,
        piscem_exec_path = join(config["piscem_path"], "target", "release", "piscem"),
        tmpdir = join(config["tmp_dir"], "k{k}_m{m}")
    shell:
        """
            export TMPDIR=/scratch0
            ulimit -n 2048
            /usr/bin/time -o {params.time_ind} {params.piscem_exec_path} build \
                -s {input} \
                -k {params.k} \
                -m {params.m} \
                -t {params.threads} \
                -o {params.ind_pref} \
                -w {params.tmpdir} \
                --overwrite
        """

rule run_piscem_map:
    input:
        index = rules.run_pisc_index.output,
        r1 = rules.run_mason.output.read1,
        r2 = rules.run_mason.output.read2,
        r3 = rules.run_mason.output.read3   
    output:
        out_piscem = piscem_map_out_sam_file,
        out_time = time_out_piscem_map
    params:
        threads = get_qos("run_chromap_map")["cpus_per_task"],
        piscem_cpp = config["piscem_atac_path"],
        ind_pref = piscem_ind_pref,
        k = lambda wildcards:wildcards.k,
        m = lambda wildcards:wildcards.m,
        thr = lambda wildcards:wildcards.thr,
        out_dir = piscem_map_out_path_dir
    shell:
        """
            /usr/bin/time -o {output.out_time} {params.piscem_cpp} \
                --index {params.ind_pref} \
                --read1 {input.r1} \
                --read2 {input.r2} \
                --barcode {input.r3} \
                --output {params.out_dir} \
                --thr {params.thr} \
                --threads {params.threads} \
                --sam-format
        """

rule all_piscem:
    input:
        expand(rules.run_pisc_index.output, k = config["k"], m = config["m"],
            ref_genome = r_gen),
        expand(rules.run_piscem_map.output, k = config["k"], m = config["m"],
            thr = config['thr'], ref_genome = r_gen, nfrag = config["num_frags"],
            length = config["read_length"], error_rate = config["error_rate"])