bowtie2_output_path = join(config["out_path"], "bowtie2")
bowtie2_ind_path = join(bowtie2_output_path, "index")
bowtie2_ind_ref_genome_path = join(bowtie2_ind_path, "{ref_genome}")
bowtie2_time_ind_path = join(bowtie2_ind_ref_genome_path, f"time_ind_{{ref_genome}}.out")
bowtie2_ind_pref = join(bowtie2_ind_ref_genome_path, "{ref_genome}")

bowtie2_map_output_path = join(bowtie2_output_path, "map", "{ref_genome}")
map_out_bowtie2_sam_file = join(bowtie2_map_output_path, f"map_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.sam")
map_out_bowtie2_sam_mf_file = join(bowtie2_map_output_path, f"map_mf_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.sam")
time_out_bowtie2_map = join(bowtie2_map_output_path, f"time_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.out")

rule run_bowtie2_index:
    input:
        lambda wc:refs[wc.ref_genome]
    output:
        [f"{bowtie2_ind_pref}.{x}.bt2" for x in range(1,4)],
        [f"{bowtie2_ind_pref}.rev.{x}.bt2" for x in range(1,2)]
    params:
        threads = 32,
        time_ind = bowtie2_time_ind_path,
        pref = bowtie2_ind_pref,
        bowtie2_build = join(config["bowtie2_path"], "bowtie2-build")
    shell:
        """
            /usr/bin/time -o {params.time_ind} {params.bowtie2_build} \
                --threads {params.threads} {input} {params.pref}
        """

rule run_bowtie2_map:
    input:
        rules.run_bowtie2_index.output
    output:
        map_out_bowtie2_sam_file
    params:
        threads = 32,
        time_out = time_out_bowtie2_map,
        bowtie2 = join(config["bowtie2_path"], "bowtie2"),
        pref = rules.run_bowtie2_index.params.pref,
        r1 = rules.run_mason.output.read1,
        r2 = rules.run_mason.output.read2
    shell:
        """
            {params.bowtie2} -X2000 --threads {params.threads} \
            -x {params.pref} -1 {params.r1} \
            -2 {params.r2} > {output}
        """

rule run_bowtie2_mf_map:
    input:
        rules.run_bowtie2_index.output
    output:
        map_out_bowtie2_sam_mf_file
    params:
        threads = 32,
        bowtie2 = join(config["bowtie2_path"], "bowtie2"),
        pref = rules.run_bowtie2_index.params.pref,
        r1 = rules.run_mason.output.r1mf,
        r2 = rules.run_mason.output.r2mf
    shell:
        """
            {params.bowtie2} -X2000 --threads {params.threads} \
            -x {params.pref} -1 {params.r1} \
            -2 {params.r2} > {output}
        """

rule all_bowtie2:
    input:
        expand(rules.run_bowtie2_index.output, ref_genome = r_gen),
        expand(rules.run_bowtie2_map.output, ref_genome = r_gen, 
            nfrag = config["num_frags"], length = config["read_length"],
            error_rate = config["error_rate"]),
        expand(rules.run_bowtie2_mf_map.output, ref_genome = r_gen, 
            nfrag = config["num_frags"], length = config["read_length"],
            error_rate = config["error_rate"])