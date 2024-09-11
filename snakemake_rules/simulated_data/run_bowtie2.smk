bowtie2_output_path = join(config["out_path"], "bowtie2")
bowtie2_ind_path = join(bowtie2_output_path, "index")
bowtie2_ind_ref_genome_path = join(bowtie2_ind_path, "{ref_genome}")
bowtie2_time_ind_path = join(bowtie2_ind_ref_genome_path, f"time_ind_{{ref_genome}}.out")
bowtie2_ind_pref = join(bowtie2_ind_ref_genome_path, "{ref_genome}")

rule bowtie2_index:
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

rule all_bowtie2:
    input:
        expand(rules.bowtie2_index.output, ref_genome = r_gen)