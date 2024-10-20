def get_mason_fq(fname):
    # print(nfrag)
    if fname == "rsam":
        return join(mason_sim_path_ref, f"{fname}_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.sam")
    else:
        return join(mason_sim_path_ref, f"{fname}_nfrag={{nfrag}}_length={{length}}_er={{error_rate}}.fq")

mason_sim_path_ref = join(mason_sim_path, "{ref_genome}", f"read_length={{length}}")

rule run_mason:
    input:
        lambda wc:refs[wc.ref_genome]
    params:
        mason_sim = config["mason_sim_path"],
        paftools = f"{config["k8_path"]} {config["paftools_path"]}",
        seqtk = config["seqtk_path"],
        nfrags = lambda wildcards:wildcards.nfrag,
        read_length = lambda wildcards:wildcards.length,
        error_rate = lambda wildcards:wildcards.error_rate

    output:
        read1 = get_mason_fq("read1"),
        read2 = get_mason_fq("read2"),
        read3 = get_mason_fq("read3"),
        rsam = get_mason_fq("rsam"),
        read = get_mason_fq("read"),
        r1mf = get_mason_fq("r1mf"),
        r2mf = get_mason_fq("r2mf")

    shell:
        """
            ulimit -Sn 4096
            {params.mason_sim} -ir {input} \
                -n {params.nfrags} \
                -o {output.read1} \
                -or {output.read2} \
                -oa {output.rsam} \
                --illumina-read-length {params.read_length}\
                --illumina-prob-mismatch-scale {params.error_rate}
            
            {params.paftools} mason2fq {output.rsam} > {output.read}
            {params.seqtk} seq -1 {output.read} > {output.r1mf}
            {params.seqtk} seq -2 {output.read} > {output.r2mf}
            cp {output.read1} {output.read3}
        """


rule all_mason:
    input:
        expand(rules.run_mason.output.r2mf, ref_genome = r_gen, nfrag = config["num_frags"], length = config["read_length"], error_rate = config["error_rate"])
