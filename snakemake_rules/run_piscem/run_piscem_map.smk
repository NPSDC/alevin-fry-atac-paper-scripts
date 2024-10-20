### output path
map_output_path = join(pisc_output_path, "map_output")
out_data = join(map_output_path, "{data}")
out_dir_k_m = join(out_data, f"k{{k}}_m_{{m}}")
out_dir_k_m_rem = join(out_dir_k_m, f"bin-size={{bin_size}}_thr={{thr}}_orp={{orp}}")
out_rad = join(out_dir_k_m_rem, "map.rad")
out_bed = join(out_dir_k_m_rem, "map.bed")
time_out = join(out_dir_k_m_rem, "time_align.out")

# fastqs = [[join(atac_data_path, "{data}", f"{data_name}_S3_L00{x}_{d_map[key]}_001.fastq.gz") for x in Ls] for key in d_map.keys()]

rule all_piscem_map:
    input:
        expand(time_out, data = data_names, m = config["m"], k = config["k"], thr = config["thr"], 
            orp = config["kmers_orphans"], bin_size = config["bin_size"]),
        # expand(out_rad, data = data_names, m = config["m"], k = config["k"], thr = config["thr"], 
        #     orp = config["kmers_orphans"], bin_size = config["bin_size"]),
        # expand(out_bed, data = data_names, m = config["m"], k = config["k"], thr = config["thr"], 
        #     orp = config["kmers_orphans"], bin_size = config["bin_size"])
        expand(out_bed, data = data_names, m = config["m"], k = [23, 25], thr = [0.7, 1], 
            orp = ["false"], bin_size = 1000)

rule run_piscem_map:
    input:
        ind = lambda wildcards:f"{ind_k_m_pref}.sshash".format(m=wildcards.m,
                k=wildcards.k, org = data_dict[wildcards.data]["org"]),
        read1 = lambda wildcards:get_fastq(wildcards.data, "read1"),
        read2 = lambda wildcards:get_fastq(wildcards.data, "read2"),
        barcode = lambda wildcards:get_fastq(wildcards.data, "barcode")
    output:
        # out_rad = out_rad,
        time_out = time_out
    params:
        ind_pref = lambda wildcards:ind_k_m_pref.format(m = wildcards.m,
                        k = wildcards.k, org = data_dict[wildcards.data]["org"]),
        piscem_exec_path = join(piscem_exec_path, "target", "release", "piscem"),
        out_dir = out_dir_k_m_rem,
        out_rad = out_rad,
        threads = get_qos("run_piscem_map")["cpus_per_task"],
        read1 = lambda wildcards,input: ",".join(input.read1),
        read2 = lambda wildcards,input: ",".join(input.read2),
        barcode = lambda wildcards,input: ",".join(input.barcode),
        k = lambda wildcards:wildcards.k,
        m = lambda wildcards:wildcards.m,
        thr = lambda wildcards:wildcards.thr,
        bin_size = lambda wildcards:wildcards.bin_size,
        orp = lambda wildcards: "--check-kmer-orphan" if wildcards.orp=="true" else " ",
        rm = lambda wildcards: "yes" if int(wildcards.k) > 25 or int(wildcards.bin_size) > 1000 else ""
    
    shell:
        """
            /usr/bin/time -o {output.time_out} {params.piscem_exec_path} \
                map-sc-atac \
                --index {params.ind_pref} \
                --read1 {params.read1} \
                --read2 {params.read2} \
                --barcode {params.barcode} \
                --output {params.out_dir} \
                --thr {params.thr} \
                --bin-size {params.bin_size} \
                {params.orp} \
                --threads {params.threads}
            echo "{params.rm}"
            if [ "{params.rm}" == "yes" ]; then
                rm {params.out_rad}
            fi
        """

rule run_piscem_dedup:
    input:
        out_rad
    output:
        out_bed
    params:
        map_dir = out_dir_k_m_rem,
        threads = get_qos("run_piscem_map")["cpus_per_task"],
        pisc_dpath = config["piscem_dedup_path"],
        whitelist_file = lambda wc: whl_map[data_dict[wc.data]['whl_type']],
        rev_comp = lambda wc: data_dict[wc.data]['rc']
    shell:
        """
            ../bash_scripts/run_piscem_dedup.sh {params.pisc_dpath} \
                {params.map_dir} {params.whitelist_file} \
                {params.rev_comp} {params.threads} {params.map_dir}
        """    