## Need to integrate piscem_cpp with rust piscem
from os.path import join
from collections import OrderedDict
import yaml

k = config["k"]
m = 15

thr = config["thr"]

pisc_output_path = config["pisc_output_path"]

### index path
ind_output_path = join(pisc_output_path, "index")
pref_ind = config["prefix_index"]
pref_ind_k = f"{pref_ind}_k{k}"
ind_k_main_dir = join(ind_output_path, pref_ind_k)
ind_k_m_dir = join(ind_k_main_dir, f"m{m}")
ind_k_m_pref = join(ind_k_m_dir, f"{pref_ind_k}_m{m}")

### output path
map_output_path = join(pisc_output_path, "map_output")
out_data = join(map_output_path, "{data}")
out_dir_k_m = join(out_data, f"k_{k}", f"m_{m}")
out_dir_k_m_thr = join(out_dir_k_m, f"thr={{thr}}")
out_rad = join(out_dir_k_m_thr, "map.rad")
out_bed = join(out_dir_k_m_thr, "map.bed")
time_out = join(out_dir_k_m_thr, "time_align.out")

### data path
atac_data_path = config["raw_data_dir"]
data_names = config["data_names"]

### fastqs
Ls = [1,2]
data_name = "_".join("{data}".split("_")[:-1])

# fastqs = [[join(atac_data_path, "{data}", f"{data_name}_S3_L00{x}_{d_map[key]}_001.fastq.gz") for x in Ls] for key in d_map.keys()]

def get_fastq(data, key):
    d_map = OrderedDict({"read1":"R1", "read2":"R3", "barcode":"R2"})
    # data_name = "_".join(wildcards.data.split("_")[:-1])
    data_name = "_".join(data.split("_")[:-1])
    return [join(atac_data_path, data, f"{data_name}_S3_L00{x}_{d_map[key]}_001.fastq.gz") for x in Ls]

rule all_piscem_map:
    input:
        expand(time_out, data = data_names, thr = thr),
        expand(out_rad, data = data_names, thr = thr)

rule run_piscem_map:
    input:
        ind = f"{ind_k_m_pref}.sshash",
        read1 = lambda wildcards:get_fastq(wildcards.data, "read1"),
        read2 = lambda wildcards:get_fastq(wildcards.data, "read2"),
        barcode = lambda wildcards:get_fastq(wildcards.data, "barcode")
    output:
        out_rad = out_rad,
        time_out = time_out
    params:
        thr = lambda wildcards: wildcards.thr,
        ind_pref = ind_k_m_pref,
        pisc_atac = config["piscem_atac_path"],
        out_dir = out_dir_k_m_thr,
        threads = get_qos("run_piscem_map")["cpus_per_task"],
        read1 = lambda wildcards,input: ",".join(input.read1),
        read2 = lambda wildcards,input: ",".join(input.read2),
        barcode = lambda wildcards,input: ",".join(input.barcode)
    
    shell:
        """
            /usr/bin/time -o {output.time_out} {params.pisc_atac} \
                --index {params.ind_pref} \
                --read1 {params.read1} \
                --read2 {params.read2} \
                --barcode {params.barcode} \
                --output {params.out_dir} \
                --thr {params.thr} \
                --threads {params.threads}
        """