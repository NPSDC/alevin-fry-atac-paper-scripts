from os.path import join
import yaml
if not workflow.overwrite_configfiles:
    configfile: "../config.yml"

# wf_path = config["wf_local"]
# wf_local = join(wf_path, "config.yaml")
# with open(wf_local, 'r') as f:
#     qos_type = yaml.safe_load(f)

# def get_qos(rule_name):
#     return qos_type[rule_name] if rule_name in qos_type else qos_type['default-resources']


### index path
ind_output_path = join(chromap_output_path, "index")
time_ind_path = join(ind_output_path, f"time_ind_{{org}}.out")
ind_ref = join(ind_output_path, f"{pref}.ind")
# input_ref_file = config["ref_fasta_input"]

rule all_chromap_index:
    input:
        expand(ind_ref, org=input_ref_org),
        expand(time_ind_path, org=input_ref_org)

rule run_chromap_index:
    input:
        ref_file = lambda wildcards:input_ref_dict[wildcards.org]
        
    output:
        ind_ref = ind_ref,
        time_out = time_ind_path

    params:
        threads = get_qos("run_chromap_index")["cpus_per_task"],
        chromap_soft = config["chromap_path"]
    shell:
        """
            /usr/bin/time -o {output.time_out} {params.chromap_soft} -i -r {input.ref_file} -o {output.ind_ref} -t {params.threads}
        """