import sys
import os
import pickle

sys.path.insert(0, snakemake.params["pymod"])
import genotyp_output as go
import load_dump as ld

def save_stats(outstats):
    # samtools_stats=os.path.join(config['OutFolders']['main_results_dir'],config['samtools_stats'])
    samtools_stats=snakemake.output[3]
    ld.yaml_dumpDic(samtools_stats,outstats)
    return samtools_stats

def assembly_config_file(pk_refdatabase):
    # samtools_stats_yaml = os.path.join(config['OutFolders']['main_results_dir'],config['samtools_stats'])
    # readsconfig="/".join(config['readsconfig'].split("/")[1:])
    # readsconfig="/".join(snakemake.input[2].split("/")[1:])
    # asmconfig={'readsconfig':readsconfig,'samtools_stats':"/".join(samtools_stats_yaml.split("/")[1:]),'ref':"/".join(pk_refdatabase.split("/")[1:])}
    asmconfig={'readsconfig':snakemake.input[3],'samtools_stats':snakemake.output[3],'ref':snakemake.input[2]}
    # asm_configfile=os.path.join(config['OutFolders']['configs_fld'],config['assembly_config'])
    # ld.yaml_dumpDic(asm_configfile,asmconfig)
    ld.yaml_dumpDic(snakemake.output[4],asmconfig)


config=snakemake.config
config["MaxSheetNbr"] = ( config["outnbsheetperxls"] if config["outnbsheetperxls"]>=1 else 1 )

with open(snakemake.input[0], 'rb') as f:
    combined_stats=pickle.load(f)
with open(snakemake.input[1], 'rb') as f:
    reads_nbr=pickle.load(f)
pk_refdatabase=snakemake.input[2]

go.write_xls_output(
        config['outfolder'],
        config['OutFolders']['main_results_dir'],
        combined_stats,
        reads_nbr,
        pk_refdatabase,
        config["mismatchthrld"],
        # args.mismatchthrld,
        config)

go.write_genotypTXT(
        config['outfolder'],
        config['OutFolders']['main_results_dir'],
        combined_stats,
        pk_refdatabase)

samtools_stats=save_stats(combined_stats)
assembly_config_file(pk_refdatabase)

