import sys
import os
import pickle

sys.path.insert(0, snakemake.params["pymod"])
import genotyp_output as go
import load_dump as ld

def save_stats(outstats):
    # samtools_stats=os.path.join(config['OutFolders']['main_results_dir'],config['samtools_stats'])
    samtools_stats=snakemake.output["sam"]
    ld.yaml_dumpDic(samtools_stats,outstats)
    return samtools_stats

def assembly_config_file(pk_refdatabase):
    # samtools_stats_yaml = os.path.join(config['OutFolders']['main_results_dir'],config['samtools_stats'])
    # readsconfig="/".join(config['readsconfig'].split("/")[1:])
    # readsconfig="/".join(snakemake.input[2].split("/")[1:])
    # asmconfig={'readsconfig':readsconfig,'samtools_stats':"/".join(samtools_stats_yaml.split("/")[1:]),'ref':"/".join(pk_refdatabase.split("/")[1:])}
    asmconfig={'readsconfig':snakemake.input[3],'samtools_stats':snakemake.output["sam"],'ref':snakemake.input[2]}
    # asm_configfile=os.path.join(config['OutFolders']['configs_fld'],config['assembly_config'])
    # ld.yaml_dumpDic(asm_configfile,asmconfig)
    ld.yaml_dumpDic(snakemake.output["asm"],asmconfig)

def SortKey(tup):
    key, d = tup
    return d[config['SortKey']]

def write_tsv_output(outfile,results_folder,outstats,reads_nbr,pk_refdatabase,mismatchthrld,inconfig):

    global config
    config = inconfig
    
    refinfo=ld.pickle_loadDic(pk_refdatabase)
    
    outf = open(outfile, "w")
    sample_list=sorted(outstats.keys())
    samples_names=sample_list
    headers=["Sample", "Allele", "Reads mapped", 'Genotyp Score','Error rate','Mean covered Depth','Normalized covered Depth','Homolog IDs','Reference Length (bp)','bases mapped','mismatches','Region Coverage','reads <{} mismatchs ratio'.format(mismatchthrld),'Mean covered Depth mismatch ratio corrected','Normalized covered Depth mismatch ratio corrected','Warnings']
    outf.writelines("\t".join(headers) + "\n")
    
    for sample in samples_names:
        items = sorted(outstats[sample].items(), key = SortKey, reverse=config['sortOrderReverse'])
        for allele,stats in items:
            reads_ratio=float(stats['reduced_stats']['reads mapped'])/float(stats['stats']['reads mapped'])
            values=[sample,allele, stats["stats"]["reads mapped"], stats['gscore'],stats['stats']['error rate'],stats['coverage']['depth_covered_mean'],stats['coverage']['NORM_depth_covered_mean'],refinfo[allele]['grpRef'],refinfo[allele]['len_seq'],stats['stats']['bases mapped'],stats['stats']['mismatches'],stats['coverage']['covered_pct'],reads_ratio,stats['coverage']['depth_covered_mean']*reads_ratio,stats['coverage']['NORM_depth_covered_mean']*reads_ratio,""]
            outf.writelines("\t".join([str(i) for i in values]) + "\n")
    outf.close()

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

write_tsv_output(
        snakemake.output["tsv"],
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

