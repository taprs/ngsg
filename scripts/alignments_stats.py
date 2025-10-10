#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np

sys.path.insert(0, snakemake.params["pymod"])
import calculate_genotypes as cg
import genotyp_depth as gd
import load_dump as ld
import bamanalyse

def alignments_stats(bwt_results,pk_refdatabase):
    # stdout_print("Do Alignments stats...",printInLog=True)
    # tps_stat=time.time()
    outstats={}
    refinfo=ld.pickle_loadDic(pk_refdatabase)
    depthfiles={}
    reads_nbr={}
    # sample_nbr=len(bwt_results.keys())
    # c=0
    # tps1=0
    # ellapsed_times=[]
    for sample, refdatas in bwt_results.items():
        # sample=snakemake.wildcards["sample"]
        # refID=snakemake.wildcards["ref"]
        outstats[sample]={}
        depthfiles[sample]={}
        # c+=1
        # if tps1>0:
            # ellapsed_times.append(time.time()-tps1)
            # ETE=float(np.mean(ellapsed_times))*(sample_nbr-(c-1))
        # else:
            # ETE=0
        # stdout_print("\t{}/{} => do stats for {} -- ETE {} s".format(c,sample_nbr,sample,round(ETE,2)),printInLog=True)
        # tps1 = time.time()
        for refID, idval in refdatas.items():
            stats_file=idval['stats']
            reduced_stats_file=idval['reduced_stats']
    # stats_file, reduced_stats_file = snakemake.input[0:2]
    # idval["bam_sorted"], idval["reduced_sorted"] = snakemake.input[2:4]
            reflength=refinfo[refID]['len_seq']
            stats_results=parse_stats_file(stats_file)
            reduced_stats_results=parse_stats_file(reduced_stats_file)
            depthyaml=os.path.join(idval['aln_path'],'depthinfos_{}'.format(idval['aln_name']))
            reduced_depthyaml=os.path.join(idval['aln_path'],'reduced_depthinfos_{}'.format(idval['aln_name']))
            depthfile=os.path.join(idval['aln_path'],'{}_SORTED.bam.depth'.format(idval['aln_name']))
            reduced_depthfile=os.path.join(idval['aln_path'],'{}_reduced.bam.depth'.format(idval['aln_name']))
            if stats_results['reads mapped']>snakemake.config['readsMappedThrld']:
                print(idval['aln_name'])
                depthfiles[sample][refID]={'sorted':depthyaml,'reduced':reduced_depthyaml}
                outstats[sample][refID]={'stats':stats_results}
                outstats[sample][refID].update({'reduced_stats':reduced_stats_results})
                outstats[sample][refID].update(analyse_bam(idval['bam_sorted']))
                outstats[sample][refID].update({'ref_length':reflength})
                outstats[sample][refID].update({'coverage':extract_coverage(depthfile,depthyaml,reflength)})
                outstats[sample][refID].update({'reduced_coverage':extract_coverage(reduced_depthfile,reduced_depthyaml,reflength)})
                # outstats[sample][refID].update({'coverage':extract_coverage(idval['bam_sorted'],depthyaml,reflength)})
                # outstats[sample][refID].update({'reduced_coverage':extract_coverage(idval['reduced_sorted'],reduced_depthyaml,reflength)})
            if sample not in reads_nbr.keys():
                reads_nbr[sample]=get_reads_nbr(bwt_results[sample][refID]['errBWTfile'])
    # stdout_print("END Alignments stats in {} s".format(round(time.time()-tps_stat,2)),printInLog=True)
    # return outstats,depthfiles
    return outstats,depthfiles,reads_nbr

def parse_stats_file(stats_file):
    if os.path.exists(stats_file):
        filestats={}
        with open(stats_file) as stfile:
            for line in stfile:
                tmp=[v.strip() for v in line.split(':')]
                keyid=tmp[0]
                value=float(tmp[1])
                filestats[keyid]=value
            if filestats['bases mapped']>0:
                error_rate=filestats['mismatches']/filestats['bases mapped']
                return {'raw sequences':filestats['raw total sequences'],'reads mapped':filestats['reads mapped'],'bases mapped':filestats['bases mapped'],'mismatches':filestats['mismatches'],'average length':filestats['average length'],'error rate':error_rate}
    return {'reads mapped':0,'bases mapped':0,'mismatches':0,'average length':0,'error rate':0}

def analyse_bam(bamfile):
    #{'ratio_readsNoInDel':ratio_readsNoInDel,'readsNoInDelNbr':readsNoInDelNbr,'mismatchs_list':XM_list,'mismatch_mean':float(np.mean(XM_list))}

    # mismatchthrld=args.mismatchthrld
    mismatchthrld=snakemake.config["mismatchthrld"]
    
    bamdatas=bamanalyse.read_bam(bamfile)
    XM_list=bamdatas['mismatchs_list']
    
    if len(XM_list)>0:
        n_bins=max(XM_list)-min(XM_list)
        if n_bins==0: n_bins=1
        #if 0 in XM_list: print(bamfile)
        distrib_XM_values=[[int(i) for i in v] for v in np.histogram(XM_list,bins=n_bins)]
    else:
        distrib_XM_values=[]
    Nb_below_Thrld=sum([1 for v in XM_list if v<=mismatchthrld])
    ratio_below_Thrld = 0
    if len(XM_list)>0:
        ratio_below_Thrld = float(Nb_below_Thrld)/float(len(XM_list))
        
    #MAPQ: MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. A value 255 indicates that the mapping quality is not available.
    MQ_list=[v for v in bamdatas['MAPQ_list'] if v!=255]
    if len(MQ_list)>0:
        mean_MQ = float(np.mean(MQ_list))
    else:
        mean_MQ = 0
    #prob all mapping positions are wrong = (1/10)^(PHRED_SCORE/10)
    MQ_prob=0.1**(mean_MQ/10)
    
    return {'readsNb':bamdatas['readsNb'],'ratio_readsNoInDel':bamdatas['ratio_readsNoInDel'],'readsNoInDelNbr':bamdatas['readsNoInDelNbr'],'mismatch_mean':bamdatas['mismatch_mean'],'distrib_XM_values':distrib_XM_values,'MQ_prob':MQ_prob}

# def extract_coverage(bamfile,depthyaml,reflength):
def extract_coverage(depthfile,depthyaml,reflength):
    # cmd = "{path}samtools depth {bam} | cut -f 2-".format(path=samtools_path,bam=bamfile)
    # depthrstl = utils.run(args,"samtools depth", '', cmd, retOUT=True)
    # depthstdout=depthrstl[0].decode('utf-8')
    
    with open(depthfile, 'r') as file:
        depthstdout=file.readlines()

    depth_datas={}
    # for depthline in depthstdout.split('\n'):
    for depthline in depthstdout:
        depth_array=depthline.split('\t')
        if len(depth_array)==2:
            depth_datas[int(depth_array[0])]=int(depth_array[1])
            
    depth_mean=float(sum(depth_datas.values()))/float(reflength)
    depth_covered_mean=float(np.mean(list(depth_datas.values())))
    covered_pos=len(depth_datas)
    covered_pct=float(covered_pos)/float(reflength)
    min_depth=0
    max_depth=0
    if reflength==len(depth_datas):
        min_depth=min(depth_datas.values())
    if len(depth_datas.values())>0:
        max_depth=max(depth_datas.values())
    ld.yaml_dumpDic(depthyaml,{'ref_length':reflength,'depth_datas':depth_datas})
    
    return {'depth_mean':depth_mean,'depth_covered_mean':depth_covered_mean,'covered_pos':covered_pos,'covered_pct':covered_pct,'min_depth':min_depth,'max_depth':max_depth}

def get_reads_nbr(stderr_file):
    with open(stderr_file,'r') as sfile:
        line = sfile.readline()
        try:
            return int(line.split()[0])
        except:
            return 0





with open(snakemake.input[0], 'rb') as f:
    bwt_results = pickle.load(f)

config=snakemake.config

outstats,depthfiles,reads_nbr=alignments_stats(bwt_results, snakemake.input[1])
combined_stats=cg.combine_alleles_stats(outstats,snakemake.input[1])
cg.determine_genotypes(combined_stats,snakemake.input[1],
                       snakemake.config['lbda'],
                       snakemake.config['mu'],
                       snakemake.config['sigma'],
                       snakemake.config['genotyp_alleleProb_THRLD'])
gd.normalize_depth(combined_stats,snakemake.input[1])

with open(snakemake.output[0], 'wb') as f:
    pickle.dump(combined_stats,f)
with open(snakemake.output[1], 'wb') as f:
    pickle.dump(reads_nbr,f)
