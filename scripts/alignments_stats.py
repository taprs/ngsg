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
import pysam

# def alignments_stats(bwt_results,pk_refdatabase):
def alignments_stats(bwt_results):
    # stdout_print("Do Alignments stats...",printInLog=True)
    # tps_stat=time.time()
    outstats={}
    # refinfo=ld.pickle_loadDic(pk_refdatabase)
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
        # depthfiles[sample]={}
        # c+=1
        # if tps1>0:
            # ellapsed_times.append(time.time()-tps1)
            # ETE=float(np.mean(ellapsed_times))*(sample_nbr-(c-1))
        # else:
            # ETE=0
        # stdout_print("\t{}/{} => do stats for {} -- ETE {} s".format(c,sample_nbr,sample,round(ETE,2)),printInLog=True)
        # tps1 = time.time()
        # for refID, idval in refdatas.items():
        refID="ref"
        idval=refdatas[refID]
        stats_file=idval['stats']
        reduced_stats_file=idval['reduced_stats']
        stats_results=parse_stats_file(stats_file)
        reduced_stats_results=parse_stats_file(reduced_stats_file, reduced=True)
        cov=extract_coverage(stats_results)
        reduced_cov=extract_coverage(reduced_stats_results,reduced=True)
        analys=analyse_bam(idval["bam_sorted"])
        seqs = stats_results.keys()
        reads_nbr[sample]=0
        for seqID in seqs:
            seqID=seqID.split('|')[0].strip() 
            reads_nbr[sample]+=stats_results[seqID]["stats"]["reads mapped"]
            if stats_results[seqID]["stats"]["reads mapped"]<=snakemake.config["readsMappedThrld"]:
                continue
    # stats_file, reduced_stats_file = snakemake.input[0:2]
    # idval["bam_sorted"], idval["reduced_sorted"] = snakemake.input[2:4]
            # reflength=refinfo[refID]['len_seq']
            # reflength=refinfo[refID]['len_seq']
            # stats_results=parse_stats_file(stats_file)
            # reduced_stats_results=parse_stats_file(reduced_stats_file)
            # depthyaml=os.path.join(idval['aln_path'],'depthinfos_{}'.format(idval['aln_name']))
            # reduced_depthyaml=os.path.join(idval['aln_path'],'reduced_depthinfos_{}'.format(idval['aln_name']))
            # depthfile=os.path.join(idval['aln_path'],'{}_SORTED.bam.depth'.format(idval['aln_name']))
            # reduced_depthfile=os.path.join(idval['aln_path'],'{}_reduced.bam.depth'.format(idval['aln_name']))
            # if stats_results['reads mapped']>snakemake.config['readsMappedThrld']:
            #     print(idval['aln_name'])
            #     # depthfiles[sample][refID]={'sorted':depthyaml,'reduced':reduced_depthyaml}
            #     outstats[sample][refID]={'stats':stats_results}
            #     outstats[sample][refID].update({'reduced_stats':reduced_stats_results})
            #     outstats[sample][refID].update(analyse_bam(idval['bam_sorted']))
            #     outstats[sample][refID].update({'ref_length':reflength})
            #     outstats[sample][refID].update({'coverage':extract_coverage(depthfile,depthyaml,reflength)})
            #     outstats[sample][refID].update({'reduced_coverage':extract_coverage(reduced_depthfile,reduced_depthyaml,reflength)})
                # outstats[sample][refID].update({'coverage':extract_coverage(idval['bam_sorted'],depthyaml,reflength)})
                # outstats[sample][refID].update({'reduced_coverage':extract_coverage(idval['reduced_sorted'],reduced_depthyaml,reflength)})
            # if sample not in reads_nbr.keys():
                # reads_nbr[sample]=get_reads_nbr(bwt_results[sample][refID]['errBWTfile'])
                # reads_nbr[sample]=get_reads_nbr(bwt_results[sample][refID]['errBWTfile'])
            # for seqID in stats_results.keys():
            # if stats_results[seqID]['reads mapped']>snakemake.config['readsMappedThrld']:
                # print(idval['aln_name'])
                # depthfiles[sample][refID]={'sorted':depthyaml,'reduced':reduced_depthyaml}
                # outstats[sample][seqID]={'stats':stats_results[seqID]}
                # outstats[sample][seqID].update({'reduced_stats':reduced_stats_results[seqID]})
                # outstats[sample][seqID].update(analyse_bam(idval['bam_sorted']))
                # outstats[sample][seqID].update({'ref_length':reflength})
                # outstats[sample][seqID].update({'coverage':extract_coverage(depthfile,depthyaml,reflength)})
                # outstats[sample][seqID].update({'coverage':extract_coverage(depthfile,depthyaml,reflength)})
                # outstats[sample][seqID].update({'reduced_coverage':extract_coverage(reduced_depthfile,reduced_depthyaml,reflength)})
                # outstats[sample][seqID].update({'coverage':extract_coverage(idval['bam_sorted'],depthyaml,reflength)})
                # outstats[sample][seqID].update({'reduced_coverage':extract_coverage(idval['reduced_sorted'],reduced_depthyaml,reflength)})
            outstats[sample][seqID]=analys[seqID]
            outstats[sample][seqID].update(stats_results[seqID])
            outstats[sample][seqID].update(reduced_stats_results[seqID])
            outstats[sample][seqID].update(cov[seqID])
            outstats[sample][seqID].update(reduced_cov[seqID])
        # if sample not in reads_nbr.keys():
        #     reads_nbr[sample]=get_reads_nbr(bwt_results[sample][seqID]['errBWTfile'])
        #         reads_nbr[sample]=get_reads_nbr(bwt_results[sample][seqID]['errBWTfile'])
    # stdout_print("END Alignments stats in {} s".format(round(time.time()-tps_stat,2)),printInLog=True)
    # return outstats,depthfiles
    # return outstats,depthfiles,reads_nbr
    return outstats,reads_nbr

def extract_coverage(stats, reduced=False):
    suffix= ("reduced_stats" if reduced else "stats")
    suffix2= ("reduced_coverage" if reduced else "coverage")
    depth_datas={}
    for r, ss in stats.items():
        s=ss[suffix]
        reflength=float(s["len_seq"])
        depth_mean=float(s["bases mapped"])/reflength
        covered_pos=float(s["covered_pos"])
        covered_pct=covered_pos/reflength
        depth_covered_mean=( float(s["bases mapped"])/covered_pos if covered_pos>0 else 0 )
        # min_depth=0
        # max_depth=0
        # if reflength==covered_pos:
        #     min_depth=1
        # if len(depth_datas.values())>0:
        #     max_depth=max(depth_datas.values())
        # ld.yaml_dumpDic(depthyaml,{'ref_length':reflength,'depth_datas':depth_datas})
        depth_datas[r]={suffix2:{'depth_mean':depth_mean,'depth_covered_mean':depth_covered_mean,'covered_pos':covered_pos,'covered_pct':covered_pct}}
    return depth_datas

def parse_stats_file(stats_file, reduced=False):
    suffix= ("reduced_stats" if reduced else "stats")
    filestats={}
    if os.path.exists(stats_file):
        with open(stats_file) as stfile:
            for line in stfile:
                # tmp=[v.strip() for v in line.split(':')]
                _=[v.strip() for v in line.split('\t')]
                _[0]=_[0].split('|')[0].strip() 
                filestats[_[0]]={suffix:{}}
                # keyid=tmp[0]
                # value=float(tmp[1])
                # filestats[keyid]=value
                filestats[_[0]][suffix]["len_seq"]=float(_[2])
                filestats[_[0]][suffix]["raw sequences"]=float(_[3])
                filestats[_[0]][suffix]["reads mapped"]=float(_[3])
                filestats[_[0]][suffix]["bases mapped"]=float(_[5])
                filestats[_[0]][suffix]["covered_pos"]=float(_[6])
                filestats[_[0]][suffix]["mismatches"]=( float(_[7]) if len(_)==8 else 0 )
            # if filestats['bases mapped']>0:
                if filestats[_[0]][suffix]['bases mapped']>0:
                    filestats[_[0]][suffix]["error rate"]=filestats[_[0]][suffix]['mismatches']/filestats[_[0]][suffix]['bases mapped']
                else:
                    filestats[_[0]][suffix]["error rate"]=0
                    # error_rate=filestats['mismatches']/filestats['bases mapped']
                # return {'raw sequences':filestats['raw total sequences'],'reads mapped':filestats['reads mapped'],'bases mapped':filestats['bases mapped'],'mismatches':filestats['mismatches'],'average length':filestats['average length'],'error rate':error_rate}
    # return {'reads mapped':0,'bases mapped':0,'mismatches':0,'average length':0,'error rate':0}
    return filestats

def analyse_bam(bamfile):
    #{'ratio_readsNoInDel':ratio_readsNoInDel,'readsNoInDelNbr':readsNoInDelNbr,'mismatchs_list':XM_list,'mismatch_mean':float(np.mean(XM_list))}

    # mismatchthrld=args.mismatchthrld
    mismatchthrld=snakemake.config["mismatchthrld"]

    outstats={}
    
    # bamdatas=bamanalyse.read_bam(bamfile)
    bamdatas_all=read_bam(bamfile)
    for r in bamdatas_all.keys():
        bamdatas=bamdatas_all[r]
        r=r.split('|')[0].strip() 
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
        # MQ_list=[v for v in bamdatas['MAPQ_list'] if v!=255]
        # if len(MQ_list)>0:
        #     mean_MQ = float(np.mean(MQ_list))
        # else:
        #     mean_MQ = 0
        # #prob all mapping positions are wrong = (1/10)^(PHRED_SCORE/10)
        # MQ_prob=0.1**(mean_MQ/10)
        
        # return {'readsNb':bamdatas['readsNb'],'ratio_readsNoInDel':bamdatas['ratio_readsNoInDel'],'readsNoInDelNbr':bamdatas['readsNoInDelNbr'],'mismatch_mean':bamdatas['mismatch_mean'],'distrib_XM_values':distrib_XM_values,'MQ_prob':MQ_prob}
        outstats[r]={'readsNb':bamdatas['readsNb'],'ratio_readsNoInDel':bamdatas['ratio_readsNoInDel'],'readsNoInDelNbr':bamdatas['readsNoInDelNbr'],'mismatch_mean':bamdatas['mismatch_mean'],'distrib_XM_values':distrib_XM_values}
    return outstats

# def extract_coverage(bamfile,depthyaml,reflength):
# def extract_coverage(depthfile,depthyaml,reflength):
#     # cmd = "{path}samtools depth {bam} | cut -f 2-".format(path=samtools_path,bam=bamfile)
#     # depthrstl = utils.run(args,"samtools depth", '', cmd, retOUT=True)
#     # depthstdout=depthrstl[0].decode('utf-8')
#     
#     with open(depthfile, 'r') as file:
#         depthstdout=file.readlines()
#
#     depth_datas={}
#     # for depthline in depthstdout.split('\n'):
#     for depthline in depthstdout:
#         depth_array=depthline.split('\t')
#         if len(depth_array)==2:
#             depth_datas[int(depth_array[0])]=int(depth_array[1])
#             
#     depth_mean=float(sum(depth_datas.values()))/float(reflength)
#     depth_covered_mean=float(np.mean(list(depth_datas.values())))
#     covered_pos=len(depth_datas)
#     covered_pct=float(covered_pos)/float(reflength)
#     min_depth=0
#     max_depth=0
#     if reflength==len(depth_datas):
#         min_depth=min(depth_datas.values())
#     if len(depth_datas.values())>0:
#         max_depth=max(depth_datas.values())
#     ld.yaml_dumpDic(depthyaml,{'ref_length':reflength,'depth_datas':depth_datas})
    
    # return {'depth_mean':depth_mean,'depth_covered_mean':depth_covered_mean,'covered_pos':covered_pos,'covered_pct':covered_pct,'min_depth':min_depth,'max_depth':max_depth}

# def get_reads_nbr(stderr_file):
#     with open(stderr_file,'r') as sfile:
#         line = sfile.readline()
#         try:
#             return int(line.split()[0])
#         except:
#             return 0

def read_bam(bamfilename):
    """
    This function extract informations from input bam file
    Input: bam file path
    Return a dictionary of per-chromosome dictionaries:
    'readsNb': reads number in bam
    'ratio_readsNoInDel': ratio of reads contains no InDel
    'readsNoInDelNbr': number of reads contains no InDel
    'mismatchs_list': list with all numbers of mismatches in the alignment
    'mismatch_mean': mean of mismatches in the alignment
    """
    inbam=pysam.AlignmentFile(bamfilename,'rb')
    #MD_list=[]

    #XM tag from bowtie2 manual:
    #XM:i:<N> The number of mismatches in the alignment. Only present if SAM record is for an aligned read.

    #MAPQ: MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. A value 255 indicates that the mapping quality is not available.
    XM_list={}
    XM_Pct={}
    MQ_list={}
    readsNoInDelNbr={}
    ratio_readsNoInDel={}
    readsNb={}
    outstats={}
    for r in inbam.references:
        XM_list[r]=[]
        XM_Pct[r]=[]
        MQ_list[r]=[]
        readsNoInDelNbr[r]=0
        ratio_readsNoInDel[r]=0
        readsNb[r]=0
    # XM_list=[]
    # XM_Pct=[]
    # MQ_list=[]
    # readsNoInDelNbr=0
    # ratio_readsNoInDel=0
    # readsNb=0
    # read_list=[r for r in inbam]
    # if len(read_list)>0:
        # for read in read_list:
    for read in inbam:
        r=read.reference_name
        tags=bamanalyse.tagsToDict(read.tags)
        #print(read.reference_start,tags['XM'],tags['MD'], read.query_name,read.cigarstring,containsInDel(read))
        # readsNb+=1
        readsNb[r]+=1
        if not bamanalyse.containsInDel(read):
            #MD_list.append(tags['MD'])
            XM_list[r].append(tags['XM'])
            XM_Pct[r].append(float(tags['XM'])/float(read.reference_length))
            MQ_list[r].append(read.mapping_quality)
    for r in inbam.references:
        if readsNb[r]==0:
            readsNoInDelNbr[r]=0
            ratio_readsNoInDel[r]=0
        else:
            readsNoInDelNbr[r]=len(XM_list[r])
            ratio_readsNoInDel[r]=float(readsNoInDelNbr[r])/float(readsNb[r])
            #print(readsNb,readsNoInDelNbr)
            # mutatedfraction={m:float(Counter(XM_list[r])[m])/float(readsNoInDelNbr[r]) for m in sorted(Counter(XM_list[r]).keys())}
        outstats[r]={'readsNb':readsNb[r],'ratio_readsNoInDel':ratio_readsNoInDel[r],'readsNoInDelNbr':readsNoInDelNbr[r],'mismatchs_list':XM_list[r],'MAPQ_list':MQ_list[r],'mismatch_mean':float(np.mean(XM_list[r]))}
    inbam.close()
    # return {'readsNb':readsNb,'ratio_readsNoInDel':ratio_readsNoInDel,'readsNoInDelNbr':readsNoInDelNbr,'mismatchs_list':XM_list,'MAPQ_list':MQ_list,'mismatch_mean':float(np.mean(XM_list))}
    return outstats





with open(snakemake.input[0], 'rb') as f:
    bwt_results = pickle.load(f)

config=snakemake.config

outstats,reads_nbr=alignments_stats(bwt_results)
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
