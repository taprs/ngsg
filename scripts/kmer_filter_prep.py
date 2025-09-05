#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Adapted from mathieu-genete/NGSgenotyp commit e1376e8 by taprs.
No lines removed, only commented out and added some new code.
Uses code from the genotyp.py script.
Takes one entry from the reads list and applies kmer filter from a premade pickle file to it.
"""

import sys
import os
import pickle

sys.path.insert(0, snakemake.params["pymod"])
import load_dump as ld

def get_fileNoExt(filename,extlist):
    filebname=os.path.basename(filename)
    for ext in extlist:
        ext_index=filebname.upper().rfind(ext.upper())
        if ext_index>0:
            filebname=filebname[:ext_index]
    return filebname

def get_paired_reads_from_list(reads_list):
    PairedList=[]
    pair=[None,None]
    for c,reads in enumerate(reads_list):
        pair[c%2]=reads
        if (c%2)==1:
            PairedList.append(pair)
            pair=[None,None]
    return PairedList

def get_split_reads(infile,idx):
    bnameInfile=os.path.basename(infile)
    return os.path.join(outFQsplited,"{}_R{}_splited{}".format(bnameInfile[:bnameInfile.rfind('.')],idx,filesuffix))

with open(snakemake.input[0], 'rb') as file:
    readsConfig = pickle.load(file)

outFQfiltered=snakemake.config["OutFolders"]["out_filtered"] + "/need_split"
outFQsplited=snakemake.config['OutFolders']['out_splited']
filesuffix=".fastq"
# yaml_resume_file=os.path.join(outFQfiltered,config['KRF_yamlFile'])

out_template_FWD="{}FWD_filtered.fastq"
out_template_REV="{}REV_filtered.fastq"
out_template_Single="{}_filtered.fastq"

f = open(snakemake.output[0], "w")

# out_log_PK_KRF=os.path.join(config['log_folder'],"kmerRefFilter_pickle_{dt}.txt".format(dt=config['file_time']))
# out_log_KRF=os.path.join(config['log_folder'],"kmerRefFilter_log_{dt}.txt".format(dt=config['file_time']))
# jblst = utils.Jobslist("kmerRefFilter for all")
for sample,values in readsConfig.items():
# outKRF_dic[sample]={}
#Check if read already filtered is set
    extlist=[".fastq",".fq",".gz",".bz2"]
    if not values['params']['filtered']:
    # if not snakemake.params.filtered:
        readsConfig[sample]['params']['filtered']=True
        sample_reads=readsConfig[sample]['reads']
        readsConfig[sample]['reads']=[]
        # sample_reads=snakemake.input[1:]
        if values['params']['format']=='paired':
        # if snakemake.params.format=="paired":
            PairedList=get_paired_reads_from_list(sample_reads)
            for pair in PairedList:
                fq_fwd_bname=get_fileNoExt(pair[0],extlist)
                fq_rev_bname=get_fileNoExt(pair[1],extlist)

                out_filt_fwd=os.path.join(outFQsplited,out_template_FWD.format(fq_fwd_bname))
                out_filt_rev=os.path.join(outFQsplited,out_template_REV.format(fq_rev_bname))

                readsConfig[sample]['reads'].append(out_filt_fwd)
                readsConfig[sample]['reads'].append(out_filt_rev)
                if values["params"]["split"]:
                    readsConfig[sample]['reads']=[]
                    readsConfig[sample]['reads'].append(get_split_reads(out_filt_fwd,1))
                    readsConfig[sample]['reads'].append(get_split_reads(out_filt_rev,1))
                    readsConfig[sample]['reads'].append(get_split_reads(out_filt_fwd,2))
                    readsConfig[sample]['reads'].append(get_split_reads(out_filt_rev,2))
                # cmd="{pth}kmerRefFilter.py -nolog -y -k {ksize} -1 {FwdFq} -2 {RevFq} -i {pkref} -o {outdir} 2>>{log} 1>>{yamlKmer}".format(pth=kmerrefilter_path,FwdFq=pair[0],RevFq=pair[1],outdir=outFQfiltered,pkref=config['kmer_Ref_PickleFile'],log=out_log_KRF,yamlKmer=yaml_resume_file,ksize=kmerSize)
                cmd="{pth}kmerRefFilter.py -nolog -y -k {ksize} -1 {FwdFq} -2 {RevFq} -i {pkref} -o {outdir}\n".format(
                        # pth=kmerrefilter_path,
                        pth=snakemake.params.kmerrefilter_path,
                        FwdFq=pair[0],RevFq=pair[1],
                        # outdir=outFQfiltered,
                        outdir=outFQfiltered if values["params"]["split"] else outFQsplited,
                        pkref=snakemake.input[1],
                        # log=out_log_KRF,
                        # yamlKmer=yaml_resume_file,
                        ksize=snakemake.config["kmerSize"])
                # addTextToLogFile("\tkmerRefFilter paired: {}".format(cmd))
                # jblst.add_a_job(cmd,"kmerRefFilter paired {}".format(sample),out_filt_fwd)
                # print(cmd)
                f.write(cmd)
                # os.system(cmd)
                
        # So far single mode runs every read file separately, might be suboptimal?
        elif values['params']['format']=='single':
        # elif snakemake.params.format=="single":
            for fqfile in sample_reads:
                fq_bname=get_fileNoExt(fqfile,extlist)
                out_filt_single=os.path.join(outFQfiltered,out_template_Single.format(fq_bname))

                readsConfig[sample]['reads'].append(out_filt_single)
                if values["params"]["split"]:
                    readsConfig[sample]['reads']=[]
                    readsConfig[sample]['reads'].append(get_split_reads(out_filt_single,1))
                    readsConfig[sample]['reads'].append(get_split_reads(out_filt_single,2))
                # cmd="{pth}kmerRefFilter.py -nolog -y -k {ksize} -f {Fq} -i {pkref} -o {outdir} 2>>{log} 1>>{yamlKmer}".format(pth=kmerrefilter_path,Fq=fqfile,outdir=outFQfiltered,pkref=config['kmer_Ref_PickleFile'],log=out_log_KRF,yamlKmer=yaml_resume_file,ksize=kmerSize)
                cmd="{pth}kmerRefFilter.py -nolog -y -k {ksize} -f {Fq} -i {pkref} -o {outdir}\n".format(
                        # pth=kmerrefilter_path,
                        pth=snakemake.params.kmerrefilter_path,
                        Fq=fqfile,
                        # outdir=outFQfiltered,
                        outdir=outFQfiltered + ("/need_split" if values["params"]["split"] else ""),
                        # pkref=config['kmer_Ref_PickleFile'],
                        pkref=snakemake.input[1],
                        # log=out_log_KRF,
                        # yamlKmer=yaml_resume_file,
                        ksize=snakemake.config["kmerSize"])
                addTextToLogFile("\tkmerRefFilter single: {}".format(cmd))
                # jblst.add_a_job(cmd,"kmerRefFilter paired {}".format(sample),out_filt_single)
                # print(cmd)
                f.write(cmd)
                # os.system(cmd)

f.close()
# thread_list=[]
# if args.tinyverbose:
# t1 = threading.Thread(target=progress_bar,args=['KmerRefFilter progress',len(jblst.get_joblist()),thread_list])
# t1.start()
#
# aunch kmerRefFilter process
# utils.trun(args, jblst,thrd=thread_list)
#
# if args.tinyverbose:
# wait_progressbar_finished(t1)

# stdout_print(" -- END kmerRefFilter --",printInLog=True)
print(" -- END kmerRefFilter --",file=sys.stdout)

# return readsConfig
with open(snakemake.output[1], 'wb') as f:
    pickle.dump(readsConfig,f)

ld.yaml_dumpDic(snakemake.output[2],readsConfig)

