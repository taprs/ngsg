#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Adapted from mathieu-genete/NGSgenotyp commit e1376e8 by taprs.
No lines removed, only commented out and added some new code.
Uses code from genotyp.py, load_dump.py, and checkDB.py scripts.
Checks fasta reference (1st argument) and writes a pickle file (2nd argument) for kmer filtering.
"""

import sys
import os


# def run_kmerRefFilter(kmerSize,readsConfig):
def run_kmerRefFilter(kmerSize, refdatabase, outfile):
    # readsConfig = load_readsList(readsinfo)
    print(" -- kmerRefFilter --", file=sys.stdout)

    # outKRF_dic={}
    
    # outFQfiltered=config['OutFolders']['out_filtered']
    # yaml_resume_file=os.path.join(outFQfiltered,config['KRF_yamlFile'])
    
    # out_template_FWD="{}FWD_filtered.fastq"
    # out_template_REV="{}REV_filtered.fastq"
    # out_template_Single="{}_filtered.fastq"

    # out_log_PK_KRF=os.path.join(config['log_folder'],"kmerRefFilter_pickle_{dt}.txt".format(dt=config['file_time']))
    # out_log_KRF=os.path.join(config['log_folder'],"kmerRefFilter_log_{dt}.txt".format(dt=config['file_time']))

    #Create kmerRefFilter pickle database
    # stdout_print("\tconstruct kmer DB for filtering : {}".format(config['kmer_Ref_PickleFile']),printInLog=True)
    print("construct kmer DB for filtering", file=sys.stdout)
    # cmd="{pth}kmerRefFilter.py -nolog -P {cpu} -p {outpickle} -k {ksize} -r {fastaref} -o . 2>{log} 1>&2".format(pth=kmerrefilter_path,outpickle=config['kmer_Ref_PickleFile'],cpu=config['MaxParallelsJobs'],fastaref=args.refdatabase,log=out_log_PK_KRF,ksize=kmerSize)
    cmd="{pth}kmerRefFilter.py -nolog -P {cpu} -p {outpickle} -k {ksize} -r {fastaref} -o .".format(
            pth=snakemake.params.kmerrefilter_path,
            outpickle=outfile,
            cpu=snakemake.config['MaxParallelsJobs'],
            fastaref=refdatabase,
            ksize=snakemake.config["kmerSize"])
    # if not os.path.exists(config['kmer_Ref_PickleFile']) or args.force:
    # if not snakemake.config["force"]:
    os.system(cmd)
    

run_kmerRefFilter(snakemake.config["kmerSize"], snakemake.input[0], snakemake.output[0])
