#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Adapted from mathieu-genete/NGSgenotyp commit e1376e8 by taprs.
No lines removed, only commented out and added some new code.
Uses code from genotyp.py.
Takes fasta reference (1st argument) and writes a directory of single-sequence fasta files (2nd argument) for future indexing and mapping.
Maybe replace with seqkit split?..
"""

import os
import gzip
import pickle
import sys
sys.path.insert(0, snakemake.params["pymod"])
import load_dump as ld

def index_references(pk_refdatabase, outdir):
    os.makedirs(outdir, exist_ok=True)
    # stdout_print(" -- Index reference file --",printInLog=True)
    print(" -- Split reference file --",file=sys.stdout)
    index_dic={}
    bwt_out = os.path.join(config['log_folder'],config['bwt_build_out'])
    yaml_out=snakemake.output[1]

    # if os.path.exists(yaml_out) and not args.force:
        # return yaml_out
    
    #Create single fasta files to index
    refDB_dic=ld.pickle_loadDic(pk_refdatabase)
    for seqID,seqVal in refDB_dic.items():
        singlefasta=os.path.join(config['OutFolders']['single_fasta_ref_fld'],"{}.fa".format(seqID))
        singlefasta=os.path.join(outdir,"{}.fa".format(seqID))
        indexfile=os.path.join(config['OutFolders']['bwt_IndexdirResult'],seqID)
        if not seqVal['combined']:
            index_dic[seqID]={'fasta':singlefasta,'index':indexfile}
            
        if (not os.path.exists(singlefasta) or args.force) and not seqVal['combined']:
            with open(singlefasta,'w') as outfasta:
                outfasta.write(">{}\n{}".format(seqID,seqVal['seq']))
                
    #Index all fasta files
    # jblst = utils.Jobslist("bowtie2-build for all")
    # for seqID,indexval in index_dic.items():
    #     target = "{}{}".format(indexval['index'],'.1.bt2')
    #     cmd = '{bpath}bowtie2-build {reffile} {refname} >>{bwtout}'.format(bpath=bowtie2_path,reffile = indexval['fasta'],refname = indexval['index'],bwtout = bwt_out)
    #     addTextToLogFile("\tBowtie2-build: {}".format(cmd))
    #     jblst.add_a_job(cmd,"bowtie2-build {}".format(seqID),target)
    ld.yaml_dumpDic(yaml_out,index_dic)
    # utils.trun(args, jblst)
    # return yaml_out
    #

config=snakemake.config

index_references(snakemake.input[0], snakemake.output[0])
