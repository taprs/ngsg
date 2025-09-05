#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pickle

sys.path.insert(0, snakemake.params["pymod"])
import load_dump as ld


# def run_bowtie2(readsConfig,yaml_index,refexclude,alignmentMode,alignmentSensitivity,pairedAln):
def run_bowtie2(readsConfig,alignmentMode,alignmentSensitivity):
    # addTextToLogFile("-- ALIGNMENTS WITH BOWTIE2 --")
    # index_dic=ld.yaml_loadDic(yaml_index)
    index_dic=ld.yaml_loadDic(snakemake.input[1])
    aln_name_template="Aln_{refID}_VS_{sample}"
    align_outdic={}
    errBWTfile_list=[]
    # jblst = utils.Jobslist("bowtie2 for all")

    # if len(refexclude)>0:
        # stdout_print("Exclude reference file set",printInLog=True)
    # for xref in refexclude:
        # stdout_print("\t exclude reference: {}".format(xref),printInLog=True)
        
    f=open(snakemake.output[0], 'w')
    for sample,values in readsConfig.items():
        pairedAln = values["params"]["format"]=="paired"
        sample_bwt_path=os.path.join(config['OutFolders']['bwt_dirResult'],sample)
        # dir_create(sample_bwt_path)
        os.makedirs(sample_bwt_path, exist_ok=True)
        
        #Paired or single alignments (ONLY SINGLE FOR NOW)
        # if pairedAln and len(values['reads'])>1 and len(values['reads'])%2==0:
        #     # stdout_print("Do paired alignments",printInLog=True)
        #     forward_files=[s for i,s in enumerate(values['reads']) if i%2==0]
        #     reverse_files=[s for i,s in enumerate(values['reads']) if i%2!=0]
        #     samples="-1 {fwd} -2 {rev}".format(fwd=','.join(forward_files),rev=','.join(reverse_files))
        #     flags="-f 0x3"
        # else:
        # stdout_print("Do single alignments",printInLog=True)
        samples = "-U {}".format(','.join(values['reads']))
        flags="-F 0x4"
            
        align_outdic[sample]={}
        for refID,indexval in index_dic.items():
            # if refID not in refexclude:
            aln_name=aln_name_template.format(sample=sample,refID=refID)
            aln_path=os.path.join(sample_bwt_path,aln_name)
            # dir_create(aln_path)
            # os.makedirs(aln_path, exist_ok=True)

            errBWTfile=os.path.join(aln_path, "stderr_{}".format(aln_name))
                    # errBWTfile_list.append(errBWTfile)
            outbam=os.path.join(aln_path,"{}.bam".format(aln_name))
            reducedbam=os.path.join(aln_path,"{}_reduced.bam".format(aln_name))
            reducedbam_index=os.path.join(aln_path,"{}_reduced.bam.bai".format(aln_name))
            outbam_sorted=os.path.join(aln_path,"{}_SORTED.bam".format(aln_name))
            outbam_index=os.path.join(aln_path,"{}_SORTED.bam.bai".format(aln_name))
            mapped_unmapped_out = os.path.join(aln_path, "mapped_{}".format(aln_name))
            outstats=os.path.join(aln_path,"{}_samtools_stats".format(aln_name))
            reduced_outstats=os.path.join(aln_path,"{}_samtools_stats_reduced".format(aln_name))
                        
                    # bwt_cmd="{bpath}bowtie2 --phred33 --{alnmode} --{alnsens} -x {RefIndex} {samples} 2>{err_bowtie} | {path}samtools view -S -b {flags} - > {bamFile}".format(bpath=bowtie2_path,path=samtools_path,RefIndex=indexval['index'],samples=samples, err_bowtie=errBWTfile, bamFile=outbam,flags=flags,alnmode=alignmentMode,alnsens=alignmentSensitivity)
            bwt_cmd="mkdir -p {alnpath} && {bpath}bowtie2 --phred33 --{alnmode} --{alnsens} -x {RefIndex} {samples} 2>{err_bowtie} | {path}samtools view -S -b {flags} - | {path}samtools sort - > {bamFile}".format(
                    alnpath=aln_path,
                    bpath=snakemake.params.bowtie2_path,
                    path=snakemake.params.samtools_path,
                    RefIndex=indexval["index"],
                    samples=samples, 
                    err_bowtie=errBWTfile, 
                    # bamFile=outbam,
                    bamFile=outbam_sorted,
                    flags=flags,
                    alnmode=config["alignmentMode"],
                    alnsens=config["alignmentSensitivity"])
            f.write(bwt_cmd + "\n")
            align_outdic[sample][refID]={'aln_name':aln_name,'aln_path':aln_path,'aln_path':aln_path,'errBWTfile':errBWTfile,'bam':outbam,'orig_sorted':outbam_sorted,'reduced_sorted':reducedbam,'reduced_index':reducedbam_index,'bam_sorted':outbam_sorted,'bam_index':outbam_index,'orig_index':outbam_index,'bwt_cmd':bwt_cmd,'stats':outstats,'orig_stats':outstats,'reduced_stats':reduced_outstats,'align_ok':False}
            # From read_filter_bam function
            if snakemake.config["statsmismatchthrld"]:
                align_outdic[sample][refID]['bam_sorted']=reducedbam
                align_outdic[sample][refID]['bam_index']=reducedbam_index
                align_outdic[sample][refID]['stats']=reduced_outstats
                    # addTextToLogFile("bowtie2 command for {} : {}".format(aln_name,bwt_cmd))
                    # jblst.add_a_job(bwt_cmd,"bowtie2 unpaired {}".format(aln_name),outbam)
    f.close()

    # align_nbr=sum(len(v) for v in align_outdic.values())

    # stdout_print("Launch {} alignments".format(align_nbr),printInLog=True)

    #Launch all bowtie2 commands

    # thread_list=[]
    # #show progress bar if verbose option is set
    # if args.tinyverbose:
    #     t1 = threading.Thread(target=progress_bar,args=['Alignments progress',len(errBWTfile_list),thread_list])
    #     t1.start()
    #
    # utils.trun(args, jblst,thrd=thread_list)
    #
    # if args.tinyverbose:
    #     wait_progressbar_finished(t1)
    #
    # #check alignments completion
    # errnbr,align_outdic=check_alignments(align_outdic)
    # 
    # stdout_print("Alignments results: {} completed / {} misaligned".format(align_nbr-errnbr,errnbr),printInLog=True)
    # if errnbr>0:
    #     stdout_print("\tSee log file for alignments errors".format(align_nbr-errnbr,errnbr))
    #
    # stdout_print("\nSort, index and stats for bam files")
    #
    # list_stats=[]
    # list_stats_reduced=[]
    # list_index=[]
    # list_index_reduced=[]
    # list_sort=[]
    # reads_nbr={}
    # for sample,refAlnDatas in align_outdic.items():
    #     for refID, aligndatas in refAlnDatas.items():
    #         if align_outdic[sample][refID]['align_ok']:
    #             list_stats.append((align_outdic[sample][refID]['bam_sorted'],align_outdic[sample][refID]['stats']))
    #             list_stats_reduced.append((align_outdic[sample][refID]['reduced_sorted'],align_outdic[sample][refID]['reduced_stats']))
    #             list_index.append((align_outdic[sample][refID]['bam_sorted'],align_outdic[sample][refID]['bam_index']))
    #             list_index_reduced.append((align_outdic[sample][refID]['reduced_sorted'],align_outdic[sample][refID]['reduced_index']))
    #             list_sort.append((align_outdic[sample][refID]['bam'],align_outdic[sample][refID]['bam_sorted']))
                # if sample not in reads_nbr.keys():
                    # reads_nbr[sample]=get_reads_nbr(align_outdic[sample][refID]['errBWTfile'])

#                 
#     launch_samtools_sort(list_sort)
#     launch_samtools_index(list_index)
#
#     launch_samtools_stats(list_stats)
#     
#     #Filtrer avant samtools stats => à revoir
#     align_outdic=read_filter_bam(align_outdic,args.mismatchthrld)
#     
#     launch_samtools_index(list_index_reduced)
#                 
#     launch_samtools_stats(list_stats_reduced)
#     
    with open(snakemake.output[1], 'wb') as f:
        pickle.dump(align_outdic, f)
    # with open(snakemake.output[2], 'wb') as f:
    #     pickle.dump(reads_nbr, f)
    # return align_outdic,reads_nbr
#
# def get_reads_nbr(stderr_file):
#     with open(stderr_file,'r') as sfile:
#         line = sfile.readline()
#         try:
#             return int(line.split()[0])
#         except:
#             return 0
#     
# def read_filter_bam(align_outdic,mutThrld):
#     stdout_print("Start reduce",printInLog=True)
#     for sample, refdatas in align_outdic.items():
#         for refID, idval in refdatas.items():
#             bam_file=idval['bam_sorted']
#             reduced_bam=idval['reduced_sorted']
#             
#             if not os.path.exists(reduced_bam):
#                 filtInDel=(not args.keepindel)
#                 bamanalyse.filter_bam(bam_file,reduced_bam,mutThrld,filterInDel=filtInDel)
#             
#             if args.statsmismatchthrld:
#                 align_outdic[sample][refID]['bam_sorted']=reduced_bam
#                 align_outdic[sample][refID]['bam_index']=idval['reduced_index']
#                 align_outdic[sample][refID]['stats']=idval['reduced_stats']
#
#     stdout_print("END reduce",printInLog=True)
#                 
#     return align_outdic
#
# def launch_samtools_stats(list_stats,title=""):
#     stats_jblst = utils.Jobslist("bam stats for all")
#     for bamfile,statsfile in list_stats:
#         #stats command
#         stats_cmd = "{path}samtools stats {inbam} | grep ^SN | cut -f2-3 > {outstats}".format(path=samtools_path,inbam=bamfile,outstats=statsfile)
#         stats_jblst.add_a_job(stats_cmd,"samtools stats {}".format(bamfile),statsfile)
#
#     thread_list=[]
#     if args.tinyverbose:
#         t3 = threading.Thread(target=progress_bar,args=['Samtools stats progress {}'.format(title),len(stats_jblst.get_joblist()),thread_list])
#         t3.start()
#         
#     stdout_print("Launch samtools stats {}".format(title),printInLog=True)
#     utils.trun(args, stats_jblst, thrd=thread_list)
#
#     if args.tinyverbose:
#         wait_progressbar_finished(t3)
#
# def launch_samtools_index(list_index,title=""):
#     index_jblst = utils.Jobslist("bam index for all")
#     
#     for bamfile,indexfile in list_index:
#         #index command
#         index_cmd= "{path}samtools index {inbam}".format(path=samtools_path,inbam=bamfile)
#         index_jblst.add_a_job(index_cmd,"samtools index {}".format(bamfile),indexfile)
#
#     thread_list=[]
#     if args.tinyverbose:
#         t2 = threading.Thread(target=progress_bar,args=['Samtools index progress {}'.format(title),len(index_jblst.get_joblist()),thread_list])
#         t2.start()
#         
#     stdout_print("Launch samtools index {}".format(title),printInLog=True)
#     utils.trun(args, index_jblst, thrd=thread_list)
#
#     if args.tinyverbose:
#         wait_progressbar_finished(t2)
#         
# def launch_samtools_sort(list_sort,title=""):
#     sort_jblst = utils.Jobslist("bam sort for all")
#
#     for bamfile,sortfile in list_sort:
#         #sort command
#         sort_cmd = "{path}samtools sort {inbam} -o {outbam}".format(path=samtools_path,inbam=bamfile,outbam=sortfile)
#         sort_jblst.add_a_job(sort_cmd,"samtools sort {}".format(bamfile),sortfile)
#
#     thread_list=[]
#     if args.tinyverbose:
#         t1 = threading.Thread(target=progress_bar,args=['Samtools sort progress {}'.format(title),len(sort_jblst.get_joblist()),thread_list])
#         t1.start()
#         
#     stdout_print("Launch samtools sort {}".format(title),printInLog=True)
#     utils.trun(args, sort_jblst, thrd=thread_list)
#
#     if args.tinyverbose:
#         wait_progressbar_finished(t1)
#     
with open(snakemake.input[0], 'rb') as file:
    readsConfig = pickle.load(file)

config=snakemake.config

run_bowtie2(readsConfig, snakemake.config["alignmentMode"], snakemake.config["alignmentSensitivity"])
