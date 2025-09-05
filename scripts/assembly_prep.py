import sys
import os
import pickle

sys.path.insert(0, snakemake.params["pymod"])
import load_dump as ld

# def get_reads_data(assembly_config,genotypfolder_sub):
# def get_reads_data(readsconfig,genotypfolder_sub):
#     # reads_data = ld.yaml_loadDic(os.path.join(genotypfolder_sub,assembly_config['readsconfig']))
#     # reads_data = ld.yaml_loadDic(os.path.join(genotypfolder_sub,readsconfig))
#     reads_data = ld.yaml_loadDic(readsconfig)
#     # if args.indivlist:
#     if snakemake.config["indivlist"]:
#         indiv_to_assemble=[]
#         with open(args.indivlist,'r') as indivfile:
#             for line in indivfile:
#                 if line!="":
#                     indiv_to_assemble.append(line.strip())
#         ret_reads_datas={indiv:datas for indiv,datas in reads_data.items() if indiv in indiv_to_assemble}
#
#         if len(ret_reads_datas.keys())>0:
#             return ret_reads_datas
#     return reads_data

# def get_reference(assembly_config,genotypfolder_sub):
    # return os.path.basename(assembly_config['ref']), ld.pickle_loadDic(os.path.join(genotypfolder_sub,assembly_config['ref']))
    # return os.path.basename(assembly_config['ref']), ld.pickle_loadDic(os.path.join(genotypfolder_sub,assembly_config['ref']))

def write_reference_fasta(out_fasta,phylo_ref,reference,add_paralogs=False):
    #'seq', 'len_seq', 'combined', 'Paralog', 'specie', 'grpRef'
    #'allelePart', 'grpRefPart', 'grpIDPart', 'haploIDPart', 'grpID', 'haploID'
    with open(out_fasta,'w') as f_alleles:
        with open(phylo_ref,'w') as f_alleles_phylo:
            for allele,datas in reference.items():
                if (not datas['Paralog'] or add_paralogs) and not datas['combined']:
                    f_alleles.write(">{}\n{}\n".format(allele,datas['seq']))
                    allele+="|Group={}|haploID={}".format(datas['grpID'],datas['haploID'])
                    if datas['Paralog']: allele+="_Paralog"
                    f_alleles_phylo.write(">{}\n{}\n".format(allele,datas['seq']))

# def assembly(reads_data,assembler,config,main_configs):
def assembly(reads_data,assembler,config):
    #usedKmerSize = "21,41,81"
    #cmd = "{path}spades.py -t {thrds} --careful -k {kmerS} -o {outf} {fastq} > {logf}".format(path=spades_folder,outf=spadesOutFolder,kmerS=usedKmerSize,fastq=inputSpadesFQ,thrds=spThrd,logf=spadesLogFile)
    # stdout_print("Launch Assembly using '{}'".format(assembler),printInLog=True)
    out_asm={}
    # cur_dir=os.getcwd()
    cur_dir="."
    asm_config=config['ASMtools'][assembler]
    # tools_folder=main_configs['tools_folder']
    # tools_folder = snakemake.params["tools_folder"]

    # asm_path=os.path.join(tools_folder,asm_config['path'])
    asm_path=asm_config['path']
    # if os.path.exists(asm_path):
        # c=0
    f=open(snakemake.output[0], "w")
    for indiv,datas in reads_data.items():
        out_assembly=os.path.join(config['outfolder'],indiv,"Assembly_{}".format(assembler))
        #log_assembly=os.path.join(out_assembly,"log_assembly")
        log_assembly=os.path.join("log_assembly")
        # dir_create(out_assembly)
        # os.makedirs(out_assembly, exists_ok=True)
        # if args.paired and len(datas['reads'])%2==0:
        if snakemake.config["paired"] and len(datas['reads'])%2==0:
            fwd_reads=asm_config['fastq_sep'].join(["{}{}".format(asm_config['options']['fwd'],os.path.join(cur_dir,fq)) for c,fq in enumerate(datas['reads']) if c%2==0])
            rev_reads=asm_config['fastq_sep'].join(["{}{}".format(asm_config['options']['rev'],os.path.join(cur_dir,fq)) for c,fq in enumerate(datas['reads']) if c%2==1])
            fastq_str=fwd_reads+asm_config['fastq_sep']+rev_reads
        else:
            fastq_str=asm_config['fastq_sep'].join(["{}{}".format(asm_config['options']['single'],os.path.join(cur_dir,fq)) for fq in datas['reads']])
            
        # options = asm_config['cmdline'].format(thrds=args.threaded,kmerS=asm_config['kmers'],outf=out_assembly,fastq=fastq_str,logf=log_assembly)
        options = asm_config['cmdline'].format(thrds=snakemake.config["threaded"],kmerS=asm_config['kmers'],outf=out_assembly,fastq=fastq_str,logf=log_assembly)
        # cmd = asm_path + " " + options
        #print(cmd)
        
        # out_contigs=os.path.join(out_assembly.format(assembler),asm_config['out_contigs'])
        out_contigs=os.path.join(out_assembly,asm_config['out_contigs'])
        
        #Create contig simlink:
        out_contigs_local=os.path.join("Assembly_{}".format(assembler),asm_config['out_contigs'])
        out_contigs_link="contigs.fa"

        # os.chdir(os.path.join(cur_dir,config['outfolder'],indiv))
        # if os.path.islink(out_contigs_link):
            # os.remove(out_contigs_link)
        # os.symlink(out_contigs_local,out_contigs_link)
        # os.chdir(cur_dir)
        out_contigs_link_abs=os.path.join(config['outfolder'],indiv,out_contigs_link)
        
        cmd = "mkdir -p " + out_assembly + " && " + \
              asm_path + " " + options + \
              " && ln -sf " + out_contigs_local + " " + os.path.join(out_asembly, indiv, out_contigs_link)
        f.writelines(cmd+"\n")

        # if not os.path.exists(out_contigs):
        #     os.chdir(os.path.join(cur_dir,out_assembly))
        #     os.system(cmd)
        #     os.chdir(cur_dir)

        asm_OK=True
        # if not os.path.exists(out_contigs_link_abs):
        #     asm_OK=False

        out_asm[indiv]={'contigs':out_contigs_link_abs,'asmfolder':os.path.join(config['outfolder'],indiv),'asm_OK':asm_OK}
        
        # c+=1
        # prog_pct=100*float(c)/float(len(reads_data.keys()))
        # progress_bar_pct("Assembly progress",prog_pct)
    # else:
        # exitProg("Tool not found: {}".format(asm_path))
    # print("")
    f.close()
    return out_asm

assembly_config=snakemake.config
# readsconfig=snakemake.config["readsconfig"]
readsconfig=snakemake.input[1]
genotypfolder_sub=os.path.dirname(snakemake.config["genotypfolder"])

config=snakemake.config
# out_asembly=os.path.join(args.genotypfolder,"Assembly_{}".format(args.asmsuffix))
out_asembly=os.path.join(config["genotypfolder"],"Assembly_{}".format(config["asmsuffix"]))
config['outfolder']=out_asembly
config['log_folder']=os.path.join(config['outfolder'],'logs')
config['reffolder']=os.path.join(config['outfolder'],"refs")
os.makedirs(config['outfolder'], exist_ok=True)
os.makedirs(config['reffolder'], exist_ok=True)


# pkfilename,reference_dic = get_reference(assembly_config,genotypfolder_sub)

reference_dic=ld.pickle_loadDic(snakemake.input[0])

# reads_data = get_reads_data(assembly_config,genotypfolder_sub)
# reads_data = get_reads_data(readsconfig,genotypfolder_sub)

reads_data=ld.yaml_loadDic(snakemake.input[1])

# ref_fasta=os.path.join(config['reffolder'],"{}.fasta".format(pkfilename))
ref_fasta=snakemake.output[2]
phylo_ref=snakemake.output[3]
write_reference_fasta(ref_fasta,phylo_ref,reference_dic,add_paralogs=snakemake.config["includeParalogsyass"])
# asm_contigs=assembly(reads_data,snakemake.config["assembler"],config,main_configs)
asm_contigs=assembly(reads_data,snakemake.config["assembler"],config)

with open(snakemake.output[1], "wb") as f:
    pickle.dump(asm_contigs,f)
