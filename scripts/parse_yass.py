import os
import pickle
import re
from Bio import SeqIO

sys.path.insert(0, snakemake.params["pymod"])
import load_dump as ld

def parse_contigs_infos(asm_contigs,assembler):
    out_infos={}
    cov_regex=config['ASMtools'][assembler]['regex_cov']
    contig_nbr_regex=config['ASMtools'][assembler]['regex_contig_nbr']
    for sample,datas in asm_contigs.items():
        out_infos[sample]={}
        contigs_file=datas['contigs']
        fasta=SeqIO.parse(contigs_file,'fasta')
        for rec in fasta:
            seq_id=str(rec.description)
            COV_rslt=re.search(cov_regex,seq_id)
            contig_nbr=re.search(contig_nbr_regex,seq_id).group(0)
            if COV_rslt is not None:
                cov_val=float(COV_rslt.group(0)) 
            else:
                cov_val=0.0
                
            out_infos[sample][contig_nbr]={'description':seq_id,'length':len(rec.seq),'cov':cov_val}
    return out_infos

# def parse_yass_results(yass_results,contigs_infos,reference_dic,assembler,filter_paralogs=True):
def parse_yass_results(contigs_infos,reference_dic,assembler,filter_paralogs=True):
    out_yass={}
    # outfile_yass=os.path.join(config['outfolder'],"raw_yass_datas.pk")
    strand_val={1:"+",-1:"-"}
    # for sample,yass_stdout in yass_results.items():
    sample=snakemake.wildcards["sample"]
    yass_stdout=snakemake.input[0]
    yass_rslt=parse_yass_stdout(yass_stdout,assembler)
    out_rslt=[]
    for datas in yass_rslt:
        ref_len=reference_dic[datas['Subject_id']]['len_seq']
        datas['ref_length']=ref_len
        ident_cov=datas['percent_identity']*(datas['alignment_length']/ref_len)
        datas['coverage']=(datas['alignment_length']/ref_len)
        datas['ident_cov']=ident_cov
        datas['contig_depth']=contigs_infos[sample][str(datas['Query_id'])]['cov']
        datas['strand']=strand_val[int((datas['s_end']-datas['s_start'])/abs(datas['s_end']-datas['s_start']))]

        #Apply Thresholds
        if all([eval(str(datas[param])+thrld) for param,thrld in config['YASSstep']['Thresholds'].items()]):
            out_rslt.append(datas)
    out_yass[sample]=out_rslt

    # ld.pickle_dumpDic(outfile_yass,out_yass)
    # return outfile_yass
    return out_yass

def parse_yass_stdout(yass_stdout,assembler):
    #query => contigs
    #subject => ref allele
    out_rslt=[]
    fields=['Query_id','Subject_id', 'percent_identity', 'alignment_length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'e-value', 'bit_score']
    contig_nbr_regex=config['ASMtools'][assembler]['regex_contig_nbr']
    with open(yass_stdout,'r') as yassfile:
        for line in yassfile:
            if not line.startswith("#"):
                tmp=line.strip().split('\t')
                tmp[0]=re.search(contig_nbr_regex,tmp[0]).group(0)
                out_rslt.append({fields[i]:string_to_number(val) for i,val in enumerate(tmp)})
    return out_rslt

def string_to_number(text):
    if bool(re.search("^[+-]?([0-9]*[.])?[0-9]+",text)):
        if "." not in text and "e" not in text:
            return int(text)
        return float(text)
    return text

# def get_oriented_contigs(raw_yass_file,contigs_infos,asm_contigs):
def get_oriented_contigs(yass_datas,contigs_infos,asm_contigs):
    # yass_datas=ld.pickle_loadDic(raw_yass_file)
    out_yass_contigs={}
    for sample,all_datas in yass_datas.items():
        # contigs= {rec.description:rec.seq for rec in SeqIO.parse(asm_contigs[sample]['cl_contigs'],'fasta')}
        contigs= {rec.description:rec.seq for rec in SeqIO.parse(snakemake.input[1],'fasta')}
        # stdout_print("\n====> Extract contigs from yass -- sample: {} -- total number of contigs: {}".format(sample,len(contigs)),printInLog=True)
        rc_count=0
        unique_contigs_list=set([(str(datas['Query_id']),datas['strand'],datas['contig_depth']) for datas in all_datas])
        phylo_cut={}
        for datas in all_datas:
            q_id=str(datas['Query_id'])
            if q_id not in phylo_cut.keys():
                phylo_cut[q_id]={'start':-1,'end':0}
            q_start=datas['q_start']-1
            q_end=datas['q_end']

            if q_start<phylo_cut[q_id]['start'] or phylo_cut[q_id]['start']==-1:
                phylo_cut[q_id]['start']=q_start

            if q_end>phylo_cut[q_id]['end']:
                phylo_cut[q_id]['end']=q_end
        extracted_contigs={}
        for contig_id,strand,contig_depth in unique_contigs_list:
            contig_description=contigs_infos[sample][contig_id]['description']
            contig_seq=str(contigs[contig_description])
            phylo_seq=contig_seq[phylo_cut[contig_id]['start']:phylo_cut[contig_id]['end']]
            if strand=='-':
                rc_count+=1
                # stdout_print("\t{} reverse complement: {}".format(rc_count,contig_description),printInLog=True)
                contig_seq=reverse_complement(contig_seq)
                phylo_seq=reverse_complement(phylo_seq)
            extracted_contigs_id="{}|contig={}|l={}|d={}|s={}".format(sample,contig_id,len(contig_seq),contig_depth,strand)
            if strand=="-":
                extracted_contigs_id+="_RC"
            extracted_contigs[extracted_contigs_id]={'c_seq':contig_seq,'phylo_seq':phylo_seq}
        # stdout_print("{} contigs extracted <====".format(len(extracted_contigs)),printInLog=True)
        # outfasta=os.path.join(asm_contigs[sample]['asmfolder'],"yass_oriented_contigs.fa")
        # outfasta_phylo=os.path.join(asm_contigs[sample]['asmfolder'],"phylo_truncated_contigs.fa")
        outfasta, outfasta_phylo=snakemake.output[0:2]
        # out_yass_contigs[sample]={'contigs':outfasta,'c_phylo':outfasta_phylo,'nb_contigs':len(extracted_contigs.keys())}
        with open(outfasta_phylo,'w') as outphylo:
            with open(outfasta,'w') as outf:
                for c_id,seq in extracted_contigs.items():
                    outf.write(">{}\n{}\n".format(c_id,seq['c_seq']))
                    outphylo.write(">{}\n{}\n".format(c_id,seq['phylo_seq']))
    # return out_yass_contigs

def reverse_complement(seq):
    rev_seq=seq[::-1].upper()
    bases={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    rev_comp=""
    for b in rev_seq:
        if b in bases.keys():
            rev_comp+=bases[b]
        else:
            rev_comp+="N"
    return rev_comp


reference_dic = ld.pickle_loadDic(snakemake.input[3])

config=snakemake.config
with open(snakemake.input[2], 'rb') as f:
    asm_contigs=pickle.load(f)

contigs_infos=parse_contigs_infos(asm_contigs,snakemake.config["assembler"])
raw_yass_file=parse_yass_results(contigs_infos,reference_dic,snakemake.config["assembler"])
oriented_contigs=get_oriented_contigs(raw_yass_file,contigs_infos,asm_contigs)
