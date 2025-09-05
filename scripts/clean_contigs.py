import os
from Bio import SeqIO
import pickle

def clean_contigs_file(asm_contigs):
    out_asm_contigs={}
    for sample,datas in asm_contigs.items():
        contigs_file=datas['contigs']
        cleaned_contigs=os.path.join(datas['asmfolder'],"cleaned_contigs.fa")
        datas['cl_contigs']=cleaned_contigs
        if os.path.exists(contigs_file):
            out_asm_contigs[sample]=datas
            if not os.path.exists(cleaned_contigs):
                fasta = SeqIO.parse(contigs_file,'fasta')
                with open(cleaned_contigs,'w') as outcontigs:
                    for rec in fasta:
                        # if len(rec.seq)>=args.contigminlen and (not args.contigmaxlen or len(rec.seq)<=contigmaxlen):
                        if len(rec.seq)>=snakemake.config["contigminlen"] and (not snakemake.config["contigmaxlen"] or len(rec.seq)<=snakemake.config["contigmaxlen"]):
                            outcontigs.write(rec.format('fasta'))
    return out_asm_contigs

with open(snakemake.input[0], "rb") as f:
    asm_contigs=pickle.load(f)

out_asm_contigs = clean_contigs_file(asm_contigs)

with open(snakemake.output[0], "wb") as f:
    pickle.dump(out_asm_contigs, f)
