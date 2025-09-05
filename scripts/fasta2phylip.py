import sys
from Bio import SeqIO

inFasta = snakemake.input[0]
outPhy = snakemake.output[0]
fastaFile = SeqIO.parse(inFasta,"fasta")
phyOut = ""
seqLen=-1
N = 0
for rec in fastaFile:
        N += 1
        seq = str(rec.seq)
        phyOut += "{}\t{}\n".format(rec.id,seq)
        if seqLen>=0 and len(seq)!=seqLen:
                # addTextToLogFile("ERROR fasta sequences lentgh should be identical...")
                print("ERROR fasta sequences lentgh should be identical...", file=sys.stderr)
                # return False
                raise SystemExit(1)
        seqLen = len(seq)

with open(outPhy,"w") as of:
    of.write("{} {}\n".format(N,seqLen))
    of.write(phyOut)

