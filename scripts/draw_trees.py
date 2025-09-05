from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

tree = Phylo.read(snakemake.input[0],'newick')
Phylo.draw(tree)
plt.axis('off')
fig=plt.gcf()
fig.set_size_inches(tree.total_branch_length()*3, len(tree.get_terminals())*0.3)
plt.savefig(snakemake.output[0])
