from ete3 import Tree

out_nearest={}
# for sample,phyml_tree in phyml_samples.items():
phyml_tree=snakemake.input[0]
sample=snakemake.wildcards["sample"]
tree = Tree(phyml_tree)
out_nearest[sample]={}
for leaf in tree:
    if 'contig' in leaf.name:
        prev_nod=leaf
        end_while=False
        out_nearest[sample][leaf.name]=[]
        while True:
            prev_nod=prev_nod.up
            for child_leaf in prev_nod:
                if 'contig' not in child_leaf.name:
                    out_nearest[sample][leaf.name].append((child_leaf.name,tree.get_distance(leaf,child_leaf)))
                    end_while=True

            if prev_nod.is_root() or end_while:
                break

# out_txt="#Contigs closest Alleles\n"
# out_txt+="#\tContig\tAllele\tdistance\n"
out_txt=""
for sample,contigs in out_nearest.items():
    out_txt+="sample: {}\n".format(sample)
    for contig,alleles in contigs.items():
        if len(alleles)>0:
            closest_allele=sorted(alleles, key = lambda x:x[1])[0]
            out_txt+="\t{}\t{}\t{}\n".format(contig,closest_allele[0],closest_allele[1])

# out_closest_alleles=os.path.join(config['outfolder'],"contigs_closest_alleles.txt")
out_closest_alleles=snakemake.output[0]
# stdout_print("Write contigs closests alleles file: {}".format(out_closest_alleles),printInLog=True)
with open(out_closest_alleles,'w') as outc:
    outc.write(out_txt)
