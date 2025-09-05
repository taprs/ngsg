import re
import os
import sys

#link conda libraries
os.environ["LD_LIBRARY_PATH"] = os.environ["CONDA_PREFIX"] + "/lib:" + ( os.environ["LD_LIBRARY_PATH"] if "LD_LIBRARY_PATH" in os.environ.keys() else "" )

# Supply configfile via --configfile or uncomment the line below
configfile: workflow.basedir + "/config.yaml"

# Snakefile must stay in the ngsg folder!

App_Folder = workflow.basedir
Features_Path = os.path.join(App_Folder,"features/")
pymod_path = os.path.join(App_Folder,"pymod/")
tools_folder=os.path.join(App_Folder,'tools')
readsconfig=os.path.join(config['OutFolders']['configs_fld'],config['readsconfig'])

fastarefbase=os.path.splitext(os.path.basename(config["fastaref"]))[0]

outdir="Assembly_" + config["asmsuffix"] + "/"

# Parse sample names from reads file
samp=[]
with open(config["readsinfo"], "r") as f:
  for i,l in enumerate(f.readlines()):
    if (i==0):
      continue
    samp.append(l.split(",")[0])

### genotyp

rule genotyp:
  input:
    "Results/xls_samples_index.txt",
    "Results/putative_alleles.txt",
    "Results/putative_groups.txt",
    "Results/samtools_stats.yaml",
    "configs/asm_config.yaml"

rule kmer_filter_index:
  input:
    config["fastaref"]
  params:
    kmerrefilter_path=os.path.join(tools_folder, config["kmerRefFilter"]["bin_folder"])
  output:
    "Refs/KRFDB_" + fastarefbase + "_k" + str(config["kmerSize"]) + ".pk"
  script:
    "scripts/kmer_filter_index.py"

rule kmer_filter_prep:
  input:
    "readsConfig.pk",
    "Refs/KRFDB_" + fastarefbase + "_k" + str(config["kmerSize"]) + ".pk"
  params:
    kmerrefilter_path=os.path.join(tools_folder, config["kmerRefFilter"]["bin_folder"]),
    pymod=pymod_path
  output:
    "kmerRefFilter_commands.sh",
    "readsConfig_filt.pk",
    readsconfig
  script:
    "scripts/kmer_filter_prep.py"

rule parse_readsinfo:
  input:
    config["readsinfo"]
  output:
    "readsConfig.pk"
  script:
    "scripts/parse_readsinfo.py"

checkpoint kmer_filter_run:
  input:
    "kmerRefFilter_commands.sh"
  output:
    directory("FQ_filtered"),
    directory("FQ_splited"),
    directory("FQ_filtered/need_split")
  threads: workflow.cores * 0.75
  shell:
    "mkdir -p {output} && parallel -j {threads} -a {input}"

rule split_reads:
  input:
    "FQ_filtered/need_split/{sample}.fastq"
  params:
    compress=False
  output:
    "FQ_splited/{sample}_R1_splited.fastq",
    "FQ_splited/{sample}_R2_splited.fastq"
  script:
    "scripts/split_reads.py"

rule parse_reference_db:
  input:
    config["fastaref"]
  output:
    "Refs/" + fastarefbase + ".pk"
  script:
    "scripts/parse_reference_db.py"
    
checkpoint write_single_fastas:
  input:
    "Refs/" + fastarefbase + ".pk" 
  output: 
    directory("Refs/Single_Fastas"), 
    os.path.join(config['OutFolders']['configs_fld'],config['index_cfg'])
  params:
    pymod=pymod_path
  script:
    "scripts/write_single_fastas.py"

def aggregate_index_single_fastas(wildcards):
  chk_output = checkpoints.write_single_fastas.get(**wildcards).output[0]
  file_names = expand("Refs/BWT_RefIndexes/{ref}.{ext}.bt2",
                      ref=glob_wildcards(os.path.join(chk_output, "{ref}.fa")).ref,
                      ext=["1","2","3","4","rev.1","rev.2"])
  return file_names

def aggregate_kmer_filter_run(wildcards):
  chk_output = checkpoints.kmer_filter_run.get(**wildcards).output
  file_names = expand("FQ_splited/{sample}_R{ind}_splited.fastq",
                      sample = glob_wildcards(os.path.join(chk_output[2], "{sample}.fastq")).sample,
                      ind = ["1", "2"]) + \
               expand("FQ_splited/{sample}_filtered.fastq",
                      sample = glob_wildcards(os.path.join(chk_output[1], "{sample}_filtered.fastq")).sample,
                      ind = ["1", "2"])
  return file_names

rule index_single_fastas:
  input:
    "Refs/Single_Fastas/{ref}.fa"
  output:
    multiext("Refs/BWT_RefIndexes/{ref}.", 
             "1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2")
  shell:
    "bowtie2-build {input} Refs/BWT_RefIndexes/{wildcards.ref}"

rule bowtie2_prep:
  input:
    "readsConfig_filt.pk",
    os.path.join(config['OutFolders']['configs_fld'],config['index_cfg']),
    aggregate_index_single_fastas,
  params:
    pymod=pymod_path
  output:
    "bowtie2_commands.sh",
    "Results/align.pk",
  script:
    "scripts/bowtie2_prep.py"

rule samtools_index_sorted:
  input:
    "{file}_SORTED.bam"
  output:
    "{file}_SORTED.bam.bai"
  shell:
    "samtools index {input}"

rule samtools_index_reduced:
  input:
    "{file}_reduced.bam"
  output:
    "{file}_reduced.bam.bai"
  shell:
    "samtools index {input}"

rule samtools_stats_sorted:
  input:
    multiext("{file}_SORTED", ".bam", ".bam.bai")
  output:
    "{file}_samtools_stats",
  shell:
    "samtools stats {input[0]} | grep ^SN | cut -f2-3 > {output}"

rule samtools_stats:
  input:
    multiext("{file}_reduced", ".bam", ".bam.bai")
  output:
    "{file}_samtools_stats_reduced",
  shell:
    "samtools stats {input[0]} | grep ^SN | cut -f2-3 > {output}"

rule filter_bam:
  input:
    "{file}_SORTED.bam"
  output:
    "{file}_reduced.bam",
    "{file}_reduced_OTHER.bam",
  script:
    "scripts/filter_bam.py"

rule samtools_depth_reduced:
  input:
    "{file}_reduced.bam",
    "{file}_reduced.bam.bai"
  output:
    "{file}_reduced.bam.depth"
  shell:
    "samtools depth {input[0]} | cut -f2- > {output}"

rule samtools_depth_sorted:
  input:
    "{file}_SORTED.bam",
    "{file}_SORTED.bam.bai"
  output:
    "{file}_SORTED.bam.depth"
  shell:
    "samtools depth {input[0]} | cut -f2- > {output}"

def aggregate_bowtie2_run(wildcards):
  chk_output = checkpoints.bowtie2_run.get(**wildcards).output[0]
  file_names = expand(
    "BWT_Results/{pat}{ext}",
    pat = glob_wildcards(os.path.join(chk_output, "{pat}_SORTED.bam")).pat,
    ext=[
      "_reduced.bam.depth", "_SORTED.bam.depth", 
      "_samtools_stats", "_samtools_stats_reduced"
    ])
  return file_names

def aggregate_bowtie2_map(wildcards):
  chk_output = checkpoints.bowtie2_run.get(**wildcards).output[0]
  file_names = expand(
    "BWT_Results/{pat}{ext}",
    pat = glob_wildcards(os.path.join(chk_output, "{pat}_SORTED.bam")).pat,
    ext=[
      ".bam"
    ])
  return file_names

rule alignments_stats:
  input:
    "Results/align.pk",
    "Refs/" + fastarefbase + ".pk",
    aggregate_bowtie2_run,
  output:
    "Results/align_stats.pk",
    "Results/reads_nbr.pk",
  params:
    pymod = pymod_path,
  script:
    "scripts/alignments_stats.py"

rule write_results:
  input:
    "Results/align_stats.pk",
    "Results/reads_nbr.pk",
    "Refs/" + fastarefbase + ".pk",
    readsconfig
  output:
    "Results/xls_samples_index.txt",
    "Results/putative_alleles.txt",
    "Results/putative_groups.txt",
    "Results/samtools_stats.yaml",
    "configs/asm_config.yaml"
  params:
    pymod = pymod_path,
  script:
    "scripts/write_results.py"

checkpoint bowtie2_run:
  input:
    "bowtie2_commands.sh",
    aggregate_kmer_filter_run
  output:
    directory("BWT_Results")
  threads: 
    workflow.cores * 0.75
  shell:
    "parallel -j {threads} -a {input[0]}"

### haploasm

rule haploasm:
  input:
    outdir + "contigs_closest_alleles.txt",
    outdir + "All_samples_phylogeny.pdf",

rule assembly_prep:
  input:
    "Refs/" + fastarefbase + ".pk",
    readsconfig,
  output:
    outdir + "assembly_commands.sh",
    outdir + "asm_contigs.pk",
    outdir + "refs/" + fastarefbase + ".pk.fasta",
    outdir + "refs/" + fastarefbase + "_phylo.fasta",
  params:
    pymod=pymod_path,
  script:
    "scripts/assembly_prep.py"

rule assembly_run:
  input:
    outdir + "assembly_commands.sh",
    outdir + "refs/" + fastarefbase + ".pk.fasta",
    aggregate_kmer_filter_run
  output:
    expand(outdir + "{sample}/contigs.fa", sample=samp) 
  params:
    jobs=int( workflow.cores / config["threaded"])
  shell:
    "parallel -j {params.jobs} -a {input[0]}"

rule clean_contigs:
  input:
    outdir + "asm_contigs.pk",
    expand(outdir + "{sample}/contigs.fa", sample=samp) 
  output:
    outdir + "asm_contigs_clean.pk",
    expand(outdir + "{sample}/cleaned_contigs.fa", sample=samp) 
  script:
    "scripts/clean_contigs.py"

rule do_yass:
  input:
    outdir + "{sample}/cleaned_contigs.fa",
    outdir + "refs/" + fastarefbase + ".pk.fasta",
  output:
    outdir + "{sample}/yass_stdout"
  params:
    yass_params=config["YASSstep"]["yass_params"]
  shell:
    "yass -d 2 {params.yass_params} {input[0]} {input[1]} > {output[0]}"

rule parse_yass:
  input:
    outdir + "{sample}/yass_stdout",
    outdir + "{sample}/cleaned_contigs.fa",
    outdir + "asm_contigs.pk",
    "Refs/" + fastarefbase + ".pk",
  output:
    outdir + "{sample}/yass_oriented_contigs.fa",
    outdir + "{sample}/phylo_truncated_contigs.fa",
  params:
    pymod=pymod_path,
  script:
    "scripts/parse_yass.py"

rule get_phylo_fasta:
  input:
    expand(outdir + "{sample}/phylo_truncated_contigs.fa", sample=samp),
    outdir + "refs/" + fastarefbase + "_phylo.fasta",
  output:
    outdir + "phylo_all_samples.fa"
  shell:
    "cat {input} > {output}"
    
rule get_phylo_contigs:
  input:
    outdir + "{sample}/phylo_truncated_contigs.fa",
    outdir + "refs/" + fastarefbase + "_phylo.fasta",
  output:
    outdir + "{sample}/phylo_all_contigs.fa"
  shell:
    "cat {input} > {output}"

rule do_muscle:
  input:
    "{prefix}/phylo_all_{postfix}.fa",
  output:
    "{prefix}/MUSCLE_phylo_all_{postfix}.fa",
  shell:
    "muscle -in {input} -out {output}"

rule fasta2phylip:
  input:
    "{prefix}/MUSCLE_{postfix}.fa",
  output:
    "{prefix}/PHYLIP_{postfix}.phylip",
  script:
    "scripts/fasta2phylip.py"

rule do_phyml:
  input:
    "{prefix}/PHYLIP_{postfix}.phylip",
  output:
    "{prefix}/PHYLIP_{postfix}.phylip_phyml_stats.txt",
    "{prefix}/PHYLIP_{postfix}.phylip_phyml_tree.txt",
  shell:
    "phyml --quiet --no_memory_check -i {input} -d nt"

rule draw_trees_samples:
  input:
    outdir + "{sample}/PHYLIP_phylo_all_contigs.phylip_phyml_tree.txt",
  output:
    outdir + "{sample}/{sample}_phylogeny.pdf",
  script:
    "scripts/draw_trees.py"

rule draw_trees_all:
  input:
    outdir + "PHYLIP_phylo_all_samples.phylip_phyml_tree.txt",
  output:
    outdir + "All_samples_phylogeny.pdf",
  script:
    "scripts/draw_trees.py"

rule list_nearest_nodes:
  input:
    outdir + "{sample}/PHYLIP_phylo_all_contigs.phylip_phyml_tree.txt",
  output:
    outdir + "{sample}/contigs_closest_alleles.txt",
  script:
    "scripts/list_nearest_nodes.py"

rule join_node_lists:
  input:
    expand(outdir + "{sample}/contigs_closest_alleles.txt", sample=samp)
  output:
    outdir + "contigs_closest_alleles.txt"
  shell:
    "echo $'#Contigs closest Alleles\n#\tContig\tAllele\tdistance\n' > {output}; cat {input} >> {output}"

