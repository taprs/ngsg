# `ngsg`: rewrite of `NGSgenotyp2` for Snakemake

## Purpose and motivation

Some genomic regions can be too polymorphic for determining their alleles by mapping short reads to a single reference, e.g. due to insertions of transposable elements. Following routine works better here:

1. Match reads against a reference panel of alleles;
2. Find the best-matching allele, quantify allelic dosage using relative sequencing depth;
3. Optionally, assemble the reads mapped to the closest allele if the sample could present a novel allele.

The great and long-established [NGSgenotyp2 pipeline](https://github.com/mathieu-genete/NGSgenotyp) implements this very logic and proves very useful in determining alleles of the self-incompatibility locus (a.k.a. S-locus) in crucifers. Here, we try to marry it with the more recent advances of bioinformatics tooling, especially Snakemake as a time-saving alternative to homebrew pipeline execution code.

`ngsg` shares a good deal of code with its predecessor but:

- is easier to install and run,
- is faster due to smart job scheduling,
- is easier to maintain due to smaller codebase,
- tries to make use of the newer/updated tools,
- tries to be more intuitive.

As a result, it is more straighforward to apply `ngsg` for genotyping other "difficult" regions beyond the S-locus.

Please refer to the [NGSgenotyp2 README](README_NGSgenotyp.md) for more details on the internals of the pipeline.

## Installation

```bash
git clone https://github.com/taprs/ngsg.git
cd ngsg
conda env create -n ngsg_env -f environment.yaml
```

## Execution

The input data provided to the pipeline are:

1. Sequencing reads in gzipped FASTQ format;
2. CSV-like table with sample names and paths to the reads files;
3. Reference panel of alleles in FASTA format.

(Structure of input files is described in detail in the [NGSgenotyp2 README](README_NGSgenotyp.md).)

Having the inputs set up, the easiest way to run the pipeline is as follows:

```bash
conda activate ngsg_env
mkdir ngsg_output
cd ngsg_output
snakemake --cores <n_cores> -s /path/to/ngsg/Snakefile -R haploasm --config readsinfo=/path/to/readsinfo fastaref=/path/to/fastaref
```

Remove `-R haploasm` if only mapping and genotyping (as in `NGSgenotyp genotyp`) is required. 

A more involved way is to copy and modify `config.yaml` and replace the `--config` part with `--configfile /path/to/new_config.yaml`. 
