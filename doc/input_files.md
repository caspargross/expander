# File input

The application can run with or without target annotation.

## Without target file:

The script screens the current directory for folders with the prefix `Sample_` and runs the pipeline on the folder content

## With target file:

The pipeline parses a FOFN (File of filenames) with the default name `sample_targets.tsv` in the working directory. Different FOFN file can be specfied in `config.yaml` or by adding ` --config sample_targets=mysamplefile.tsv` as snakemake run parameters
## PacBio

For PacBio Hifi data we need the following input files:
- sample.bam
- sample.bam.pbi
- sample.consensusreadset.xml`

## ONT Nanopore

Currently not yet implemented