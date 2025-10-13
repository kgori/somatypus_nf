somatypus_nf
============

## Summary

This is an implementation of the SNV and Indel variant calling pipeline used by
the [Transmissible Cancer Group](https://www.tcg.vet.cam.ac.uk/) at the University of Cambridge.

## TL;DR

Get the pipeline:

    git clone https://github.com/kgori/somatypus_nf.git
    
Get a container to handle the dependencies (Docker or Singularity):

    docker pull kgori/somatypus-dev:latest
    singularity pull somatypus_container.sif docker://kgori/somatypus-dev:latest

Run:

    nextflow run somatypus_nf/main.nf \
        --inputDir  <PATH_TO_BAM/CRAM_FILES> \
        --outputDir <PATH_TO_PUT_OUTPUT_FILES_IN> \
        --reference <FASTA_REFERENCE>

## Description

The somatypus_nf pipeline is derived from Adrian Baez-Ortega's [Somatypus](https://github.com/baezortega/somatypus) pipeline.
Its main component is [Platypus](https://github.com/andyrimmer/Platypus).
Specifically, somatypus_nf uses my [fork](https://github.com/kgori/Platypus) of the Platypus software,
which is able to read VCF files with `csi`-format indices. This is needed to work with genomes
with chromosomes larger than 2^29 bases (~537Mb), such as that of the 
[Tasmanian Devil](https://www.ensembl.org/Sarcophilus_harrisii/Info/Index).

[Nextflow](https://www.nextflow.io/) handles coordinating the steps of the analysis.
There are two example nextflow config files included: `nextflow.config.sanger`, which is the one I
use to run the pipeline on the [Sanger](https://www.sanger.ac.uk/) Farm, and `nextflow.config`, which
is more generic. It should be possible to get the pipeline to run on any other platform by tweaking 
this config file.

## Software Dependencies

 - python 2.7 (for Platypus)
   - pysam 0.20.0
   - cython 0.29.15
   - enum34 1.1.10
 - htslib
 - samtools
 - bcftools
 - vcftools

There's a Dockerfile in this repo that will build a container with all the dependencies installed.
You can [download the container](https://hub.docker.com/repository/docker/kgori/somatypus-dev/general) from Docker hub, e.g.

    docker pull kgori/somatypus-dev:latest
    singularity pull somatypus_container.sif docker://kgori/somatypus-dev:latest

## Program options

    nextflow run somatypus_nf/main.nf \
        --inputDir       <PATH_TO_INPUT_BAM/CRAM_FILES> \
        --outputDir      <FOLDER_WHERE_RESULTS_WILL_GO> \
        --reference      <FASTA_REF_WITH_FAI_INDEX> \
        --regions        [OPTIONAL_LIST_OF_REGIONS_TO_WORK_ON] \
        --extra          [OPTIONAL_LIST_OF_PLATYPUS_ARGS] \
        --vcfSplitSize   <DEFAULT=500000> \
        --indelSplitSize <DEFAULT=500000> \
        
- `inputDir` a folder containing the indexed BAM or CRAM files to analyse 
- `outputDir` the output folder, which will be created if it doesn't exist
- `reference` the indexed Fasta reference file. Must be the same as used to align the inputs.
- `regions` restrict the analysis to this list of regions, in CHR:START-END format. This can
  be a comma-separated list, or a text file with one region per line.
- `extra` any extra parameters to pass to Platypus, as a comma-separated list.
- `--vcfSplitSize=N` and --`indelSplitSize=N` parallel genotyping jobs with be launched
  for each N variants. Default is 500000.
