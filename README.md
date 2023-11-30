somatypus_nf
============

*Summary*

This is an implementation of the SNV and Indel variant calling pipeline used by
the [Transmissible Cancer Group](https://www.tcg.vet.cam.ac.uk/) at the University of Cambridge.

*Description*

The somatypus_nf pipeline is derived from Adrian Baez-Ortega's Somatypus pipeline.
Its main component is [Platypus](https://github.com/andyrimmer/Platypus).
Specifically, somatypus_nf uses my [fork](https://github.com/kgori/Platypus) of the Platypus software,
which is able to read VCF files with `csi`-format indices. This is needed to work with genomes
with chromosomes larger than 2^29 bases (~537Mb), such as that of the 
[Tasmanian Devil](https://www.ensembl.org/Sarcophilus_harrisii/Info/Index).

[Nextflow](https://www.nextflow.io/) handles coordinating the steps of the analysis.
This repository contains an example nextflow config, which is the one I use to run the pipeline on
the [Sanger](https://www.sanger.ac.uk/) Farm. It should be possible to get the pipeline to run on
any other platform with some tweaking of this config file, but I can't guarantee this.

*Software Dependencies*
 - python 2.7 (for Platypus)
   - pysam 0.20.0
   - cython 0.29.15
   - enum34 1.1.10
 - htslib
 - samtools
 - bcftools
 - vcftools

 There's a Dockerfile in this repo that will build a container with all the dependencies installed.
 You can download the container from Docker hub, e.g.
     docker pull kgori/somatypus-dev:latest
     singularity pull somatypus_container.sif docker://kgori/somatypus-dev:latest
