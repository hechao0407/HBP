HBP
=====

Hi-C Analysis Software


Citation
========

HBP: an optimized and flexible pipeline for the interaction analysis with the Hi-C datasets
Chao He.  Accepted. 


HBP Installation
==================

1. HBP depends on the following R packages.

 a) optparse
 b) grid
 c) lattice
 d) IDPmisc
 e) OmicCircos
 f) stringr
 g) ggplot2
 h) igraph
 i) reshape2
 j) pgirmess
 k) coin
 l) multcomp
 m) flexclust
 n) HiTC
 o) rtracklayer
 p) gplots


They can be installed through bioconductor. For example to install the package 'OmicCircos' open R and type the following

::

  source("http://www.bioconductor.org/biocLite.R")
　biocLite(c("optparse","IDPmisc","OmicCircos","stringr","ggplot2","igraph","reshape2","pgirmess","coin","multcomp","flexclust","HiTC","rtracklayer","gplots"))


2. HBP depends on the following software pacakges which should be installed and included in the system PATH prior to using HBP.

 a) Bowtie2     (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
 b) samtools(>0.1.19)   (http://samtools.sourceforge.net/)
 c) HiC-Pro     (http://github.com/nservant/HiC-Pro)
 d) Python(>2.7) with *pysam*, *bx*, *numpy*, and *scipy* libraries

3. Once dependencies are installed HBP can be installed from the command line using the following command.

::

  R CMD INSTALL HBP_0.1.0.tar.gz

Features
========

The implemented pipeline HBP (Hi-C BED file analysis Pipeline) integrates existing pipelines focusing on individual steps of Hi-C data processing into an all-in-one package with adjustable parameters to infer the consensus 3D structure of genome from raw Hi-C sequencing data. What’s more, HBP could assign statistical confidence estimation for chromatin interactions, and clustering interaction loci according to enrichment tracks or topological structure automatically.


Usage of HBP
==============


Example for regular interactions calling
:: 

 generate_enzyme_file(enzyme="DpnII",enzymesite="GATC",chrom_file="chrom_dm3.sizes",enzymedir="annotation",enzymeoverhangs5=0,genomeName="dm3",resolution=1)
 
 run_hicpro(hicpro_path="HiC-Pro",inputfile="rawdata",configfile="config-hicpro.txt",outdir="dm3")
 
 generate_matrix(all_hic_file="SRR389764_1000_iced.matrix",all_bed_file="SRR389764_1000_abs.bed",outputpdf="FALSE",matrix_dir="dm3",resolution=1,chrom_file="chrom_dm3.sizes")
 
 if_distribution_analysis(all_hic_file="SRR1658648_50000_iced.matrix",all_bed_file="SRR1658648_50000_abs.bed",bedFile="L1_human_only.bed",matrix_dir="L1_GM12878",resolution=50,chrom_file="chrom_hg19.sizes")
 
 network_analysis(bedFile="CTCF_hg19_encodeCluster_GM12878.bed",matrix_dir="CTCF_GM12878",chrom="chr1",chrstart=145000000,chrend=150000000,resolution=50)
 
 circos_plot(bedFile="CTCF_hg19_encodeCluster_GM12878.bed",wig_dir="wig",matrix_dir="CTCF_GM12878",chrom="chr8",chrstart=69300000,chrend=74300000,resolution=50)
 
 statistical_analysis(bedFile="CTCF_hg19_encodeCluster_GM12878.bed",wig_dir="wig",matrix_dir="CTCF_GM12878",chrom="chr8",chrstart=69300000,chrend=74300000,resolution=50)
 
 
Some parameters were set by default, users can learn them more in the help file