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


They can be installed through bioconductor. For example to install the package 'OmicCircos' open R and type the following

::

  source("http://www.bioconductor.org/biocLite.R")
　biocLite("OmicCircos")


2. HBP depends on the following software pacakges which should be installed and included in the system PATH prior to using HBP.

 a) Bowtie2     (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
 b) samtools(>0.1.19)   (http://samtools.sourceforge.net/)
 c) HiC-Pro     (http://github.com/nservant/HiC-Pro)
 d) Python(>2.7) with *pysam*, *bx*, *numpy*, and *scipy* libraries

3. Once dependencies are installed HBP can be used from the command line using the following command.

::

  Rscript HBP.R --stages 1:4

Features
========

HBP uses a Hi-C raw dataset or an interaction matrix file, a BED format file containing specific sites information, and some tracks files containing histone modification or motif enrichment level information (which can be downloaded from the UCSC or made by the user’s own files) as input files. Users can get a network graph, and get a list of name, degree, closeness, betweenness, local cluster coefficient, eigenvector centrality, and other information of the network node. HBP will also plot a scatter diagram to point out the location of specific sites. And then make a circos picture with tracks, make clusters according to the tracks or according to the topological structure and do a statistical analysis.


Usage of HBP
==============

::

  Rscript HBP.R [--options]

Example for regular interactions calling
:: 

 Rscript HBP.R --inFastqDir fastq --bedFile dm3_CP190_int.bed --wigFile dm3_BEAF.wig --wigFile2 dm3_CTCF.wig --wigFile3 dm3_suHw.wig --stages 1:4 --genome_db /data2/my/dm3/ --genomeName dm3 --bowtiePath /usr/bin/bowtie2 --bowtieIndex /data2/my/bowtie2_index/dm3 --enzyme DpnII  



Parameters
----------


ALL STAGES
~~~~~~~~~~


``stages``
 stages of the pipeline to execute.  stage can be either a single stage (e.g 1 or a range of stagnes e.g 1:4). default = 1:4

``tmpDir``
 this will contain up to 3X the size of the largest input .sra file
 
``genomeName``
 The name of genome
 
``outputpdf``
 output pdf format, if choose false, it will output jpeg format. Default = FALSE

``chrom``
 the chrom to invertigate. default = all
 
``chrstart``
 the chrom location to start. default = 0

``chrend``
 the chrom location to end. default = 0

``resolution``
 Hi-C heatmap resolution. default=100(kb)


STAGE 1 PARAMETERS
~~~~~~~~~~

``input_dir``
 fastq or fastq.gz input dir

``enzymesite``
 Restriction Enzyme sites

``enzymeoverhangs5``
 restriction enzyme overhangs 5'

``hicpro_config``
 hic_pro config file

``hicpro_path``
 the path of HiC-Pro
 
``iced_normalize``
 when set true, the matrix will be normalized by the iced algorithm. default = TRUE
 
``enzyme``
 Restriction Enzyme Name
 

STAGE 2 PARAMETERS
~~~~~~~~~~

``bedFile``
 the path of the specific sites file( BED format)
 
``bedWindow``
 the window of the peak. default = 0
 
``netthreshold``
 the threshold to identify an interaction. default = 0
 
``netplot``
 draw the network plot,when node number is too much,it may be not working. default = TRUE
 
``NetClusterType``
 the method of topological clustering, can be choose from NULL,multileve,edgeBetweenness,walktrap,labelPropagation. default = multileve
 
``NetVertexSize``
 network vertex size. default = 2
 
``NetVertexChangeSize``
 the parameter to change the node vertex, can be choose from NULL,degree,closeness,betweenness,Local_cluster_coefficient,Eigenvector_centrality. default = degree
 
``NetVertexLableDist``
 distance between label and vertex in the network. default = 0.1
 
``NetVertexColor``
 the color of vertex. default = red
 
``NetVertexLabelCex``
 the size of network vertex label. default = 0.3


STAGE 3 PARAMETERS
~~~~~~~~~~

``wigFile``
 the path of track file

``wigFile2``
 the path of track file2

``wigFile3``
 the path of track file3

``circosHmSize``
 track file heatmap line size. default = 0.1

``circosThreshold``
 the threshold to identify an interaction. default = 0
 
``circosCmpSize``
comparision sites line size. default = 5
 
``circosLineWidth``
 circos line width. default = 0.01
 
``circosLinecolor``
 circos line color, if choose rainbow will use random color. default = rainbow


STAGE 4 PARAMETERS
~~~~~~~~~~

``groupNum``
 the random group number to determine the statistical difference. default = 100

``dist_method``
 the method to calculater the statistical distance, can be chosen from manhattan,euclidean,minkowski,chebyshev,mahalanobis,canberra. default = euclidean

``clust_method``
 the method to make tracks enrichment clusters , can be chosen from average,centroid,median,complete,single,ward.D. default = complete
    
``clust_k``
 the cluster number to make. default = 4

``clust_label``
 add labels at the cluster tree picture or not. default = FALSE
    
``threshold``
 analysis threshold. default = 0
    

Output Files
============

``chr*.matrix``

the matrix file of each chromosome

``chr*_network.csv``

this file contains the network information of each node. such as node location infromation, degree and so on. and if we make topological clusters, there will have the cluter information of each node.

``chr*_cluster.csv``

this file contains the tracks information of each node, and if we make tracks enrichment cluster, there will have the cluter information of each node.

``chr*_*****.jpeg/pdf``

the picture of different function. the details are in the supplement.








