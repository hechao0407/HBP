suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("lattice"))
suppressPackageStartupMessages(library("IDPmisc"))
suppressPackageStartupMessages(library("OmicCircos"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("igraph"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("pgirmess"))
suppressPackageStartupMessages(library("coin"))
suppressPackageStartupMessages(library("multcomp"))
suppressPackageStartupMessages(library("flexclust"))
suppressPackageStartupMessages(library("HiTC"))
suppressPackageStartupMessages(library("rtracklayer"))



print ("Starting HBP analysis tool")
Sys.time()
set.seed(1234)

##################################### read commandline paramters #####################################

# read in parameters
option_list <- list(
  
  #---------- GENERAL PARAMETERS ----------#
  make_option(c("--stages"),  default="2:4",help="stages of the pipeline to execute"),
  make_option(c("--tmpDir"),   default="/data1/tmp",help="this will contain up to 3X the size of the largest input .sra file"),
  make_option(c("--genomeName"),  default="hg19",help=" "),
  make_option(c("--outputpdf"),  default="FALSE",help="output pdf format,if choose false,it will output jpeg format"),
  make_option(c("--chrom"),  default="all",help="the chrom to invertigate"),
  make_option(c("--chrstart"),  default="0",help="the chrom location to start"),
  make_option(c("--chrend"),  default="0",help="the chrom location to end"),
  make_option(c("--resolution"),  default="100" ,help="Hi-C heatmap resulution"),
  
  
  
  #---------- STAGE 1 PARAMETERS ----------#
  make_option(c("--input_dir"), default="SRR1658673",help="input dir"),
  make_option(c("--enzyme"),  default="MboI" ,help="Restriction Enzyme Name"),
  make_option(c("--enzymesite"),  default="GATC" ,help="restriction enzyme sites"),
  make_option(c("--enzymeoverhangs5"),  default="0" ,help="restriction enzyme overhangs 5'"),
  make_option(c("--hicpro_config"),  default="config-hicpro.txt" ,help="hic_pro config file"),
  make_option(c("--hicpro_path"),  default="HiC-Pro" ,help="the path of HiC-Pro"),
  make_option(c("--iced_normalize"),  default="TRUE" ,help="when set true, the matrix will be normalized by the iced algorithm"),
  make_option(c("--chrom_file"),  default="chrom_hg19.sizes" ,help="contains the info of chrom, just like the example"),
  
  
  
  
  
  #---------- STAGE 2 PARAMETERS ----------#
  make_option(c("--bedFile"),  default="dm3_BEAF_int.bed",help="bed file"),
  make_option(c("--bedWindow"),  default="0",help="the window of the peak"),
  make_option(c("--netplot"),  default="TRUE",help="draw the network plot,when node number is too much,it may be not working"),
  make_option(c("--NetClusterType"),  default="multileve",help="can be choose from NULL,multileve,edgeBetweenness,walktrap,labelPropagation"),
  make_option(c("--NetVertexSize"),  default="2",help="network vertex size"),
  make_option(c("--NetVertexChangeSize"),  default="degree",help="can be choose from NULL,degree,closeness,betweenness,Local_cluster_coefficient,Eigenvector_centrality"),
  make_option(c("--NetVertexLableDist"),  default="0.1",help="distance between label and vertex in the network"),
  make_option(c("--NetVertexColor"),  default="red",help="the color of vertex"),
  make_option(c("--NetVertexLabelCex "),  default="0.3",help="the size of network vertex label"),
  
  
  
  
  
  
  #---------- STAGE 3 PARAMETERS ----------#
  make_option(c("--wigFile"),  default="dm3_H3K4me1.wig",help="wig file"),
  make_option(c("--wigFile2"),  default="dm3_H3K4me3_male.wig",help="wig file"),
  make_option(c("--wigFile3"),  default="dm3_H3K9Ac_male.wig",help="wig file"),
  make_option(c("--circosHmSize"),  default="0.1",help="heatmap line size"),
  make_option(c("--circosCmpSize"),  default="5",help="comparision sites line size"),
  make_option(c("--circosLineWidth"),  default="0.01",help="circos line width"),
  make_option(c("--circosLinecolor"),  default="rainbow",help="circos line color"),
  
  
  
  #---------- STAGE 4 PARAMETERS ----------#
  make_option(c("--groupNum"),  default="100",help="analysis group number"),
  make_option(c("--dist_method"),  default="euclidean",help="manhattan,euclidean,minkowski,chebyshev,mahalanobis,canberra"),
  make_option(c("--clust_method"),  default="complete",help="average,centroid,median,complete,single,ward.D,density"),
  make_option(c("--clust_k"),  default="5",help="the cluster number we get"),
  make_option(c("--threshold"),  default="0",help="analysis threshold")
  
  
  
)



# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

#print(opt)

if (grepl( ":",opt$stages))
{
  stages = strsplit(opt$stages,split=":")[[1]]
  if (length(stages) == 2)
  {
    opt$stages = seq(as.numeric(stages[1]),as.numeric(stages[2]))
  }
  if (length(stages) == 1)
  {
    opt$stages = as.numeric(stages[1])
  }
}


# gather arguments
input_dir=as.character(opt["input_dir"])
tmpDir = as.character(opt["tmpDir"])
genomeName = as.character(opt["genomeName"])
enzyme = as.character(opt["enzyme"])
resolution = as.numeric(as.character(opt["resolution"]))
bedFile=as.character(opt["bedFile"])
bedWindow = as.numeric(as.character(opt["bedWindow"]))
NetClusterType=as.character(opt["NetClusterType"])
wigFile=as.character(opt["wigFile"])
wigFile3=as.character(opt["wigFile3"])
wigFile2=as.character(opt["wigFile2"])
netplot=as.character(opt["netplot"])
groupNum = as.numeric(as.character(opt["groupNum"]))
circosHmSize = as.numeric(as.character(opt["circosHmSize"]))
circosLineWidth = as.numeric(as.character(opt["circosLineWidth"]))
circosCmpSize = as.numeric(as.character(opt["circosCmpSize"]))
circosLinecolor=as.character(opt["circosLinecolor"])
outputpdf=as.character(opt["outputpdf"])
chrom=as.character(opt["chrom"])
chrstart = as.numeric(as.character(opt["chrstart"]))
chrend = as.numeric(as.character(opt["chrend"]))
threshold = as.numeric(as.character(opt["threshold"]))
NetVertexSize = as.numeric(as.character(opt["NetVertexSize"]))
NetVertexLabelCex = as.numeric(as.character(opt["NetVertexLabelCex"]))
NetVertexLableDist = as.numeric(as.character(opt["NetVertexLableDist"]))
NetVertexColor=as.character(opt["NetVertexColor"])
NetVertexChangeSize=as.character(opt["NetVertexChangeSize"])
dist_method=as.character(opt["dist_method"])
clust_method=as.character(opt["clust_method"])
clust_k=as.numeric(as.character(opt["clust_k"]))
enzymesite=as.character(opt["enzymesite"])
enzymeoverhangs5=as.numeric(as.character(opt["enzymeoverhangs5"]))
hicpro_config=as.character(opt["hicpro_config"])
hicpro_path=as.character(opt["hicpro_path"])
iced_normalize=as.character(opt["iced_normalize"])
chrom_file=as.character(opt["chrom_file"])







##########################step1############################


if (1 %in% opt$stages)
{
  chrom_info=read.table(chrom_file,fill=TRUE, stringsAsFactors=FALSE)
  all_genomedb=read.csv("all_genome_db.csv",header=FALSE)
  enzymesitesfile=paste(enzyme,"_resfrag_",genomeName,".bed",sep="")
  enzyme_dir=list.files(path="annotation",full.names=F,pattern=enzymesitesfile)
  if(length(enzyme_dir)==0)
  {
    for(i in 1:length(all_genomedb))
    {
      dbfind=-1
      dbrequire=FALSE
      r_genomedb=NULL
      dbfind=regexpr(genomeName,all_genomedb[1,i])
      if(dbfind!=-1)
      {
        r_genomedb=as.character(all_genomedb[1,i])
        dbrequire=require(all_genomedb[1,i],character.only = TRUE)
        break
      }
    }
    res_ok=FALSE
    if(!is.null(r_genomedb))
    {
      if(dbrequire==FALSE)
      {
        print(paste("please install the R package ",r_genomedb,sep=""))
      }
      else
      {
        all_chr <- chrom_info[,1]
        resFrag <- getRestrictionFragmentsPerChromosome(resSite=enzymesite, chromosomes=all_chr, overhangs5=enzymeoverhangs5, genomePack=r_genomedb)
        allRF <- do.call("c",resFrag)
        names(allRF) <- unlist(sapply(resFrag, function(x){paste0("HIC_", seqlevels(x), "_", 1:length(x))}))
        export(allRF, format="bed", con=paste("annotation/",enzymesitesfile,sep=""))
        res_ok=TRUE
      }
    }else
    {
      print(paste("can not find the genome file",sep=""))
    }
  }
  else
  {
    print(paste("enzymesites file ",enzymesitesfile," is already existed, skip to generate it"))
    res_ok=TRUE
  }
  
  
  if(res_ok==TRUE)
  {
    #rrr=as.character(resolution*1000)
    rrr=paste(resolution,"000",sep = "")
    print(rrr)
    all_bed_file=paste(genomeName,"/hic_results/matrix/",input_dir,"/raw/",rrr,"/",input_dir,"_",rrr,"_abs.bed",sep="")
    if((iced_normalize=="TRUE")&&(iced_normalize=="true"))
    {
      all_hic_file=paste(genomeName,"/hic_results/matrix/",input_dir,"/iced/",rrr,"/",input_dir,"_",rrr,"_iced.matrix",sep="")
      hic_path=paste(genomeName,"/hic_results/matrix/",input_dir,"/iced/",rrr,sep="")
      hic_file=paste(genomeName,"_",rrr,"_iced.matrix",sep="")
    }
    else
    {
      all_hic_file=paste(genomeName,"/hic_results/matrix/",input_dir,"/raw/",rrr,"/",input_dir,"_",rrr,".matrix",sep="")
      hic_path=paste(genomeName,"/hic_results/matrix/",input_dir,"/raw/",rrr,sep="")
      hic_file=paste(genomeName,"_",rrr,".matrix",sep="")
    }
    hic_data_dir=list.files(path=hic_path,full.names=F,pattern=hic_file)
    if(length(hic_data_dir)==0)
    {
      hiccmd=paste(hicpro_path," -i rawdata -o ",genomeName," -c ",hicpro_config,sep="")
      print(hiccmd)
      system(hiccmd)
    }else
    {
      print(paste("hic file ",hic_file," is already existed, skip to generate it"))
    }
    all_bed_data=read.table(file=all_bed_file, fill=TRUE, stringsAsFactors=FALSE)
    all_hic_data=read.table(file=all_hic_file,fill=TRUE, stringsAsFactors=FALSE)
    chr_num=dim(chrom_info)[1]
    for(tttt in 1:chr_num)
    {
      chrname=chrom_info[tttt,1]
      chr_bed=NULL
      chr_bed=rbind(chr_bed,all_bed_data[all_bed_data[,1]==chrname,])
      chr_max_length=dim(chr_bed)[1]
      chr_hic=NULL
      chr_hic=all_hic_data[all_hic_data[,1]>=chr_bed[1,4],]
      chr_hic=chr_hic[chr_hic[,1]<=chr_bed[chr_max_length,4],]
      chr_hic=chr_hic[chr_hic[,2]>=chr_bed[1,4],]
      chr_hic=chr_hic[chr_hic[,2]<=chr_bed[chr_max_length,4],]
      chr_hic_data=matrix(data=0, nrow = chr_max_length, ncol = chr_max_length)
      for(i in 1:dim(chr_hic)[1])
      {
        chr_hic_data[(chr_hic[i,1]-chr_bed[1,4]),(chr_hic[i,2]-chr_bed[1,4])]=chr_hic[i,3]
        chr_hic_data[(chr_hic[i,2]-chr_bed[1,4]),(chr_hic[i,1]-chr_bed[1,4])]=chr_hic[i,3]
      }
      tmpfilename=paste(genomeName,"/",chrname,".matrix",sep="")
      write.table(chr_hic_data,tmpfilename,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
      
      if((outputpdf=="TRUE")||(outputpdf=="true"))
      {
        pdf(paste(genomeName,"/",chrname,"_heatmap.pdf",sep=""))
        
      }else
      {
        jpeg(paste(genomeName,"/",chrname,"_heatmap.jpeg",sep=""),width=1000,height=1000,quality = 100)
      }
      hm_dim=dim(chr_hic_data)[1]
      hm_mean=mean(chr_hic_data)
      for(n in 1:hm_dim)
      {
        
        for(nn in 1:hm_dim)
        {
          if(chr_hic_data[n,nn]>10*hm_mean)
            chr_hic_data[n,nn]=10*hm_mean
        }
      }
      m_color=heat.colors(100)
      n_color=m_color
      for(i in 1:100)
      {
        m_color[i]=n_color[101-i]
      }
      tmphm=levelplot(chr_hic_data,col.regions=m_color,xlim=c(0,hm_dim),ylim=c(0,hm_dim),xlab="chrom",ylab="chrom")
      plot(tmphm)
      dev.off()
      
      
      rm(chr_hic_data)
      gc()
      print(paste(chrname," finish!",sep=""))
    }
  }
  
}


#########################step2############################

#find overlap
if (2 %in% opt$stages)
{
  gc()
  source("calbed.R")
  matrix_dir=list.files(path=genomeName,full.names=F,pattern=".matrix")
  matrix_full_dir=list.files(path=genomeName,full.names=T,pattern=".matrix")
  m_bed=load_bed(bedFile)
  chrNum=length(matrix_dir)
  if(chrom=="all")
  {
    for (i in 1:chrNum)
    {
      
      
      tmpNum=regexpr(".matrix",matrix_dir[i])
      chrName=substr(matrix_dir[i],1,tmpNum-1)
      chrBed=choose_chr_bed(m_bed,chrName)
      chrBedNum=dim(chrBed)[1]
      print(matrix_full_dir[i])
      if(chrBedNum>0)
      {
        print(chrName)
        chrCmap=read.table(file=matrix_full_dir[i], fill=TRUE, stringsAsFactors=FALSE)
        chrTotSize=dim(chrCmap)[1]
        chrBedBin=check_bed_bin(chrBed,resolution*1000)
        chrBedMatrix=convert_bed_to_matrix(chrBed,resolution*1000,chrName,chrTotSize,bedWindow)
        chrBedToBedInter=find_bed_to_bed_interaction(chrCmap,chrBedMatrix,chrBedBin,chrBed,chrName,chrCmap)
        
        if((netplot=="true")||(netplot=="TRUE"))
        {
          netgraph=data.frame("p1"=character(dim(chrBedToBedInter)[1]),"p2"=character(dim(chrBedToBedInter)[1]),"weight"=numeric(dim(chrBedToBedInter)[1]))
          netgraph[,1]=as.data.frame(paste(chrBedToBedInter[,2],":",chrBedToBedInter[,3],"-",chrBedToBedInter[,4],sep=""))
          netgraph[,2]=as.data.frame(paste(chrBedToBedInter[,8],":",chrBedToBedInter[,9],"-",chrBedToBedInter[,10],sep=""))
          netgraph[,3]=as.data.frame(chrBedToBedInter[,13])
          
          
          set.seed(1234)
          #pdf(paste(genomeName,"/",chrName,"_netplot.pdf",sep=""))
          
          if((outputpdf=="TRUE")||(outputpdf=="true"))
          {
            pdf(paste(genomeName,"/",chrName,"_netplot.pdf",sep=""),width = 8,height = 8)
            
          }else
          {
            jpeg(paste(genomeName,"/",chrName,"_netplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
          }
          
          g = graph.data.frame(netgraph,directed = F)
          
          
          netnodename=names(V(g))
          netnodechrnum=regexpr(":",netnodename)
          netnodechr=substr(netnodename,1,netnodechrnum-1)
          netnodestartnum=regexpr("-",netnodename)
          netnodestart=as.numeric(substr(netnodename,netnodechrnum+1,netnodestartnum-1))
          netnodeendnum=nchar(netnodename)
          netnodeend=as.numeric(substr(netnodename,netnodestartnum+1,netnodeendnum))
          
          netcsv=data.frame("chrom"=character(0),"start"=numeric(0),"end"=numeric(0),"degree"=numeric(0),"closeness"=numeric(0),"betweenness"=numeric(0),"Local_cluster_coefficient"=numeric(0),"Eigenvector_centrality"=numeric(0),"membership"=numeric(0),stringsAsFactors=FALSE)
          netdegree=degree(g)
          netcloseness=closeness(g)
          netbetweenness=betweenness(g)
          netcoefficient=transitivity(g, type="local")
          netcentrality=evcent(g)$vector
          
          netcsv[1:(length(netdegree)),1]=netnodechr
          netcsv[1:(length(netdegree)),2]=netnodestart
          netcsv[1:(length(netdegree)),3]=netnodeend
          
          netcsv[1:(length(netdegree)),4]=as.data.frame(netdegree)
          netcsv[1:(length(netdegree)),5]=as.data.frame(netcloseness)
          netcsv[1:(length(netdegree)),6]=as.data.frame(netbetweenness)
          netcsv[1:(length(netdegree)),7]=as.data.frame(netcoefficient)
          netcsv[1:(length(netdegree)),8]=as.data.frame(netcentrality)
          
          if(NetVertexChangeSize=="degree")
          {
            V(g)$deg<-netcsv[,4]
            deg_range=range(netcsv[,4])
            V(g)$size=NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
          }else if(NetVertexChangeSize=="closeness")
          {
            V(g)$deg<-netcsv[,5]
            deg_range=range(netcsv[,5])
            V(g)$size=NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
          }else if(NetVertexChangeSize=="betweenness")
          {
            V(g)$deg<-netcsv[,6]
            deg_range=range(netcsv[,6])
            V(g)$size=NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
          }else if(NetVertexChangeSize=="Local_cluster_coefficient")
          {
            V(g)$deg<-netcsv[,7]
            deg_range=range(netcsv[,7])
            V(g)$size=NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
          }else if(NetVertexChangeSize=="Eigenvector_centrality")
          {
            V(g)$deg<-netcsv[,8]
            deg_range=range(netcsv[,8])
            V(g)$size=NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
            V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
          }else
          {
            V(g)$size=NetVertexSize
          }
          
          
          
          
          if(NetClusterType=="NULL")
          {
            plot(g,layout=layout.fruchterman.reingold, vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
            netcsv[1:(length(netdegree)),9]="NULL"
          }
          if(NetClusterType=="edgeBetweenness")
          {
            system.time(ec <- edge.betweenness.community(g)) 
            print(modularity(ec)) 
            netcsv[1:(length(netdegree)),9]=ec$membership
            plot(ec, g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
          }
          if(NetClusterType=="walktrap")
          {
            system.time(wc <- walktrap.community(g)) 
            netcsv[1:(length(netdegree)),9]=wc$membership
            
            print(modularity(wc)) 
            plot(wc , g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex) 
            
          }
          if(NetClusterType=="multileve")
          {
            system.time(mc <- multilevel.community(g, weights=NA)) 
            print(modularity(mc)) 
            plot(mc, g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex) 
            netcsv[1:(length(netdegree)),9]=mc$membership
            
          }
          if(NetClusterType=="labelPropagation")
          {
            system.time(lc <- label.propagation.community(g)) 
            print(modularity(lc)) 
            plot(lc , g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex) 
            netcsv[1:(length(netdegree)),9]=lc$membership
            
          }
          dev.off()
          write.csv(netcsv,paste(genomeName,"/",chrName,"_network.csv",sep=""),row.names = FALSE)
        }
        
        
        bedIplot=cbind(rbind(chrBedToBedInter[,14],chrBedToBedInter[,15]),rbind(chrBedToBedInter[,15],chrBedToBedInter[,14]))
        hm_dim=dim(chrCmap)[1]
        chrCmap=as.matrix(chrCmap)
        
        
        hm_mean=mean(chrCmap)
        for(j in 1:hm_dim)
        {
          
          chrCmap[j,which(chrCmap[j,]>5*hm_mean)]=5*hm_mean
        }
        chrhmCmap=melt(chrCmap)
        print(paste("plot ",chrName,"bed picture",sep=""))
        
        if((outputpdf=="TRUE")||(outputpdf=="true"))
        {
          pdf(paste(genomeName,"/",chrName,"_bedplot.pdf",sep=""),width = 8,height = 8)
          
        }else
        {
          jpeg(paste(genomeName,"/",chrName,"_bedplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
        }
        
        grid.newpage()
        heatmapViewport <- viewport(height=0.5, width=0.5, x=0.25,y=0.5) 
        scatterViewport <- viewport(height=0.5, width=0.5, x=0.75,y=0.5)
        densityViewport <- viewport(height=0.25,width=0.5, x=0.75,y=0.125)
        hmdensityViewport <- viewport(height=0.25,width=0.5,x=0.25,y=0.125)
        jit=position_jitter(width=0.5)
        hmrange=range(chrhmCmap[,1])
        bedIplot=bedIplot+hmrange[1]
        chrhm = ggplot(chrhmCmap, aes(x=Var1, y=Var2, fill=value))+scale_y_discrete(breaks=seq(0, 10, 5))+xlab('chrom')+ylab("chrom")+scale_fill_gradient(low='white', high='red')+geom_tile()+guides(fill=FALSE)
        chrbedplot=qplot(bedIplot[1,],bedIplot[2,],alpha=I(1/10),size=I(1))+xlab('chrom')+ylab("chrom")+geom_jitter(position=jit,colour="black",alpha=1/100)
        
        chrbeddensitydata=as.data.frame(c(chrBedToBedInter[,14],chrBedToBedInter[,15]))
        colnames(chrbeddensitydata)="bed"
        chrbeddensity=ggplot(chrbeddensitydata)+geom_density(aes(x=bed))
        chrhmdensitydata=NULL
        for(iiii in 1:chrTotSize)
        {
          #pp[iiii,1]=length(which(chrCmap[,iiii]>0))
          chrhmdensitydata=c(chrhmdensitydata,which(chrCmap[,iiii]>0))
        }
        chrhmdensitydata=as.data.frame(chrhmdensitydata)
        colnames(chrhmdensitydata)="chrom"
        chrcmapdensity=ggplot(chrhmdensitydata)+geom_density(aes(x=chrom))
        print(chrhm,vp=heatmapViewport)
        print(chrbedplot,vp=scatterViewport)
        print(chrbeddensity,vp=densityViewport)
        print(chrcmapdensity,vp=hmdensityViewport)
        
        dev.off()
        
        write.table(chrBedToBedInter,file=paste(genomeName,"/",chrName,"_BedToBedInter.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
        rm(chrCmap)
        rm(chrBed)
        rm(chrBedBin)
        rm(chrBedMatrix)
        rm(chrBedToBedInter)
        gc()
      }
    }
  }else
  {
    for (i in 1:chrNum)
    {
      
      
      tmpNum=regexpr(".matrix",matrix_dir[i])
      chrName=substr(matrix_dir[i],1,tmpNum-1)
      if(chrName==chrom)
      {
        if(chrend<=0)
        {
          chrBed=choose_chr_bed(m_bed,chrName)
          chrBedNum=dim(chrBed)[1]
          print(matrix_full_dir[i])
          if(chrBedNum>0)
          {
            print(chrName)
            chrCmap=read.table(file=matrix_full_dir[i], fill=TRUE, stringsAsFactors=FALSE)
            chrTotSize=dim(chrCmap)[1]
            chrBedBin=check_bed_bin(chrBed,resolution*1000)
            chrBedMatrix=convert_bed_to_matrix(chrBed,resolution*1000,chrName,chrTotSize,bedWindow)
            chrBedToBedInter=find_bed_to_bed_interaction(chrCmap,chrBedMatrix,chrBedBin,chrBed,chrName,chrCmap)
            
            if((netplot=="true")||(netplot=="TRUE"))
            {
              netgraph=data.frame("p1"=character(dim(chrBedToBedInter)[1]),"p2"=character(dim(chrBedToBedInter)[1]),"weight"=numeric(dim(chrBedToBedInter)[1]))
              netgraph[,1]=as.data.frame(paste(chrBedToBedInter[,2],":",chrBedToBedInter[,3],"-",chrBedToBedInter[,4],sep=""))
              netgraph[,2]=as.data.frame(paste(chrBedToBedInter[,8],":",chrBedToBedInter[,9],"-",chrBedToBedInter[,10],sep=""))
              netgraph[,3]=as.data.frame(chrBedToBedInter[,13])
              
              set.seed(1234)
              if((outputpdf=="TRUE")||(outputpdf=="true"))
              {
                pdf(paste(genomeName,"/",chrName,"_netplot.pdf",sep=""),width = 8,height = 8)
                
              }else
              {
                jpeg(paste(genomeName,"/",chrName,"_netplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
              }
              g = graph.data.frame(netgraph,directed = F)
              
              netnodename=names(V(g))
              netnodechrnum=regexpr(":",netnodename)
              netnodechr=substr(netnodename,1,netnodechrnum-1)
              netnodestartnum=regexpr("-",netnodename)
              netnodestart=as.numeric(substr(netnodename,netnodechrnum+1,netnodestartnum-1))
              netnodeendnum=nchar(netnodename)
              netnodeend=as.numeric(substr(netnodename,netnodestartnum+1,netnodeendnum))
              
              netcsv=data.frame("chrom"=character(0),"start"=numeric(0),"end"=numeric(0),"degree"=numeric(0),"closeness"=numeric(0),"betweenness"=numeric(0),"Local_cluster_coefficient"=numeric(0),"Eigenvector_centrality"=numeric(0),"membership"=numeric(0),stringsAsFactors=FALSE)
              netdegree=degree(g)
              netcloseness=closeness(g)
              netbetweenness=betweenness(g)
              netcoefficient=transitivity(g, type="local")
              netcentrality=evcent(g)$vector
              
              netcsv[1:(length(netdegree)),1]=netnodechr
              netcsv[1:(length(netdegree)),2]=netnodestart
              netcsv[1:(length(netdegree)),3]=netnodeend
              
              netcsv[1:(length(netdegree)),4]=as.data.frame(netdegree)
              netcsv[1:(length(netdegree)),5]=as.data.frame(netcloseness)
              netcsv[1:(length(netdegree)),6]=as.data.frame(netbetweenness)
              netcsv[1:(length(netdegree)),7]=as.data.frame(netcoefficient)
              netcsv[1:(length(netdegree)),8]=as.data.frame(netcentrality)
              
              if(NetVertexChangeSize=="degree")
              {
                V(g)$deg<-netcsv[,4]
                deg_range=range(netcsv[,4])
                V(g)$size=NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
              }else if(NetVertexChangeSize=="closeness")
              {
                V(g)$deg<-netcsv[,5]
                deg_range=range(netcsv[,5])
                V(g)$size=NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
              }else if(NetVertexChangeSize=="betweenness")
              {
                V(g)$deg<-netcsv[,6]
                deg_range=range(netcsv[,6])
                V(g)$size=NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
              }else if(NetVertexChangeSize=="Local_cluster_coefficient")
              {
                V(g)$deg<-netcsv[,7]
                deg_range=range(netcsv[,7])
                V(g)$size=NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
              }else if(NetVertexChangeSize=="Eigenvector_centrality")
              {
                V(g)$deg<-netcsv[,8]
                deg_range=range(netcsv[,8])
                V(g)$size=NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
              }else
              {
                V(g)$size=NetVertexSize
              }
              
              
              
              
              if(NetClusterType=="NULL")
              {
                plot(g,layout=layout.fruchterman.reingold, vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                netcsv[1:(length(netdegree)),9]="NULL"
              }
              if(NetClusterType=="edgeBetweenness")
              {
                system.time(ec <- edge.betweenness.community(g)) 
                print(modularity(ec)) 
                netcsv[1:(length(netdegree)),9]=ec$membership
                plot(ec, g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              }
              if(NetClusterType=="walktrap")
              {
                system.time(wc <- walktrap.community(g)) 
                netcsv[1:(length(netdegree)),9]=wc$membership
                
                print(modularity(wc)) 
                plot(wc , g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex) 
                
              }
              if(NetClusterType=="multileve")
              {
                system.time(mc <- multilevel.community(g, weights=NA)) 
                print(modularity(mc)) 
                plot(mc, g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex) 
                netcsv[1:(length(netdegree)),9]=mc$membership
                
              }
              if(NetClusterType=="labelPropagation")
              {
                system.time(lc <- label.propagation.community(g)) 
                print(modularity(lc)) 
                plot(lc , g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex) 
                netcsv[1:(length(netdegree)),9]=lc$membership
                
              }
              dev.off()
              write.csv(netcsv,paste(genomeName,"/",chrName,"_network.csv",sep=""),row.names = FALSE)
              
            }
            bedIplot=cbind(rbind(chrBedToBedInter[,14],chrBedToBedInter[,15]),rbind(chrBedToBedInter[,15],chrBedToBedInter[,14]))
            
            hm_dim=dim(chrCmap)[1]
            chrCmap=as.matrix(chrCmap)
            hm_mean=mean(chrCmap)
            for(j in 1:hm_dim)
            {
              
              chrCmap[j,which(chrCmap[j,]>5*hm_mean)]=5*hm_mean
            }
            chrhmCmap=melt(chrCmap)
            print(paste("plot ",chrName,"bed picture",sep=""))
            if((outputpdf=="TRUE")||(outputpdf=="true"))
            {
              pdf(paste(genomeName,"/",chrName,"_bedplot.pdf",sep=""),width = 8,height = 8)
              
            }else
            {
              jpeg(paste(genomeName,"/",chrName,"_bedplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
            }
            
            grid.newpage()
            heatmapViewport <- viewport(height=0.5, width=0.5, x=0.25,y=0.5) 
            scatterViewport <- viewport(height=0.5, width=0.5, x=0.75,y=0.5)
            densityViewport <- viewport(height=0.25,width=0.5, x=0.75,y=0.125)
            hmdensityViewport <- viewport(height=0.25,width=0.5,x=0.25,y=0.125)
            jit=position_jitter(width=0.5)
            hmrange=range(chrhmCmap[,1])
            bedIplot=bedIplot+hmrange[1]
            chrhm = ggplot(chrhmCmap, aes(x=Var1, y=Var2, fill=value))+scale_y_discrete(breaks=seq(0, 10, 5))+xlab('chrom')+ylab("chrom")+scale_fill_gradient(low='white', high='red')+geom_tile()+guides(fill=FALSE)
            chrbedplot=qplot(bedIplot[1,],bedIplot[2,],alpha=I(1/10),size=I(1))+xlab('chrom')+ylab("chrom")+geom_jitter(position=jit,colour="black",alpha=1/100)
            
            chrbeddensitydata=as.data.frame(c(chrBedToBedInter[,14],chrBedToBedInter[,15]))
            colnames(chrbeddensitydata)="bed"
            chrbeddensity=ggplot(chrbeddensitydata)+geom_density(aes(x=bed))
            
            chrhmdensitydata=NULL
            for(iiii in 1:chrTotSize)
            {
              #pp[iiii,1]=length(which(chrCmap[,iiii]>0))
              chrhmdensitydata=c(chrhmdensitydata,which(chrCmap[,iiii]>0))
            }
            chrhmdensitydata=as.data.frame(chrhmdensitydata)
            colnames(chrhmdensitydata)="chrom"
            
            
            chrcmapdensity=ggplot(chrhmdensitydata)+geom_density(aes(x=chrom))
            
            
            
            print(chrhm,vp=heatmapViewport)
            print(chrbedplot,vp=scatterViewport)
            print(chrbeddensity,vp=densityViewport)
            print(chrcmapdensity,vp=hmdensityViewport)
            
            dev.off()
            write.table(chrBedToBedInter,file=paste(genomeName,"/",chrName,"_BedToBedInter.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
            rm(chrCmap)
            rm(chrBed)
            rm(chrBedBin)
            rm(chrBedMatrix)
            rm(chrBedToBedInter)
            gc()
          }
          
          
          
        }else
        {
          chrBed=choose_chr_bed(m_bed,chrName)
          chrBed=chrBed[which(chrBed[,3]<chrend),]
          chrBed[,2]=chrBed[,2]-chrstart
          chrBed[,3]=chrBed[,3]-chrstart
          chrBed=chrBed[which(chrBed[,2]>0),]
          chrBedNum=dim(chrBed)[1]
          
          print(matrix_full_dir[i])
          if(chrBedNum>0)
          {
            
            print(chrName)
            chrCmap=read.table(file=matrix_full_dir[i], fill=TRUE, stringsAsFactors=FALSE)
            
            tmpCmapStart=abs(ceiling((chrstart/(resolution*1000))))
            tmpCmapEnd=abs(ceiling((chrend/(resolution*1000))))
            
            chrTotSize=dim(chrCmap)[1]
            if(tmpCmapEnd<tmpCmapStart)
            {
              print("please input correct start and end number")
              break
            }
            if(tmpCmapEnd>chrTotSize)
            {
              tmpCmapEnd=chrTotSize
            }
            if(tmpCmapStart>chrTotSize)
            {
              tmpCmapStart=chrTotSize
            }
            chrCmap=chrCmap[tmpCmapStart:tmpCmapEnd,tmpCmapStart:tmpCmapEnd]
            chrTotSize=dim(chrCmap)[1]
            
            
            
            chrBedBin=check_bed_bin(chrBed,resolution*1000)
            chrBedMatrix=convert_bed_to_matrix(chrBed,resolution*1000,chrName,chrTotSize,bedWindow)
            
            chrBedToBedInter=find_bed_to_bed_interaction(chrCmap,chrBedMatrix,chrBedBin,chrBed,chrName,chrCmap)
            
            if((netplot=="true")||(netplot=="TRUE"))
            {
              netgraph=data.frame("p1"=character(dim(chrBedToBedInter)[1]),"p2"=character(dim(chrBedToBedInter)[1]),"weight"=numeric(dim(chrBedToBedInter)[1]))
              netgraph[,1]=as.data.frame(paste(chrBedToBedInter[,2],":",(chrBedToBedInter[,3]+chrstart),"-",(chrBedToBedInter[,4]+chrstart),sep=""))
              netgraph[,2]=as.data.frame(paste(chrBedToBedInter[,8],":",(chrBedToBedInter[,9]+chrstart),"-",(chrBedToBedInter[,10]+chrstart),sep=""))
              netgraph[,3]=as.data.frame(chrBedToBedInter[,13])
              
              set.seed(1234)
              
              if((outputpdf=="TRUE")||(outputpdf=="true"))
              {
                pdf(paste(genomeName,"/",chrName,"_netplot.pdf",sep=""),width = 8,height = 8)
                
              }else
              {
                jpeg(paste(genomeName,"/",chrName,"_netplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
              }
              
              g = graph.data.frame(netgraph,directed = F)
              
              netnodename=names(V(g))
              netnodechrnum=regexpr(":",netnodename)
              netnodechr=substr(netnodename,1,netnodechrnum-1)
              netnodestartnum=regexpr("-",netnodename)
              netnodestart=as.numeric(substr(netnodename,netnodechrnum+1,netnodestartnum-1))
              netnodeendnum=nchar(netnodename)
              netnodeend=as.numeric(substr(netnodename,netnodestartnum+1,netnodeendnum))
              
              netcsv=data.frame("chrom"=character(0),"start"=numeric(0),"end"=numeric(0),"degree"=numeric(0),"closeness"=numeric(0),"betweenness"=numeric(0),"Local_cluster_coefficient"=numeric(0),"Eigenvector_centrality"=numeric(0),"membership"=numeric(0),stringsAsFactors=FALSE)
              netdegree=degree(g)
              netcloseness=closeness(g)
              netbetweenness=betweenness(g)
              netcoefficient=transitivity(g, type="local")
              netcentrality=evcent(g)$vector
              
              netcsv[1:(length(netdegree)),1]=netnodechr
              netcsv[1:(length(netdegree)),2]=netnodestart
              netcsv[1:(length(netdegree)),3]=netnodeend
              
              netcsv[1:(length(netdegree)),4]=as.data.frame(netdegree)
              netcsv[1:(length(netdegree)),5]=as.data.frame(netcloseness)
              netcsv[1:(length(netdegree)),6]=as.data.frame(netbetweenness)
              netcsv[1:(length(netdegree)),7]=as.data.frame(netcoefficient)
              netcsv[1:(length(netdegree)),8]=as.data.frame(netcentrality)
              
              if(NetVertexChangeSize=="degree")
              {
                V(g)$deg<-netcsv[,4]
                deg_range=range(netcsv[,4])
                V(g)$size=NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
              }else if(NetVertexChangeSize=="closeness")
              {
                V(g)$deg<-netcsv[,5]
                deg_range=range(netcsv[,5])
                V(g)$size=NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
              }else if(NetVertexChangeSize=="betweenness")
              {
                V(g)$deg<-netcsv[,6]
                deg_range=range(netcsv[,6])
                V(g)$size=NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
              }else if(NetVertexChangeSize=="Local_cluster_coefficient")
              {
                V(g)$deg<-netcsv[,7]
                deg_range=range(netcsv[,7])
                V(g)$size=NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
              }else if(NetVertexChangeSize=="Eigenvector_centrality")
              {
                V(g)$deg<-netcsv[,8]
                deg_range=range(netcsv[,8])
                V(g)$size=NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])/5)]$size=2*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*2/5)]$size=3*NetVertexSize
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*3/5)]$size=4*NetVertexSize 
                V(g)[deg>=(deg_range[1]+(deg_range[2]-deg_range[1])*4/5)]$size=5*NetVertexSize
              }else
              {
                V(g)$size=NetVertexSize
              }
              
              
              
              
              if(NetClusterType=="NULL")
              {
                plot(g,layout=layout.fruchterman.reingold, vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                netcsv[1:(length(netdegree)),9]="NULL"
              }
              if(NetClusterType=="edgeBetweenness")
              {
                system.time(ec <- edge.betweenness.community(g)) 
                print(modularity(ec)) 
                netcsv[1:(length(netdegree)),9]=ec$membership
                plot(ec, g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              }
              if(NetClusterType=="walktrap")
              {
                system.time(wc <- walktrap.community(g)) 
                netcsv[1:(length(netdegree)),9]=wc$membership
                
                print(modularity(wc)) 
                plot(wc , g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex) 
                
              }
              if(NetClusterType=="multileve")
              {
                system.time(mc <- multilevel.community(g, weights=NA)) 
                print(modularity(mc)) 
                plot(mc, g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex) 
                netcsv[1:(length(netdegree)),9]=mc$membership
                
              }
              if(NetClusterType=="labelPropagation")
              {
                system.time(lc <- label.propagation.community(g)) 
                print(modularity(lc)) 
                plot(lc , g,vertex.size=NetVertexSize,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex) 
                netcsv[1:(length(netdegree)),9]=lc$membership
                
              }
              dev.off()
              write.csv(netcsv,paste(genomeName,"/",chrName,"_network.csv",sep=""),row.names = FALSE)
              
            }
            
            
            bedIplot=cbind(rbind(chrBedToBedInter[,14],chrBedToBedInter[,15]),rbind(chrBedToBedInter[,15],chrBedToBedInter[,14]))
            
            hm_dim=dim(chrCmap)[1]
            chrCmap=as.matrix(chrCmap)
            hm_mean=mean(chrCmap)
            for(j in 1:hm_dim)
            {
              
              chrCmap[j,which(chrCmap[j,]>5*hm_mean)]=5*hm_mean
            }
            chrhmCmap=melt(chrCmap)
            print(paste("plot ",chrName,"bed picture",sep=""))
            if((outputpdf=="TRUE")||(outputpdf=="true"))
            {
              pdf(paste(genomeName,"/",chrName,"_bedplot.pdf",sep=""),width = 8,height = 8)
              
            }else
            {
              jpeg(paste(genomeName,"/",chrName,"_bedplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
            }
            
            grid.newpage()
            heatmapViewport <- viewport(height=0.5, width=0.5, x=0.25,y=0.5) 
            scatterViewport <- viewport(height=0.5, width=0.5, x=0.75,y=0.5)
            densityViewport <- viewport(height=0.25,width=0.5, x=0.75,y=0.125)
            hmdensityViewport <- viewport(height=0.25,width=0.5,x=0.25,y=0.125)
            jit=position_jitter(width=0.5)
            hmrange=range(chrhmCmap[,1])
            bedIplot=bedIplot+hmrange[1]
            chrhm = ggplot(chrhmCmap, aes(x=Var1, y=Var2, fill=value))+scale_y_discrete(breaks=seq(0, 10, 5))+xlab('chrom')+ylab("chrom")+scale_fill_gradient(low='white', high='red')+geom_tile()+guides(fill=FALSE)
            chrbedplot=qplot(bedIplot[1,],bedIplot[2,],alpha=I(1/10),size=I(1))+xlab('chrom')+ylab("chrom")+geom_jitter(position=jit,colour="black",alpha=1/100)
            
            chrbeddensitydata=as.data.frame(c(chrBedToBedInter[,14],chrBedToBedInter[,15]))
            colnames(chrbeddensitydata)="bed"
            chrbeddensity=ggplot(chrbeddensitydata)+geom_density(aes(x=bed))
            
            chrhmdensitydata=NULL
            for(iiii in 1:chrTotSize)
            {
              #pp[iiii,1]=length(which(chrCmap[,iiii]>0))
              chrhmdensitydata=c(chrhmdensitydata,which(chrCmap[,iiii]>0))
            }
            chrhmdensitydata=as.data.frame(chrhmdensitydata)
            colnames(chrhmdensitydata)="chrom"
            
            
            chrcmapdensity=ggplot(chrhmdensitydata)+geom_density(aes(x=chrom))
            
            
            
            print(chrhm,vp=heatmapViewport)
            print(chrbedplot,vp=scatterViewport)
            print(chrbeddensity,vp=densityViewport)
            print(chrcmapdensity,vp=hmdensityViewport)
            
            dev.off()
            
            
            write.table(chrBedToBedInter,file=paste(genomeName,"/",chrName,"_BedToBedInter.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
            rm(chrCmap)
            rm(chrBed)
            rm(chrBedBin)
            rm(chrBedMatrix)
            rm(chrBedToBedInter)
            gc()
          }
        }
        
      }
      
    }
  }
  
}


#########################step3############################
if (3 %in% opt$stages)
{
  source("calbed.R")
  gc()
  matrix_dir=list.files(path=genomeName,full.names=F,pattern=".matrix")
  matrix_full_dir=list.files(path=genomeName,full.names=T,pattern=".matrix")
  m_bed=load_bed(bedFile)
  chrNum=length(matrix_dir)
  
  genomeFrame=data.frame("seg.name"=character(0),"seg.start"=numeric(0),"seg.end"=numeric(0),"the.v"=character(0),"NO"=character(0),stringsAsFactors=FALSE)
  genomeFrameNum=1
  for (i in 1:chrNum)
  {
    chrCmap=read.table(file=matrix_full_dir[i], fill=TRUE, stringsAsFactors=FALSE)
    chrTotSize=dim(chrCmap)[1]
    tmpNum=regexpr(".matrix",matrix_dir[i])
    chrName=substr(matrix_dir[i],1,tmpNum-1)
    
    if(chrName==chrom)
    {
      if(chrend>0)
      {
        tmpCmapStart=abs(ceiling((chrstart/(resolution*1000))))
        tmpCmapEnd=abs(ceiling((chrend/(resolution*1000))))
        if(tmpCmapEnd<tmpCmapStart)
        {
          print("please input correct start and end number")
          break
        }
        if(tmpCmapEnd>chrTotSize)
        {
          tmpCmapEnd=chrTotSize
        }
        if(tmpCmapStart>chrTotSize)
        {
          tmpCmapStart=chrTotSize
        }
        chrCmap=chrCmap[tmpCmapStart:tmpCmapEnd,tmpCmapStart:tmpCmapEnd]
        chrTotSize=dim(chrCmap)[1]
      }
    }
    
    for (ii in 1:chrTotSize)
    {
      genomeFrame[genomeFrameNum,1]=chrName
      genomeFrame[genomeFrameNum,2]=resolution*1000*(ii-1)
      genomeFrame[genomeFrameNum,3]=resolution*1000*ii
      genomeFrame[genomeFrameNum,4]=NA
      genomeFrame[genomeFrameNum,5]=NA
      genomeFrameNum=genomeFrameNum+1
    }
    
    
    rm(chrCmap)
    gc()
  }
  
  if(chrom=="all")
  {
    for (i in 1:chrNum)
    {
      
      tmpNum=regexpr(".matrix",matrix_dir[i])
      chrName=substr(matrix_dir[i],1,tmpNum-1)
      chrBed=choose_chr_bed(m_bed,chrName)
      chrBedNum=dim(chrBed)[1]
      if(chrBedNum>0)
      {
        chrCmap=read.table(file=matrix_full_dir[i], fill=TRUE, stringsAsFactors=FALSE)
        chrTotSize=dim(chrCmap)[1]
        tmpNum2=regexpr("r",chrName)
        chrNo=substr(chrName,tmpNum2+1,nchar(chrName))
        chrBedBin=check_bed_bin(chrBed,resolution*1000)
        chrBedMatrix=convert_bed_to_matrix(chrBed,resolution*1000,chrName,chrTotSize,bedWindow)
        chrCircosMapping=calculate_omiccircos_data(chrCmap,chrBedMatrix,chrBedBin,chrBed,chrName,chrCmap)
        chrCircosMapping[,1]=chrNo
        chrCircosMapping[,4]=chrNo
        chrCircosDb=segAnglePo(genomeFrame, seg=chrName)
        seg.num<-length(unique(genomeFrame[,1]))
        colors<-rainbow(seg.num, alpha=0.5)
        if((outputpdf=="TRUE")||(outputpdf=="true"))
        {
          pdf(paste(genomeName,"/",chrName,"_circos.pdf",sep=""),width = 8,height = 8)
          
        }else
        {
          jpeg(paste(genomeName,"/",chrName,"_circos.jpeg",sep=""),width=1000,height=1000,quality = 100)
        }
        par(mar=c(2, 2, 2, 2));
        plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
        if(circosLinecolor=="rainbow")
        {
          circosLinecolor=colors
        }
        circos(R=200, cir=chrCircosDb, W=40,  mapping=chrCircosMapping, type="link2", lwd=circosLineWidth,col=circosLinecolor);
        
        print(paste("plot ",chrName,"circos picture",sep=""))
        if((is.null(wigFile))==FALSE)
        {
          m_wig=load_wig(wigFile,resolution*1000,chrName,chrTotSize,0,0)
          chrBedWig=m_wig[which(chrBedMatrix[,1]==1),]
          chrBedWig[,2]=chrBedWig[,2]*resolution*1000
          m_wig[,2]=m_wig[,2]*resolution*1000
          
          m_wig_ttest=NULL
          if((length(chrBedWig))>0)
          {
            for(iii in 1:(dim(chrBedWig)[1]))
            {
              m_wig_ttest=rbind(m_wig_ttest,t.test(m_wig[,3],mu=chrBedWig[iii,3]))
            }
            m_wig_equal=which(m_wig_ttest[,3]>0.05)
            m_wig_mean=mean(m_wig[,3])
            m_wig_nequal=which(m_wig_ttest[,3]<=0.05)
            m_wig_more=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]>m_wig_mean)]
            m_wig_less=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]<=m_wig_mean)]
            m_wig_bed=chrBedWig
            m_wig_bed[,3]=0
            m_wig_bed[m_wig_more,3]=1
            m_wig_bed[m_wig_less,3]=-1
            
            circos(R=250, cir=chrCircosDb, W=20, mapping=m_wig, col.v=3, type="heatmap", lwd=circosHmSize);
            #circos(R=230, cir=chrCircosDb, W=20, mapping=chrBedWig, col.v=3, type="heatmap", lwd=circosHmSize);
            circos(R=230, cir=chrCircosDb, W=10, mapping=m_wig_bed, col.v=3, type="heatmap", lwd=circosCmpSize);
            text(400,680,family="mono",wigFile,cex=0.7)
          }
          
          
        }
        if((is.null(wigFile2))==FALSE)
        {
          m_wig=load_wig(wigFile2,resolution*1000,chrName,chrTotSize,0,0)
          chrBedWig=m_wig[which(chrBedMatrix[,1]==1),]
          chrBedWig[,2]=chrBedWig[,2]*resolution*1000
          m_wig[,2]=m_wig[,2]*resolution*1000
          
          m_wig_ttest=NULL
          if((length(chrBedWig))>0)
          {
            for(iii in 1:(dim(chrBedWig)[1]))
            {
              m_wig_ttest=rbind(m_wig_ttest,t.test(m_wig[,3],mu=chrBedWig[iii,3]))
            }
            m_wig_equal=which(m_wig_ttest[,3]>0.05)
            m_wig_mean=mean(m_wig[,3])
            m_wig_nequal=which(m_wig_ttest[,3]<=0.05)
            m_wig_more=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]>m_wig_mean)]
            m_wig_less=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]<=m_wig_mean)]
            m_wig_bed=chrBedWig
            m_wig_bed[,3]=0
            m_wig_bed[m_wig_more,3]=1
            m_wig_bed[m_wig_less,3]=-1
            
            circos(R=310, cir=chrCircosDb, W=20, mapping=m_wig, col.v=3, type="heatmap", lwd=circosHmSize);
            #circos(R=310, cir=chrCircosDb, W=20, mapping=chrBedWig, col.v=3, type="heatmap", lwd=circosHmSize);
            circos(R=290, cir=chrCircosDb, W=10, mapping=m_wig_bed, col.v=3, type="heatmap", lwd=circosCmpSize);
            text(400,750,family="mono",wigFile,cex=0.7)
          }
          
          
        }
        if((is.null(wigFile3))==FALSE)
        {
          m_wig=load_wig(wigFile3,resolution*1000,chrName,chrTotSize,0,0)
          chrBedWig=m_wig[which(chrBedMatrix[,1]==1),]
          chrBedWig[,2]=chrBedWig[,2]*resolution*1000
          m_wig[,2]=m_wig[,2]*resolution*1000
          if((length(chrBedWig))>0)
          {
            m_wig_ttest=NULL
            for(iii in 1:(dim(chrBedWig)[1]))
            {
              m_wig_ttest=rbind(m_wig_ttest,t.test(m_wig[,3],mu=chrBedWig[iii,3]))
            }
            m_wig_equal=which(m_wig_ttest[,3]>0.05)
            m_wig_mean=mean(m_wig[,3])
            m_wig_nequal=which(m_wig_ttest[,3]<=0.05)
            m_wig_more=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]>m_wig_mean)]
            m_wig_less=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]<=m_wig_mean)]
            m_wig_bed=chrBedWig
            m_wig_bed[,3]=0
            m_wig_bed[m_wig_more,3]=1
            m_wig_bed[m_wig_less,3]=-1
            
            circos(R=370, cir=chrCircosDb, W=20, mapping=m_wig, col.v=3, type="heatmap", lwd=circosHmSize);
            #circos(R=380, cir=chrCircosDb, W=20, mapping=chrBedWig, col.v=3, type="heatmap", lwd=circosHmSize);
            circos(R=350, cir=chrCircosDb, W=10, mapping=m_wig_bed, col.v=3, type="heatmap", lwd=circosCmpSize);
            text(400,820,family="mono",wigFile,cex=0.7)
          }
          
          
          
        }
        dev.off()
        
        rm(chrCmap)
        rm(chrBed)
        rm(chrBedBin)
        rm(chrBedMatrix)
        rm(chrCircosMapping)
        gc()
      }
    }
  }else
  {
    for (i in 1:chrNum)
    {
      
      tmpNum=regexpr(".matrix",matrix_dir[i])
      chrName=substr(matrix_dir[i],1,tmpNum-1)
      if(chrName==chrom)
      {
        chrBed=choose_chr_bed(m_bed,chrName)
        if(chrend>0)
        {
          chrBed=chrBed[which(chrBed[,3]<chrend),]
          chrBed[,2]=chrBed[,2]-chrstart
          chrBed[,3]=chrBed[,3]-chrstart
          chrBed=chrBed[which(chrBed[,2]>0),]
        }
        chrBedNum=dim(chrBed)[1]
        
        
        if(chrBedNum>0)
        {
          chrCmap=read.table(file=matrix_full_dir[i], fill=TRUE, stringsAsFactors=FALSE)
          
          chrTotSize=dim(chrCmap)[1]
          if(chrend>0)
          {
            tmpCmapStart=abs(ceiling((chrstart/(resolution*1000))))
            tmpCmapEnd=abs(ceiling((chrend/(resolution*1000))))
            if(tmpCmapEnd<tmpCmapStart)
            {
              print("please input correct start and end number")
              break
            }
            if(tmpCmapEnd>chrTotSize)
            {
              tmpCmapEnd=chrTotSize
            }
            if(tmpCmapStart>chrTotSize)
            {
              tmpCmapStart=chrTotSize
            }
            chrCmap=chrCmap[tmpCmapStart:tmpCmapEnd,tmpCmapStart:tmpCmapEnd]
            chrTotSize=dim(chrCmap)[1]
          }
          
          
          
          print(paste("plot ",chrName,"circos picture",sep=""))
          
          
          
          
          
          tmpNum2=regexpr("r",chrName)
          chrNo=substr(chrName,tmpNum2+1,nchar(chrName))
          #chrNo=as.numeric(chrNo)
          chrBedBin=check_bed_bin(chrBed,resolution*1000)
          chrBedMatrix=convert_bed_to_matrix(chrBed,resolution*1000,chrName,chrTotSize,bedWindow)
          
          
          
          chrCircosMapping=calculate_omiccircos_data(chrCmap,chrBedMatrix,chrBedBin,chrBed,chrName,chrCmap)
          chrCircosMapping[,1]=chrNo
          chrCircosMapping[,4]=chrNo
          chrCircosDb=segAnglePo(genomeFrame, seg=chrName)
          seg.num<-length(unique(genomeFrame[,1]))
          colors<-rainbow(seg.num, alpha=0.5)
          
          if((outputpdf=="TRUE")||(outputpdf=="true"))
          {
            pdf(paste(genomeName,"/",chrName,"_circos.pdf",sep=""),width = 8,height = 8)
            
          }else
          {
            jpeg(paste(genomeName,"/",chrName,"_circos.jpeg",sep=""),width=1000,height=1000,quality = 100)
          }
          par(mar=c(2, 2, 2, 2));
          plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
          if(circosLinecolor=="rainbow")
          {
            circosLinecolor=colors
          }
          circos(R=200, cir=chrCircosDb, W=40,  mapping=chrCircosMapping, type="link2", lwd=circosLineWidth,col=circosLinecolor);
          
          print(paste("plot ",chrName,"circos picture",sep=""))
          if((is.null(wigFile))==FALSE)
          {
            m_wig=load_wig(wigFile,resolution*1000,chrName,chrTotSize,chrstart,chrend)
            
            
            chrBedWig=m_wig[which(chrBedMatrix[,1]==1),]
            chrBedWig[,2]=chrBedWig[,2]*resolution*1000
            m_wig[,2]=m_wig[,2]*resolution*1000
            
            m_wig_ttest=NULL
            if((length(chrBedWig))>0)
            {
              for(iii in 1:(dim(chrBedWig)[1]))
              {
                m_wig_ttest=rbind(m_wig_ttest,t.test(m_wig[,3],mu=chrBedWig[iii,3]))
              }
              m_wig_equal=which(m_wig_ttest[,3]>0.05)
              m_wig_mean=mean(m_wig[,3])
              m_wig_nequal=which(m_wig_ttest[,3]<=0.05)
              m_wig_more=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]>m_wig_mean)]
              m_wig_less=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]<=m_wig_mean)]
              m_wig_bed=chrBedWig
              m_wig_bed[,3]=0
              m_wig_bed[m_wig_more,3]=1
              m_wig_bed[m_wig_less,3]=-1
              
              circos(R=250, cir=chrCircosDb, W=20, mapping=m_wig, col.v=3, type="heatmap", lwd=circosHmSize);
              #circos(R=230, cir=chrCircosDb, W=20, mapping=chrBedWig, col.v=3, type="heatmap", lwd=circosHmSize);
              circos(R=230, cir=chrCircosDb, W=10, mapping=m_wig_bed, col.v=3, type="heatmap", lwd=circosCmpSize);
              text(400,680,family="mono",wigFile,cex=0.7)
            }
            
            
          }
          if((is.null(wigFile2))==FALSE)
          {
            m_wig=load_wig(wigFile2,resolution*1000,chrName,chrTotSize,chrstart,chrend)
            
            
            chrBedWig=m_wig[which(chrBedMatrix[,1]==1),]
            chrBedWig[,2]=chrBedWig[,2]*resolution*1000
            m_wig[,2]=m_wig[,2]*resolution*1000
            
            m_wig_ttest=NULL
            if((length(chrBedWig))>0)
            {
              for(iii in 1:(dim(chrBedWig)[1]))
              {
                m_wig_ttest=rbind(m_wig_ttest,t.test(m_wig[,3],mu=chrBedWig[iii,3]))
              }
              m_wig_equal=which(m_wig_ttest[,3]>0.05)
              m_wig_mean=mean(m_wig[,3])
              m_wig_nequal=which(m_wig_ttest[,3]<=0.05)
              m_wig_more=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]>m_wig_mean)]
              m_wig_less=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]<=m_wig_mean)]
              m_wig_bed=chrBedWig
              m_wig_bed[,3]=0
              m_wig_bed[m_wig_more,3]=1
              m_wig_bed[m_wig_less,3]=-1
              
              circos(R=320, cir=chrCircosDb, W=20, mapping=m_wig, col.v=3, type="heatmap", lwd=circosHmSize);
              #circos(R=310, cir=chrCircosDb, W=20, mapping=chrBedWig, col.v=3, type="heatmap", lwd=circosHmSize);
              circos(R=300, cir=chrCircosDb, W=10, mapping=m_wig_bed, col.v=3, type="heatmap", lwd=circosCmpSize);
              text(400,750,family="mono",wigFile2,cex=0.7)
            }
            
          }
          if((is.null(wigFile3))==FALSE)
          {
            m_wig=load_wig(wigFile3,resolution*1000,chrName,chrTotSize,chrstart,chrend)
            
            
            
            chrBedWig=m_wig[which(chrBedMatrix[,1]==1),]
            chrBedWig[,2]=chrBedWig[,2]*resolution*1000
            m_wig[,2]=m_wig[,2]*resolution*1000
            
            m_wig_ttest=NULL
            if((length(chrBedWig))>0)
            {
              for(iii in 1:(dim(chrBedWig)[1]))
              {
                m_wig_ttest=rbind(m_wig_ttest,t.test(m_wig[,3],mu=chrBedWig[iii,3]))
              }
              m_wig_equal=which(m_wig_ttest[,3]>0.05)
              m_wig_mean=mean(m_wig[,3])
              m_wig_nequal=which(m_wig_ttest[,3]<=0.05)
              m_wig_more=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]>m_wig_mean)]
              m_wig_less=m_wig_nequal[which(chrBedWig[m_wig_nequal[],3]<=m_wig_mean)]
              m_wig_bed=chrBedWig
              m_wig_bed[,3]=0
              m_wig_bed[m_wig_more,3]=1
              m_wig_bed[m_wig_less,3]=-1
              
              circos(R=390, cir=chrCircosDb, W=20, mapping=m_wig, col.v=3, type="heatmap", lwd=circosHmSize);
              #circos(R=380, cir=chrCircosDb, W=20, mapping=chrBedWig, col.v=3, type="heatmap", lwd=circosHmSize);
              circos(R=370, cir=chrCircosDb, W=10, mapping=m_wig_bed, col.v=3, type="heatmap", lwd=circosCmpSize);
              text(400,820,family="mono",wigFile3,cex=0.7)
            }
            
          }
          dev.off()
          rm(chrCmap)
          rm(chrBed)
          rm(chrBedBin)
          rm(chrBedMatrix)
          rm(chrCircosMapping)
          gc()
        }
      }
    }
  }
}



#########################step4############################
if (4 %in% opt$stages)
{
  source("calbed.R")
  gc()
  matrix_dir=list.files(path=genomeName,full.names=F,pattern=".matrix")
  matrix_full_dir=list.files(path=genomeName,full.names=T,pattern=".matrix")
  m_bed=load_bed(bedFile)
  chrNum=length(matrix_dir)
  
  
  if(chrom=="all")
  {
    for (i in 1:chrNum)
    {
      tmpNum=regexpr(".matrix",matrix_dir[i])
      chrName=substr(matrix_dir[i],1,tmpNum-1)
      chrBed=choose_chr_bed(m_bed,chrName)
      chrBedNum=dim(chrBed)[1]
      print(matrix_full_dir[i])
      if(chrBedNum>0)
      {
        print(chrName)
        chrCmap=read.table(file=matrix_full_dir[i], fill=TRUE, stringsAsFactors=FALSE)
        chrTotSize=dim(chrCmap)[1]
        chrBedBin=check_bed_bin(chrBed,resolution*1000)
        chrBedMatrix=convert_bed_to_matrix(chrBed,resolution*1000,chrName,chrTotSize,bedWindow)
        
        wig_num=0
        clust_wig=NULL
        if((is.null(wigFile))==FALSE)
        {
          wig_num=1
          chrBedWig1=load_bed_wig(wigFile,chrBed,chrName,0,0,0)
          clust_wig=cbind(clust_wig,chrBedWig1[,4])
          all_wig1=load_all_wig(wigFile)
          all_wig1=all_wig1[all_wig1$chr==chrom,]
          all_wig1[,2]=as.numeric(all_wig1[,2])
          all_wig1[,3]=as.numeric(all_wig1[,3])
          all_wig1[,4]=as.numeric(all_wig1[,4])
          
        }
        if((is.null(wigFile2))==FALSE)
        {
          wig_num=wig_num+1
          chrBedWig2=load_bed_wig(wigFile2,chrBed,chrName,0,0,0)
          clust_wig=cbind(clust_wig,chrBedWig2[,4])
          all_wig2=load_all_wig(wigFile2)
          all_wig2=all_wig2[all_wig2$chr==chrom,]
          all_wig2[,2]=as.numeric(all_wig2[,2])
          all_wig2[,3]=as.numeric(all_wig2[,3])
          all_wig2[,4]=as.numeric(all_wig2[,4])
        }
        if((is.null(wigFile3))==FALSE)
        {
          wig_num=wig_num+1
          chrBedWig3=load_bed_wig(wigFile3,chrBed,chrName,0,0,0)
          clust_wig=cbind(clust_wig,chrBedWig3[,4])
          all_wig3=load_all_wig(wigFile3)
          all_wig3=all_wig3[all_wig3$chr==chrom,]
          all_wig3[,2]=as.numeric(all_wig3[,2])
          all_wig3[,3]=as.numeric(all_wig3[,3])
          all_wig3[,4]=as.numeric(all_wig3[,4])
          
        }
        
        tt=c(1:(dim(clust_wig)[1]))
        clust_name=paste(tt,"_",chrBed[,1],":",(chrBed[,2]+chrstart),"-",(chrBed[,3]+chrstart),sep = "")
        row.names(clust_wig)=clust_name
        
        colnames(clust_wig)=c(wigFile,wigFile2,wigFile3)
        mydata = na.omit(clust_wig)
        suppressPackageStartupMessages(library("lattice"))
        suppressPackageStartupMessages(library("flexclust"))
        if(wig_num>1)
        {
          bcl <- bootFlexclust(mydata, k=2:7, nboot=50, FUN=cclust, multicore=FALSE)
          if((outputpdf=="TRUE")||(outputpdf=="true"))
          {
            pdf(paste(genomeName,"/",chrName,"_cluster_k_density.pdf",sep=""),width = 8,height = 8)
            
          }else
          {
            jpeg(paste(genomeName,"/",chrName,"_cluster_k_density.jpeg",sep=""),width=1000,height=1000,quality = 100)
          }
          plot(bcl)
          
          densityplot(bcl, from=0)
          
          dev.off()
        }
        
        out.dist=dist(mydata,method=dist_method) #manhattan,euclidean,minkowski,chebyshev,mahalanobis,canberra
        out.hclust=hclust(out.dist,method=clust_method) #average,centroid,median,complete,single,ward.D,density
        if((outputpdf=="TRUE")||(outputpdf=="true"))
        {
          pdf(paste(genomeName,"/",chrName,"_cluster_tree.pdf",sep=""),width = 8,height = 8)
          
        }else
        {
          jpeg(paste(genomeName,"/",chrName,"_cluster_tree.jpeg",sep=""),width=1000,height=1000,quality = 100)
        }
        plot(out.hclust)
        rect.hclust(out.hclust,clust_k)  
        cluster.id=cutree(out.hclust,clust_k)  
        dev.off()
        row.names(chrBed)=clust_name
        
        
        m_ttest_result=NULL
        m_ttest1_result=data.frame("wig1_pvalue"=numeric(dim(clust_wig)[1]),"wig1_difference"=character(dim(clust_wig)[1]),stringsAsFactors=FALSE)
        m_ttest2_result=data.frame("wig2_pvalue"=numeric(dim(clust_wig)[1]),"wig2_difference"=character(dim(clust_wig)[1]),stringsAsFactors=FALSE)
        m_ttest3_result=data.frame("wig3_pvalue"=numeric(dim(clust_wig)[1]),"wig3_difference"=character(dim(clust_wig)[1]),stringsAsFactors=FALSE)
        
        m_wig1_ttest=NULL
        m_wig2_ttest=NULL
        m_wig3_ttest=NULL
        
        logfile = paste(paste(genomeName,"/",chrName,"_statistic.txt",sep=""))
        if (file.exists(logfile) ==TRUE){file.remove(logfile)}
        starttime = paste("Analysis start time:" , as.character(Sys.time()))
        write(starttime,file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        
        
        if((is.null(wigFile))==FALSE)
        {
          for(iii in 1:(dim(clust_wig)[1]))
          {
            if(is.na(clust_wig[iii,1])==FALSE)
            {
              m_wig1_ttest=rbind(m_wig1_ttest,t.test(all_wig1[,4],mu=clust_wig[iii,1]))
            }else
            {
              m_wig1_ttest=rbind(m_wig1_ttest,t.test(all_wig1[,4],mu=0))
            }
          }
          m_wig1_equal=which(m_wig1_ttest[,3]>0.05)
          m_wig1_mean=mean(all_wig1[,4])
          m_wig1_nequal=which(m_wig1_ttest[,3]<=0.05)
          m_wig1_more=m_wig1_nequal[which(clust_wig[m_wig1_nequal[],1]>m_wig1_mean)]
          m_wig1_less=m_wig1_nequal[which(clust_wig[m_wig1_nequal[],1]<=m_wig1_mean)]
          m_ttest1_result[,1]=as.data.frame(as.matrix(m_wig1_ttest[,3]))
          m_ttest1_result[m_wig1_more[],2]="more"
          m_ttest1_result[m_wig1_less[],2]="less"
          m_ttest1_result[m_wig1_equal[],2]="equal"
          if(is.null(m_ttest_result))
          {
            m_ttest_result=m_ttest1_result
          }else
          {
            m_ttest_result=cbind(m_ttest_result,m_ttest1_result)
            
          }
          wig1_test=rbind(cbind(na.omit(clust_wig[,1]),1),cbind(na.omit(all_wig1[,4]),2))
          wig1_test=as.data.frame(wig1_test)
          colnames(wig1_test)=c("wig1_value","group")
          rownames(wig1_test)=c(1:(dim(wig1_test)[1]))
          wig1_test$group=as.factor(wig1_test$group)
          wig1_kruskal=kruskal.test(wig1_value~group, data=wig1_test)
          wig1_kruskalmc=kruskalmc(wig1_value~group, data=wig1_test, probs=0.05)
          wig1_mult <- oneway_test(wig1_value~group, data=wig1_test,
                                   ytrafo = function(data) trafo(data, numeric_trafo = rank),
                                   xtrafo = function(data) trafo(data, factor_trafo = function(x)
                                     model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                                   teststat = "max", distribution = approximate(B = 90000))
          wig1_pvalue=pvalue(wig1_mult, method = "single-step")
          write("the statistic test between BED WIG1 and ALL WIG1 :",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("test name : Kruskal-Wallis rank sum test",file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis chi-squared : ",wig1_kruskal$statistic,sep = ""),file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis df : ",wig1_kruskal$parameter,sep = ""),file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis p value : ",wig1_kruskal$p.value,sep = ""),file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("test name : Multiple comparison test after Kruskal-Wallis",file=logfile,append=TRUE)
          write(paste("significance level : ",wig1_kruskalmc$signif.level,sep = ""),file=logfile,append=TRUE)
          write(paste("observed difference  : ",wig1_kruskalmc$dif.com$obs.dif,sep = ""),file=logfile,append=TRUE)
          write(paste("critical difference  : ",wig1_kruskalmc$dif.com$critical.dif,sep = ""),file=logfile,append=TRUE)
          write(paste("exist difference  : ",wig1_kruskalmc$dif.com$difference,sep = ""),file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
        }
        if((is.null(wigFile2))==FALSE)
        {
          for(iii in 1:(dim(clust_wig)[1]))
          {
            if(is.na(clust_wig[iii,2])==FALSE)
            {
              m_wig2_ttest=rbind(m_wig2_ttest,t.test(all_wig2[,4],mu=clust_wig[iii,2]))
            }else
            {
              m_wig2_ttest=rbind(m_wig2_ttest,t.test(all_wig2[,4],mu=0))
            }
          }
          m_wig2_equal=which(m_wig2_ttest[,3]>0.05)
          m_wig2_mean=mean(all_wig2[,4])
          m_wig2_nequal=which(m_wig2_ttest[,3]<=0.05)
          m_wig2_more=m_wig2_nequal[which(clust_wig[m_wig2_nequal[],2]>m_wig2_mean)]
          m_wig2_less=m_wig2_nequal[which(clust_wig[m_wig2_nequal[],2]<=m_wig2_mean)]
          m_ttest2_result[,1]=as.data.frame(as.matrix(m_wig2_ttest[,3]))
          m_ttest2_result[m_wig2_more[],2]="more"
          m_ttest2_result[m_wig2_less[],2]="less"
          m_ttest2_result[m_wig2_equal[],2]="equal"
          if(is.null(m_ttest_result))
          {
            m_ttest_result=m_ttest2_result
          }else
          {
            m_ttest_result=cbind(m_ttest_result,m_ttest2_result)
            
          }
          wig2_test=rbind(cbind(na.omit(clust_wig[,2]),1),cbind(na.omit(all_wig2[,4]),2))
          wig2_test=as.data.frame(wig2_test)
          colnames(wig2_test)=c("wig2_value","group")
          rownames(wig2_test)=c(1:(dim(wig2_test)[1]))
          wig2_test$group=as.factor(wig2_test$group)
          wig2_kruskal=kruskal.test(wig2_value~group, data=wig2_test)
          wig2_kruskalmc=kruskalmc(wig2_value~group, data=wig2_test, probs=0.05)
          wig2_mult <- oneway_test(wig2_value~group, data=wig2_test,
                                   ytrafo = function(data) trafo(data, numeric_trafo = rank),
                                   xtrafo = function(data) trafo(data, factor_trafo = function(x)
                                     model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                                   teststat = "max", distribution = approximate(B = 90000))
          wig2_pvalue=pvalue(wig2_mult, method = "single-step")
          write("the statistic test between BED WIG2 and ALL WIG2 :",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("test name : Kruskal-Wallis rank sum test",file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis chi-squared : ",wig2_kruskal$statistic,sep = ""),file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis df : ",wig2_kruskal$parameter,sep = ""),file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis p value : ",wig2_kruskal$p.value,sep = ""),file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("test name : Multiple comparison test after Kruskal-Wallis",file=logfile,append=TRUE)
          write(paste("significance level : ",wig2_kruskalmc$signif.level,sep = ""),file=logfile,append=TRUE)
          write(paste("observed difference  : ",wig2_kruskalmc$dif.com$obs.dif,sep = ""),file=logfile,append=TRUE)
          write(paste("critical difference  : ",wig2_kruskalmc$dif.com$critical.dif,sep = ""),file=logfile,append=TRUE)
          write(paste("exist difference  : ",wig2_kruskalmc$dif.com$difference,sep = ""),file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
        }
        if((is.null(wigFile3))==FALSE)
        {
          for(iii in 1:(dim(clust_wig)[1]))
          {
            if(is.na(clust_wig[iii,3])==FALSE)
            {
              m_wig3_ttest=rbind(m_wig3_ttest,t.test(all_wig3[,4],mu=clust_wig[iii,3]))
            }else
            {
              m_wig3_ttest=rbind(m_wig3_ttest,t.test(all_wig3[,4],mu=0))
            }
          }
          m_wig3_equal=which(m_wig3_ttest[,3]>0.05)
          m_wig3_mean=mean(all_wig3[,4])
          m_wig3_nequal=which(m_wig3_ttest[,3]<=0.05)
          m_wig3_more=m_wig3_nequal[which(clust_wig[m_wig3_nequal[],3]>m_wig3_mean)]
          m_wig3_less=m_wig3_nequal[which(clust_wig[m_wig3_nequal[],3]<=m_wig3_mean)]
          m_ttest3_result[,1]=as.data.frame(as.matrix(m_wig3_ttest[,3]))
          m_ttest3_result[m_wig3_more[],2]="more"
          m_ttest3_result[m_wig3_less[],2]="less"
          m_ttest3_result[m_wig3_equal[],2]="equal"
          if(is.null(m_ttest_result))
          {
            m_ttest_result=m_ttest3_result
          }else
          {
            m_ttest_result=cbind(m_ttest_result,m_ttest3_result)
            
          }
          wig3_test=rbind(cbind(na.omit(clust_wig[,3]),1),cbind(na.omit(all_wig3[,4]),2))
          wig3_test=as.data.frame(wig3_test)
          colnames(wig3_test)=c("wig3_value","group")
          rownames(wig3_test)=c(1:(dim(wig3_test)[1]))
          wig3_test$group=as.factor(wig3_test$group)
          wig3_kruskal=kruskal.test(wig3_value~group, data=wig3_test)
          wig3_kruskalmc=kruskalmc(wig3_value~group, data=wig3_test, probs=0.05)
          wig3_mult <- oneway_test(wig3_value~group, data=wig3_test,
                                   ytrafo = function(data) trafo(data, numeric_trafo = rank),
                                   xtrafo = function(data) trafo(data, factor_trafo = function(x)
                                     model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                                   teststat = "max", distribution = approximate(B = 90000))
          wig3_pvalue=pvalue(wig3_mult, method = "single-step")
          write("the statistic test between BED WIG3 and ALL WIG3 :",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("test name : Kruskal-Wallis rank sum test",file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis chi-squared : ",wig3_kruskal$statistic,sep = ""),file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis df : ",wig3_kruskal$parameter,sep = ""),file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis p value : ",wig3_kruskal$p.value,sep = ""),file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("test name : Multiple comparison test after Kruskal-Wallis",file=logfile,append=TRUE)
          write(paste("significance level : ",wig3_kruskalmc$signif.level,sep = ""),file=logfile,append=TRUE)
          write(paste("observed difference  : ",wig3_kruskalmc$dif.com$obs.dif,sep = ""),file=logfile,append=TRUE)
          write(paste("critical difference  : ",wig3_kruskalmc$dif.com$critical.dif,sep = ""),file=logfile,append=TRUE)
          write(paste("exist difference  : ",wig3_kruskalmc$dif.com$difference,sep = ""),file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
        }
        
        
        
        
        
        clust_group=cbind(na.omit(cbind(chrBed,clust_wig,m_ttest_result)),cluster.id)
        write.csv(clust_group,paste(genomeName,"/",chrName,"_cluster.csv",sep=""),row.names = FALSE)
        clust_heatmap=NULL
        for( ii in 1:clust_k)
        {
          clust_heatmap=rbind(clust_heatmap,clust_group[which(clust_group[,"cluster.id"]==ii),])
        }
        rownames(clust_heatmap)=NULL
        if((outputpdf=="TRUE")||(outputpdf=="true"))
        {
          pdf(paste(genomeName,"/",chrName,"_cluster_heatmap.pdf",sep=""),width = 8,height = 8)
          
        }else
        {
          jpeg(paste(genomeName,"/",chrName,"_cluster_heatmap.jpeg",sep=""),width=1000,height=1000,quality = 100)
        }
        heatmap(as.matrix(clust_heatmap[,5:(5+wig_num-1)]),Rowv=NA,Colv=NA,cexCol = 1,labCol = "")
        dev.off()
        
        
        n_count=0
        for (ii in 1:chrTotSize)
        {
          if(chrBedMatrix[ii]!=0)
          {
            n_count=n_count+1
          }
        }
        random_group=matrix(data=0, nrow = groupNum, ncol = chrTotSize)
        random_result=matrix(data=0, nrow = groupNum , ncol = 3)
        
        for (ii in 1:groupNum)
        {
          tmp_site=sample(1:chrTotSize,size=n_count)
          random_group[ii,tmp_site]= 2
        }
        
        tmpCmap=chrCmap
        tmpCmap[lower.tri(tmpCmap)]=NA
        
        for(t in 1:groupNum)
        {
          
          random_result[t,1]=0
          random_result[t,2]=0
          for(ii in 1:chrTotSize)
          {
            if(random_group[t,ii]>0)
            {
              b=which(tmpCmap[ii,]>0)
              c=which((random_group[t,b])>0)
              random_result[t,1]=random_result[t,1]+length(c)
              random_result[t,2]=random_result[t,2]+length(b)
            }
          }
          
        }
        
        bed_count=0
        all_count=0
        bb_count=0
        bb_info=NULL
        bedTOall_info=NULL
        
        for(ii in 1:chrTotSize)
        {
          if(chrBedMatrix[ii]>0)
          {
            b=which(tmpCmap[ii,]>0)
            c=which((chrBedMatrix[b])>0)
            if(length(b)>0)
            {
              bedTOall_info=rbind(bedTOall_info,cbind(ii,b,t(tmpCmap[ii,b])))
            }
            if(length(c)>0)
            {
              bb_info=rbind(bb_info,cbind(ii,b[c],t(tmpCmap[ii,b[c]])))
            }
            bb_count=bb_count+length(c)
            bed_count=bed_count+length(b)
          }
        }
        
        
        
        if_test=rbind(cbind(bb_info[,3],1),cbind(bedTOall_info[,3],2))
        if_test=as.data.frame(if_test)
        colnames(if_test)=c("if_value","group")
        rownames(if_test)=c(1:(dim(if_test)[1]))
        if_test$group=as.factor(if_test$group)
        if_kruskal=kruskal.test(if_value~group, data=if_test)
        if_kruskalmc=kruskalmc(if_value~group, data=if_test, probs=0.05)
        mult <- oneway_test(if_value~group, data=if_test,
                            ytrafo = function(data) trafo(data, numeric_trafo = rank),
                            xtrafo = function(data) trafo(data, factor_trafo = function(x)
                              model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                            teststat = "max", distribution = approximate(B = 90000))
        if_pvalue=pvalue(mult, method = "single-step")
        
        
        dif_result=t.test(random_result[,1],mu=bb_count)
        
        
        
        write("",file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("the statistic test of interaction frequency between b2b and b2o :",file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("test name : Kruskal-Wallis rank sum test",file=logfile,append=TRUE)
        write(paste("Kruskal-Wallis chi-squared : ",if_kruskal$statistic,sep = ""),file=logfile,append=TRUE)
        write(paste("Kruskal-Wallis df : ",if_kruskal$parameter,sep = ""),file=logfile,append=TRUE)
        write(paste("Kruskal-Wallis p value : ",if_kruskal$p.value,sep = ""),file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("test name : Multiple comparison test after Kruskal-Wallis",file=logfile,append=TRUE)
        write(paste("significance level : ",if_kruskalmc$signif.level,sep = ""),file=logfile,append=TRUE)
        write(paste("observed difference  : ",if_kruskalmc$dif.com$obs.dif,sep = ""),file=logfile,append=TRUE)
        write(paste("critical difference  : ",if_kruskalmc$dif.com$critical.dif,sep = ""),file=logfile,append=TRUE)
        write(paste("exist difference  : ",if_kruskalmc$dif.com$difference,sep = ""),file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("the statistic test of interaction number between b2b and o2o :",file=logfile,append=TRUE)
        write("",file=logfile,append=TRUE)
        write("test name : t-test",file=logfile,append=TRUE)
        write(paste("numbers of random group : ",groupNum,sep = ""),file=logfile,append=TRUE)
        write(paste("95 percent confidence interval of random group : ",dif_result$conf.int[1]," ~ ",dif_result$conf.int[2],sep = ""),file=logfile,append=TRUE)
        write(paste("numbers of b2b : ",bb_count,sep = ""),file=logfile,append=TRUE)
        write(paste("t test p value : ",dif_result$p.value,sep = ""),file=logfile,append=TRUE)
        if(dif_result$p.value<0.05)
        {
          write("exist difference  : TRUE",file=logfile,append=TRUE)
          
        }else
        {
          write("exist difference  : FALSE",file=logfile,append=TRUE)
          
        }
      }
    }
  }else
  {
    for (i in 1:chrNum)
    {
      tmpNum=regexpr(".matrix",matrix_dir[i])
      chrName=substr(matrix_dir[i],1,tmpNum-1)
      if(chrName==chrom)
      {
        chrBed=choose_chr_bed(m_bed,chrName)
        if(chrend>0)
        {
          chrBed=chrBed[which(chrBed[,3]<chrend),]
          chrBed[,2]=chrBed[,2]-chrstart
          chrBed[,3]=chrBed[,3]-chrstart
          chrBed=chrBed[which(chrBed[,2]>0),]
        }
        chrBedNum=dim(chrBed)[1]
        print(matrix_full_dir[i])
        if(chrBedNum>0)
        {
          print(chrName)
          chrCmap=read.table(file=matrix_full_dir[i], fill=TRUE, stringsAsFactors=FALSE)
          chrTotSize=dim(chrCmap)[1]
          
          if(chrend>0)
          {
            tmpCmapStart=abs(ceiling((chrstart/(resolution*1000))))
            tmpCmapEnd=abs(ceiling((chrend/(resolution*1000))))
            if(tmpCmapEnd<tmpCmapStart)
            {
              print("please input correct start and end number")
              break
            }
            if(tmpCmapEnd>chrTotSize)
            {
              tmpCmapEnd=chrTotSize
            }
            if(tmpCmapStart>chrTotSize)
            {
              tmpCmapStart=chrTotSize
            }
            chrCmap=chrCmap[tmpCmapStart:tmpCmapEnd,tmpCmapStart:tmpCmapEnd]
            chrTotSize=dim(chrCmap)[1]
          }
          
          
          
          chrBedBin=check_bed_bin(chrBed,resolution*1000)
          chrBedMatrix=convert_bed_to_matrix(chrBed,resolution*1000,chrName,chrTotSize,bedWindow)
          
          
          
          
          wig_num=0
          clust_wig=NULL
          if((is.null(wigFile))==FALSE)
          {
            wig_num=1
            chrBedWig1=load_bed_wig(wigFile,chrBed,chrName,0,0,0)
            clust_wig=cbind(clust_wig,chrBedWig1[,4])
            all_wig1=load_all_wig(wigFile)
            all_wig1=all_wig1[all_wig1$chr==chrom,]
            all_wig1[,2]=as.numeric(all_wig1[,2])
            all_wig1[,3]=as.numeric(all_wig1[,3])
            all_wig1[,4]=as.numeric(all_wig1[,4])
            if(chrend>0)
            {
              all_wig1=all_wig1[which(all_wig1[,3]<chrend),]
            }
            if(chrstart>0)
            {
              all_wig1[,2]=all_wig1[,2]-chrstart
              all_wig1[,3]=all_wig1[,3]-chrstart
              all_wig1=all_wig1[which(all_wig1[,2]>=0),]
            }
            
          }
          if((is.null(wigFile2))==FALSE)
          {
            wig_num=wig_num+1
            chrBedWig2=load_bed_wig(wigFile2,chrBed,chrName,0,0,0)
            clust_wig=cbind(clust_wig,chrBedWig2[,4])
            all_wig2=load_all_wig(wigFile2)
            all_wig2=all_wig2[all_wig2$chr==chrom,]
            all_wig2[,2]=as.numeric(all_wig2[,2])
            all_wig2[,3]=as.numeric(all_wig2[,3])
            all_wig2[,4]=as.numeric(all_wig2[,4])
            if(chrend>0)
            {
              all_wig2=all_wig1[which(all_wig2[,3]<chrend),]
            }
            if(chrstart>0)
            {
              all_wig2[,2]=all_wig2[,2]-chrstart
              all_wig2[,3]=all_wig2[,3]-chrstart
              all_wig2=all_wig2[which(all_wig2[,2]>=0),]
            }
          }
          if((is.null(wigFile3))==FALSE)
          {
            wig_num=wig_num+1
            chrBedWig3=load_bed_wig(wigFile3,chrBed,chrName,0,0,0)
            clust_wig=cbind(clust_wig,chrBedWig3[,4])
            all_wig3=load_all_wig(wigFile3)
            all_wig3=all_wig3[all_wig3$chr==chrom,]
            all_wig3[,2]=as.numeric(all_wig3[,2])
            all_wig3[,3]=as.numeric(all_wig3[,3])
            all_wig3[,4]=as.numeric(all_wig3[,4])
            if(chrend>0)
            {
              all_wig3=all_wig3[which(all_wig3[,3]<chrend),]
            }
            if(chrstart>0)
            {
              all_wig3[,2]=all_wig3[,2]-chrstart
              all_wig3[,3]=all_wig3[,3]-chrstart
              all_wig3=all_wig3[which(all_wig3[,2]>=0),]
            }
          }
          
          tt=c(1:(dim(clust_wig)[1]))
          clust_name=paste(tt,"_",chrBed[,1],":",(chrBed[,2]+chrstart),"-",(chrBed[,3]+chrstart),sep = "")
          row.names(clust_wig)=clust_name
          
          colnames(clust_wig)=c(wigFile,wigFile2,wigFile3)
          mydata = na.omit(clust_wig)
          
          suppressPackageStartupMessages(library("lattice"))
          suppressPackageStartupMessages(library("flexclust"))
          if(wig_num>1)
          {
            bcl <- bootFlexclust(mydata, k=2:7, nboot=50, FUN=cclust, multicore=FALSE)
            if((outputpdf=="TRUE")||(outputpdf=="true"))
            {
              pdf(paste(genomeName,"/",chrName,"_cluster_k_density.pdf",sep=""),width = 8,height = 8)
              
            }else
            {
              jpeg(paste(genomeName,"/",chrName,"_cluster_k_density.jpeg",sep=""),width=1000,height=1000,quality = 100)
            }
            plot(bcl)
            densityplot(bcl, from=0)
            
            dev.off()
          }
          
          out.dist=dist(mydata,method=dist_method) #manhattan,euclidean,minkowski,chebyshev,mahalanobis,canberra
          out.hclust=hclust(out.dist,method=clust_method) #average,centroid,median,complete,single,ward.D,density
          if((outputpdf=="TRUE")||(outputpdf=="true"))
          {
            pdf(paste(genomeName,"/",chrName,"_cluster_tree.pdf",sep=""),width = 8,height = 8)
            
          }else
          {
            jpeg(paste(genomeName,"/",chrName,"_cluster_tree.jpeg",sep=""),width=1000,height=1000,quality = 100)
          }
          plot(out.hclust)
          rect.hclust(out.hclust,clust_k)  
          cluster.id=cutree(out.hclust,clust_k)  
          dev.off()
          row.names(chrBed)=clust_name
          
          
          m_ttest_result=NULL
          m_ttest1_result=data.frame("wig1_pvalue"=numeric(dim(clust_wig)[1]),"wig1_difference"=character(dim(clust_wig)[1]),stringsAsFactors=FALSE)
          m_ttest2_result=data.frame("wig2_pvalue"=numeric(dim(clust_wig)[1]),"wig2_difference"=character(dim(clust_wig)[1]),stringsAsFactors=FALSE)
          m_ttest3_result=data.frame("wig3_pvalue"=numeric(dim(clust_wig)[1]),"wig3_difference"=character(dim(clust_wig)[1]),stringsAsFactors=FALSE)
          
          m_wig1_ttest=NULL
          m_wig2_ttest=NULL
          m_wig3_ttest=NULL
          
          logfile = paste(paste(genomeName,"/",chrName,"_statistic.txt",sep=""))
          if (file.exists(logfile) ==TRUE){file.remove(logfile)}
          starttime = paste("Analysis start time:" , as.character(Sys.time()))
          write(starttime,file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          
          
          if((is.null(wigFile))==FALSE)
          {
            for(iii in 1:(dim(clust_wig)[1]))
            {
              if(is.na(clust_wig[iii,1])==FALSE)
              {
                m_wig1_ttest=rbind(m_wig1_ttest,t.test(all_wig1[,4],mu=clust_wig[iii,1]))
              }else
              {
                m_wig1_ttest=rbind(m_wig1_ttest,t.test(all_wig1[,4],mu=0))
              }
            }
            m_wig1_equal=which(m_wig1_ttest[,3]>0.05)
            m_wig1_mean=mean(all_wig1[,4])
            m_wig1_nequal=which(m_wig1_ttest[,3]<=0.05)
            m_wig1_more=m_wig1_nequal[which(clust_wig[m_wig1_nequal[],1]>m_wig1_mean)]
            m_wig1_less=m_wig1_nequal[which(clust_wig[m_wig1_nequal[],1]<=m_wig1_mean)]
            m_ttest1_result[,1]=as.data.frame(as.matrix(m_wig1_ttest[,3]))
            m_ttest1_result[,1]=as.numeric(m_ttest1_result[,1])
            m_ttest1_result[m_wig1_more[],2]="more"
            m_ttest1_result[m_wig1_less[],2]="less"
            m_ttest1_result[m_wig1_equal[],2]="equal"
            if(is.null(m_ttest_result))
            {
              m_ttest_result=m_ttest1_result
            }else
            {
              m_ttest_result=cbind(m_ttest_result,m_ttest1_result)
              
            }
            wig1_test=rbind(cbind(na.omit(clust_wig[,1]),1),cbind(na.omit(all_wig1[,4]),2))
            wig1_test=as.data.frame(wig1_test)
            colnames(wig1_test)=c("wig1_value","group")
            rownames(wig1_test)=c(1:(dim(wig1_test)[1]))
            wig1_test$group=as.factor(wig1_test$group)
            wig1_kruskal=kruskal.test(wig1_value~group, data=wig1_test)
            wig1_kruskalmc=kruskalmc(wig1_value~group, data=wig1_test, probs=0.05)
            wig1_mult <- oneway_test(wig1_value~group, data=wig1_test,
                                     ytrafo = function(data) trafo(data, numeric_trafo = rank),
                                     xtrafo = function(data) trafo(data, factor_trafo = function(x)
                                       model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                                     teststat = "max", distribution = approximate(B = 90000))
            wig1_pvalue=pvalue(wig1_mult, method = "single-step")
            write("the statistic test between BED WIG1 and ALL WIG1 :",file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("test name : Kruskal-Wallis rank sum test",file=logfile,append=TRUE)
            write(paste("Kruskal-Wallis chi-squared : ",wig1_kruskal$statistic,sep = ""),file=logfile,append=TRUE)
            write(paste("Kruskal-Wallis df : ",wig1_kruskal$parameter,sep = ""),file=logfile,append=TRUE)
            write(paste("Kruskal-Wallis p value : ",wig1_kruskal$p.value,sep = ""),file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("test name : Multiple comparison test after Kruskal-Wallis",file=logfile,append=TRUE)
            write(paste("significance level : ",wig1_kruskalmc$signif.level,sep = ""),file=logfile,append=TRUE)
            write(paste("observed difference  : ",wig1_kruskalmc$dif.com$obs.dif,sep = ""),file=logfile,append=TRUE)
            write(paste("critical difference  : ",wig1_kruskalmc$dif.com$critical.dif,sep = ""),file=logfile,append=TRUE)
            write(paste("exist difference  : ",wig1_kruskalmc$dif.com$difference,sep = ""),file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
          }
          if((is.null(wigFile2))==FALSE)
          {
            for(iii in 1:(dim(clust_wig)[1]))
            {
              if(is.na(clust_wig[iii,2])==FALSE)
              {
                m_wig2_ttest=rbind(m_wig2_ttest,t.test(all_wig2[,4],mu=clust_wig[iii,2]))
              }else
              {
                m_wig2_ttest=rbind(m_wig2_ttest,t.test(all_wig2[,4],mu=0))
              }
            }
            m_wig2_equal=which(m_wig2_ttest[,3]>0.05)
            m_wig2_mean=mean(all_wig2[,4])
            m_wig2_nequal=which(m_wig2_ttest[,3]<=0.05)
            m_wig2_more=m_wig2_nequal[which(clust_wig[m_wig2_nequal[],2]>m_wig2_mean)]
            m_wig2_less=m_wig2_nequal[which(clust_wig[m_wig2_nequal[],2]<=m_wig2_mean)]
            m_ttest2_result[,1]=as.data.frame(as.matrix(m_wig2_ttest[,3]))
            m_ttest2_result[,1]=as.numeric(m_ttest2_result[,1])
            m_ttest2_result[m_wig2_more[],2]="more"
            m_ttest2_result[m_wig2_less[],2]="less"
            m_ttest2_result[m_wig2_equal[],2]="equal"
            if(is.null(m_ttest_result))
            {
              m_ttest_result=m_ttest2_result
            }else
            {
              m_ttest_result=cbind(m_ttest_result,m_ttest2_result)
              
            }
            wig2_test=rbind(cbind(na.omit(clust_wig[,2]),1),cbind(na.omit(all_wig2[,4]),2))
            wig2_test=as.data.frame(wig2_test)
            colnames(wig2_test)=c("wig2_value","group")
            rownames(wig2_test)=c(1:(dim(wig2_test)[1]))
            wig2_test$group=as.factor(wig2_test$group)
            wig2_kruskal=kruskal.test(wig2_value~group, data=wig2_test)
            wig2_kruskalmc=kruskalmc(wig2_value~group, data=wig2_test, probs=0.05)
            wig2_mult <- oneway_test(wig2_value~group, data=wig2_test,
                                     ytrafo = function(data) trafo(data, numeric_trafo = rank),
                                     xtrafo = function(data) trafo(data, factor_trafo = function(x)
                                       model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                                     teststat = "max", distribution = approximate(B = 90000))
            wig2_pvalue=pvalue(wig2_mult, method = "single-step")
            write("the statistic test between BED WIG2 and ALL WIG2 :",file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("test name : Kruskal-Wallis rank sum test",file=logfile,append=TRUE)
            write(paste("Kruskal-Wallis chi-squared : ",wig2_kruskal$statistic,sep = ""),file=logfile,append=TRUE)
            write(paste("Kruskal-Wallis df : ",wig2_kruskal$parameter,sep = ""),file=logfile,append=TRUE)
            write(paste("Kruskal-Wallis p value : ",wig2_kruskal$p.value,sep = ""),file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("test name : Multiple comparison test after Kruskal-Wallis",file=logfile,append=TRUE)
            write(paste("significance level : ",wig2_kruskalmc$signif.level,sep = ""),file=logfile,append=TRUE)
            write(paste("observed difference  : ",wig2_kruskalmc$dif.com$obs.dif,sep = ""),file=logfile,append=TRUE)
            write(paste("critical difference  : ",wig2_kruskalmc$dif.com$critical.dif,sep = ""),file=logfile,append=TRUE)
            write(paste("exist difference  : ",wig2_kruskalmc$dif.com$difference,sep = ""),file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
          }
          if((is.null(wigFile3))==FALSE)
          {
            for(iii in 1:(dim(clust_wig)[1]))
            {
              if(is.na(clust_wig[iii,3])==FALSE)
              {
                m_wig3_ttest=rbind(m_wig3_ttest,t.test(all_wig3[,4],mu=clust_wig[iii,3]))
              }else
              {
                m_wig3_ttest=rbind(m_wig3_ttest,t.test(all_wig3[,4],mu=0))
              }
            }
            m_wig3_equal=which(m_wig3_ttest[,3]>0.05)
            m_wig3_mean=mean(all_wig3[,4])
            m_wig3_nequal=which(m_wig3_ttest[,3]<=0.05)
            m_wig3_more=m_wig3_nequal[which(clust_wig[m_wig3_nequal[],3]>m_wig3_mean)]
            m_wig3_less=m_wig3_nequal[which(clust_wig[m_wig3_nequal[],3]<=m_wig3_mean)]
            m_ttest3_result[,1]=as.data.frame(as.matrix(m_wig3_ttest[,3]))
            m_ttest3_result[,1]=as.numeric(m_ttest3_result[,1])
            m_ttest3_result[m_wig3_more[],2]="more"
            m_ttest3_result[m_wig3_less[],2]="less"
            m_ttest3_result[m_wig3_equal[],2]="equal"
            if(is.null(m_ttest_result))
            {
              m_ttest_result=m_ttest3_result
              
            }else
            {
              m_ttest_result=cbind(m_ttest_result,m_ttest3_result)
              
            }
            wig3_test=rbind(cbind(na.omit(clust_wig[,3]),1),cbind(na.omit(all_wig3[,4]),2))
            wig3_test=as.data.frame(wig3_test)
            colnames(wig3_test)=c("wig3_value","group")
            rownames(wig3_test)=c(1:(dim(wig3_test)[1]))
            wig3_test$group=as.factor(wig3_test$group)
            wig3_kruskal=kruskal.test(wig3_value~group, data=wig3_test)
            wig3_kruskalmc=kruskalmc(wig3_value~group, data=wig3_test, probs=0.05)
            wig3_mult <- oneway_test(wig3_value~group, data=wig3_test,
                                     ytrafo = function(data) trafo(data, numeric_trafo = rank),
                                     xtrafo = function(data) trafo(data, factor_trafo = function(x)
                                       model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                                     teststat = "max", distribution = approximate(B = 90000))
            wig3_pvalue=pvalue(wig3_mult, method = "single-step")
            write("the statistic test between BED WIG3 and ALL WIG3 :",file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("test name : Kruskal-Wallis rank sum test",file=logfile,append=TRUE)
            write(paste("Kruskal-Wallis chi-squared : ",wig3_kruskal$statistic,sep = ""),file=logfile,append=TRUE)
            write(paste("Kruskal-Wallis df : ",wig3_kruskal$parameter,sep = ""),file=logfile,append=TRUE)
            write(paste("Kruskal-Wallis p value : ",wig3_kruskal$p.value,sep = ""),file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("test name : Multiple comparison test after Kruskal-Wallis",file=logfile,append=TRUE)
            write(paste("significance level : ",wig3_kruskalmc$signif.level,sep = ""),file=logfile,append=TRUE)
            write(paste("observed difference  : ",wig3_kruskalmc$dif.com$obs.dif,sep = ""),file=logfile,append=TRUE)
            write(paste("critical difference  : ",wig3_kruskalmc$dif.com$critical.dif,sep = ""),file=logfile,append=TRUE)
            write(paste("exist difference  : ",wig3_kruskalmc$dif.com$difference,sep = ""),file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
            write("",file=logfile,append=TRUE)
          }
          
          
          
          
          clust_group=cbind(na.omit(cbind(chrBed,clust_wig,m_ttest_result)),cluster.id)
          
          write.csv(clust_group,paste(genomeName,"/",chrName,"_cluster.csv",sep=""),row.names = FALSE)
          clust_heatmap=NULL
          for( ii in 1:clust_k)
          {
            clust_heatmap=rbind(clust_heatmap,clust_group[which(clust_group[,"cluster.id"]==ii),])
          }
          rownames(clust_heatmap)=NULL
          if((outputpdf=="TRUE")||(outputpdf=="true"))
          {
            pdf(paste(genomeName,"/",chrName,"_cluster_heatmap.pdf",sep=""),width = 8,height = 8)
            
          }else
          {
            jpeg(paste(genomeName,"/",chrName,"_cluster_heatmap.jpeg",sep=""),width=1000,height=1000,quality = 100)
          }
          heatmap(as.matrix(clust_heatmap[,5:(5+wig_num-1)]),Rowv=NA,Colv=NA,cexCol = 1,labCol = "")
          dev.off()
          
          
          
          
          n_count=0
          for (ii in 1:chrTotSize)
          {
            if(chrBedMatrix[ii]!=0)
            {
              n_count=n_count+1
            }
          }
          random_group=matrix(data=0, nrow = groupNum, ncol = chrTotSize)
          random_result=matrix(data=0, nrow = groupNum , ncol = 3)
          
          for (ii in 1:groupNum)
          {
            tmp_site=sample(1:chrTotSize,size=n_count)
            random_group[ii,tmp_site]= 2
          }
          
          tmpCmap=chrCmap
          tmpCmap[lower.tri(tmpCmap)]=NA
          
          for(t in 1:groupNum)
          {
            
            random_result[t,1]=0
            random_result[t,2]=0
            for(ii in 1:chrTotSize)
            {
              if(random_group[t,ii]>0)
              {
                b=which(tmpCmap[ii,]>0)
                c=which((random_group[t,b])>0)
                random_result[t,1]=random_result[t,1]+length(c)
                random_result[t,2]=random_result[t,2]+length(b)
              }
            }
            
          }
          
          bed_count=0
          all_count=0
          bb_count=0
          bb_info=NULL
          bedTOall_info=NULL
          
          for(ii in 1:chrTotSize)
          {
            if(chrBedMatrix[ii]>0)
            {
              b=which(tmpCmap[ii,]>0)
              c=which((chrBedMatrix[b])>0)
              if(length(b)>0)
              {
                bedTOall_info=rbind(bedTOall_info,cbind(ii,b,t(tmpCmap[ii,b])))
              }
              if(length(c)>0)
              {
                bb_info=rbind(bb_info,cbind(ii,b[c],t(tmpCmap[ii,b[c]])))
              }
              bb_count=bb_count+length(c)
              bed_count=bed_count+length(b)
            }
          }
          
          if_test=rbind(cbind(bb_info[,3],1),cbind(bedTOall_info[,3],2))
          if_test=as.data.frame(if_test)
          colnames(if_test)=c("if_value","group")
          rownames(if_test)=c(1:(dim(if_test)[1]))
          if_test$group=as.factor(if_test$group)
          if_kruskal=kruskal.test(if_value~group, data=if_test)
          if_kruskalmc=kruskalmc(if_value~group, data=if_test, probs=0.05)
          mult <- oneway_test(if_value~group, data=if_test,
                              ytrafo = function(data) trafo(data, numeric_trafo = rank),
                              xtrafo = function(data) trafo(data, factor_trafo = function(x)
                                model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                              teststat = "max", distribution = approximate(B = 90000))
          if_pvalue=pvalue(mult, method = "single-step")
          dif_result=t.test(random_result[,1],mu=bb_count)
          
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("the statistic test of interaction frequency between b2b and b2o :",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("test name : Kruskal-Wallis rank sum test",file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis chi-squared : ",if_kruskal$statistic,sep = ""),file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis df : ",if_kruskal$parameter,sep = ""),file=logfile,append=TRUE)
          write(paste("Kruskal-Wallis p value : ",if_kruskal$p.value,sep = ""),file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("test name : Multiple comparison test after Kruskal-Wallis",file=logfile,append=TRUE)
          write(paste("significance level : ",if_kruskalmc$signif.level,sep = ""),file=logfile,append=TRUE)
          write(paste("observed difference  : ",if_kruskalmc$dif.com$obs.dif,sep = ""),file=logfile,append=TRUE)
          write(paste("critical difference  : ",if_kruskalmc$dif.com$critical.dif,sep = ""),file=logfile,append=TRUE)
          write(paste("exist difference  : ",if_kruskalmc$dif.com$difference,sep = ""),file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("the statistic test of interaction number between b2b and o2o :",file=logfile,append=TRUE)
          write("",file=logfile,append=TRUE)
          write("test name : t-test",file=logfile,append=TRUE)
          write(paste("numbers of random group : ",groupNum,sep = ""),file=logfile,append=TRUE)
          write(paste("95 percent confidence interval of random group : ",dif_result$conf.int[1]," ~ ",dif_result$conf.int[2],sep = ""),file=logfile,append=TRUE)
          write(paste("numbers of b2b : ",bb_count,sep = ""),file=logfile,append=TRUE)
          write(paste("t test p value : ",dif_result$p.value,sep = ""),file=logfile,append=TRUE)
          if(dif_result$p.value<0.05)
          {
            write("exist difference  : TRUE",file=logfile,append=TRUE)
          }else
          {
            write("exist difference  : FALSE",file=logfile,append=TRUE)
          }
        }
      }
    }
  }
}





