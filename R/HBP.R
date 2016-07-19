generate_enzyme_file<-function(enzyme="HindIII",enzymesite="AAGCTT",chrom_file="chrom_hg19.sizes",enzymedir="annotation",enzymeoverhangs5=1,genomeName="hg19",resolution=100)
{
  chrom_info=read.table(chrom_file,fill=TRUE, stringsAsFactors=FALSE)
  #all_genomedb=read.csv("all_genome_db.csv",header=FALSE)
  enzymesitesfile=paste(enzyme,"_resfrag_",genomeName,".bed",sep="")
  enzyme_dir=list.files(path=enzymedir,full.names=F,pattern=enzymesitesfile)
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
    if(!is.null(r_genomedb))
    {
      if(dbrequire==FALSE)
      {
        print(paste("try to install the R package ",r_genomedb,sep=""))
        source("http://www.bioconductor.org/biocLite.R")
        biocLite(r_genomedb)
      }else
      {
        all_chr <- chrom_info[,1]
        resFrag <- getRestrictionFragmentsPerChromosome(resSite=enzymesite, chromosomes=all_chr, overhangs5=enzymeoverhangs5, genomePack=r_genomedb)
        allRF <- do.call("c",resFrag)
        names(allRF) <- unlist(sapply(resFrag, function(x){paste0("HIC_", seqlevels(x), "_", 1:length(x))}))
        export(allRF, format="bed", con=paste(enzymedir,"/",enzymesitesfile,sep=""))
      }
    }else
    {
      print(paste("can not find the genome file",sep=""))
    }
  }else
  {
    print(paste("enzymesites file ",enzymesitesfile," is already existed, skip to generate it"))
  }
}

run_hicpro<-function(hicpro_path="HiC-Pro",inputfile="rawdata",configfile="config-hicpro",outdir="hg19")
{
  hiccmd=paste(hicpro_path," -i ",inputfile," -o ",outdir," -c ",configfile,sep="")
  print(hiccmd)
  system(hiccmd)
}


generate_matrix<-function(all_hic_file,all_bed_file,outputpdf=FALSE,matrix_dir="dm3",resolution=5,chrom_file="chrom_dm3.sizes")
{
  all_bed_data=read.table(file=all_bed_file, fill=TRUE, stringsAsFactors=FALSE)
  all_hic_data=read.table(file=all_hic_file,fill=TRUE, stringsAsFactors=FALSE)
  chrom_info=read.table(chrom_file,fill=TRUE, stringsAsFactors=FALSE)
  chr_num=dim(chrom_info)[1]
  inter_interaction=NULL
  if(file.exists(matrix_dir)==FALSE)
  {
    dir.create(matrix_dir)
  }
  for(tttt in 1:chr_num)
  {
    chrname=chrom_info[tttt,1]
    chr_bed=NULL
    chr_bed=rbind(chr_bed,all_bed_data[all_bed_data[,1]==chrname,])
    chr_max_length=dim(chr_bed)[1]
    chr_hic=NULL
    chr_hic=all_hic_data[all_hic_data[,1]>=chr_bed[1,4],]
    chr_hic=chr_hic[chr_hic[,1]<=chr_bed[chr_max_length,4],]
    inter_interaction=rbind(inter_interaction,chr_hic[chr_hic[,2]<chr_bed[1,4],])
    inter_interaction=rbind(inter_interaction,chr_hic[chr_hic[,2]>chr_bed[chr_max_length,4],])
    chr_hic=chr_hic[chr_hic[,2]>=chr_bed[1,4],]
    chr_hic=chr_hic[chr_hic[,2]<=chr_bed[chr_max_length,4],]
    chr_hic_data=matrix(data=0, nrow = chr_max_length, ncol = chr_max_length)
    for(i in 1:dim(chr_hic)[1])
    {
      chr_hic_data[(chr_hic[i,1]-chr_bed[1,4]),(chr_hic[i,2]-chr_bed[1,4])]=chr_hic[i,3]
      chr_hic_data[(chr_hic[i,2]-chr_bed[1,4]),(chr_hic[i,1]-chr_bed[1,4])]=chr_hic[i,3]
    }
    tmpfilename=paste(matrix_dir,"/",chrname,".matrix",sep="")
    write.table(chr_hic_data,tmpfilename,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

    if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
    {
      pdf(paste(matrix_dir,"/",chrname,"_heatmap.pdf",sep=""))

    }else
    {
      jpeg(paste(matrix_dir,"/",chrname,"_heatmap.jpeg",sep=""),width=1000,height=1000,quality = 100)
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
    mt=date()
    print(paste(mt," ",chrname," finish!",sep=""))
  }
  tmpfilename=paste(matrix_dir,"/inter_chrom.iam",sep="")
  write.table(inter_interaction,tmpfilename,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}



if_distribution_analysis<-function(all_hic_file,all_bed_file,bedFile,inter_chromfile=NULL,groupNum=50,random_analysis=TRUE,threshold_percent=0.005,if_bin_number=20,outputpdf=FALSE,matrix_dir="dm3",resolution=5,chrom_file="chrom_dm3.sizes",slide_window=FALSE)
{
  all_bed_data=read.table(file=all_bed_file, fill=TRUE, stringsAsFactors=FALSE)
  all_hic_data=read.table(file=all_hic_file,fill=TRUE, stringsAsFactors=FALSE)
  colnames(all_hic_data)=c("loc1","loc2","value")
  colnames(all_bed_data)=c("chrom","start","end","id")
  chrom_info=read.table(chrom_file,fill=TRUE, stringsAsFactors=FALSE)
  chr_num=dim(chrom_info)[1]
  m_bed_data=read.table(file=bedFile,fill=TRUE, stringsAsFactors=FALSE)

  if(is.null(inter_chromfile))
  {
    inter_interaction=read.table(file=paste(matrix_dir,"/inter_chrom.iam",sep=""), fill=TRUE, stringsAsFactors=FALSE)
  }else
  {
    inter_interaction=read.table(file=inter_chromfile, fill=TRUE, stringsAsFactors=FALSE)
  }
  order_hic_data<-all_hic_data[with(all_hic_data, order(value,decreasing = TRUE)), ]
  sum_counts=sum(order_hic_data[,3])
  th_1=mean(inter_interaction[,3])
  if(threshold_percent==0)
  {
    th_2=order_hic_data[1,3]
  }else
  {
    th_2=order_hic_data[ceiling(dim(order_hic_data)[1]*threshold_percent),3]
  }
  if_bin_size=(th_2-th_1)/if_bin_number
  order_hic_data=NULL
  rm(order_hic_data)
  gc()
  all_hic_data=all_hic_data[which(all_hic_data[,3]>=th_1),]
  all_hic_data=all_hic_data[which(all_hic_data[,3]<=th_2),]
  all_hic_data=all_hic_data[which(all_hic_data[,1]!=all_hic_data[,2]),]

  m_bed_loc=NULL
  for(nn in 1:dim(m_bed_data)[1])
  {
    tmp_bed_loc=all_bed_data[which(all_bed_data[,1]==m_bed_data[nn,1]),]
    tmp_bed_loc=tmp_bed_loc[which(tmp_bed_loc[,2]>(m_bed_data[nn,2]-resolution*1000)),]
    tmp_bed_loc=tmp_bed_loc[which(tmp_bed_loc[,3]<(m_bed_data[nn,3]+resolution*1000)),]
    m_bed_loc=rbind(m_bed_loc,tmp_bed_loc)
  }
  m_bed_loc=m_bed_loc[!duplicated(m_bed_loc),]

  all_coverage=matrix(data=0, nrow = chr_num, ncol = if_bin_number)
  all_actionnum=matrix(data=0,nrow=chr_num,ncol=if_bin_number)
  all_bb_ia=matrix(data=0,nrow=chr_num,ncol = if_bin_number)

  random_coverage=array(0,c(chr_num,if_bin_number,groupNum))
  random_bb_ia=array(0,c(chr_num,if_bin_number,groupNum))
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
    order_chr_hic<-chr_hic[with(chr_hic, order(value,decreasing = TRUE)), ]
    chr_bed_loc=as.vector(m_bed_loc[which(m_bed_loc[,1]==chrname),4])
    for(mm in 1:if_bin_number)
    {
      if_bin_start=th_1+(mm-1)*if_bin_size
      if_bin_end=th_1+mm*if_bin_size
      if_bin_hic=chr_hic[which(chr_hic[,3]>=if_bin_start),]
      if_bin_hic=if_bin_hic[which(if_bin_hic[,3]<if_bin_end),]
      hic_overlap=rbind(subset(if_bin_hic,loc1%in%chr_bed_loc),subset(if_bin_hic,loc2%in%chr_bed_loc))
      hic_overlap=hic_overlap[!duplicated(hic_overlap),]
      hic_bb_ia_1=subset(if_bin_hic,loc1%in%chr_bed_loc)
      hic_bb_ia_1=subset(hic_bb_ia_1,loc2%in%chr_bed_loc)
      hic_bb_ia_2=subset(if_bin_hic,loc2%in%chr_bed_loc)
      hic_bb_ia_2=subset(hic_bb_ia_2,loc1%in%chr_bed_loc)
      hic_bb_ia=rbind(hic_bb_ia_1,hic_bb_ia_2)
      hic_bb_ia=hic_bb_ia[!duplicated(hic_bb_ia),]
      all_actionnum[tttt,mm]=dim(if_bin_hic)[1]
      all_coverage[tttt,mm]=dim(hic_overlap)[1]
      all_bb_ia[tttt,mm]=dim(hic_bb_ia)[1]
      if((random_analysis==TRUE)||(random_analysis=="TRUE")||(random_analysis=="true"))
      {
        min_loc=min(chr_bed[,4])
        max_loc=max(chr_bed[,4])
        loc_number=length(chr_bed_loc)
        for(bb in 1:groupNum)
        {
          random_site=sample(min_loc:max_loc,size=loc_number)
          random_site=sort(random_site)
          random_overlap=rbind(subset(if_bin_hic,loc1%in%random_site),subset(if_bin_hic,loc2%in%random_site))
          random_overlap=random_overlap[!duplicated(random_overlap),]
          n_bb_ia_1=subset(if_bin_hic,loc1%in%random_site)
          n_bb_ia_1=subset(n_bb_ia_1,loc2%in%random_site)
          n_bb_ia_2=subset(if_bin_hic,loc2%in%random_site)
          n_bb_ia_2=subset(n_bb_ia_2,loc1%in%random_site)
          n_bb_ia=rbind(n_bb_ia_1,n_bb_ia_2)
          n_bb_ia=n_bb_ia[!duplicated(n_bb_ia),]
          random_coverage[tttt,mm,bb]=dim(random_overlap)[1]
          random_bb_ia[tttt,mm,bb]=dim(random_overlap)[1]
        }
      }
    }
  }

  coverage_percent=matrix(data=0,nrow=chr_num,ncol=if_bin_number)
  bb_ia_percent=matrix(data=0,nrow=chr_num,ncol=if_bin_number)
  random_cp=array(0,c(chr_num,if_bin_number,groupNum))
  random_bb=array(0,c(chr_num,if_bin_number,groupNum))
  sum_cp=matrix(data=0,nrow=1,ncol=if_bin_number)
  random_sum_cp=matrix(data=0,nrow=groupNum,ncol=if_bin_number)
  for(i1 in 1:if_bin_number)
  {
    if(sum(all_actionnum[,i1])!=0)
    {
      sum_cp[1,i1]=(sum(all_coverage[,i1]))/(sum(all_actionnum[,i1]))
    }
    if(sum(all_actionnum[,i1])==0)
    {
      sum_cp[1,i1]=0
    }
    if((random_analysis==TRUE)||(random_analysis=="TRUE")||(random_analysis=="true"))
    {
      for(i2 in 1:groupNum)
      {
        if(sum(all_actionnum[,i1])!=0)
        {
          random_sum_cp[i2,i1]=(sum(random_coverage[,i1,i2]))/(sum(all_actionnum[,i1]))
        }
        if(sum(all_actionnum[,i1])==0)
        {
          random_sum_cp[i2,i1]=0
        }
      }
    }
  }

  for(i1 in 1:chr_num)
  {
    for(i2 in 1:if_bin_number)
    {
      if(all_actionnum[i1,i2]!=0)
      {
        coverage_percent[i1,i2]=all_coverage[i1,i2]/all_actionnum[i1,i2]
        bb_ia_percent[i1,i2]=all_bb_ia[i1,i2]/all_actionnum[i1,i2]
      }
      if(all_actionnum[i1,i2]==0)
      {
        coverage_percent[i1,i2]=0
        bb_ia_percent[i1,i2]=0
      }
      #coverage_percent[i1,i2]=all_coverage[i1,i2]/all_actionnum[i1,i2]
      if((random_analysis==TRUE)||(random_analysis=="TRUE")||(random_analysis=="true"))
      {
        for(i3 in 1:groupNum)
        {
          if(all_actionnum[i1,i2]!=0)
          {
            random_cp[i1,i2,i3]=random_coverage[i1,i2,i3]/all_actionnum[i1,i2]
            random_bb[i1,i2,i3]=random_bb_ia[i1,i2,i3]/all_actionnum[i1,i2]
          }
          if(all_actionnum[i1,i2]==0)
          {
            random_cp[i1,i2,i3]=0
            random_bb[i1,i2,i3]=0
          }
        }
      }
    }
  }


  ttest_result=array(0,c(chr_num,if_bin_number,5))
  bb_ttest_result=array(0,c(chr_num,if_bin_number,5))
  sum_ttest=matrix(data=0,nrow = 5,ncol=if_bin_number)
  for(i1 in 1:if_bin_number)
  {
    tmpttest=t.test(random_sum_cp[,i1],mu=sum_cp[1,i1])
    sum_ttest[1,i1]=tmpttest$conf.int[1]
    sum_ttest[2,i1]=tmpttest$conf.int[2]
    sum_ttest[3,i1]=tmpttest$estimate
    sum_ttest[4,i1]=sum_cp[1,i1]
    sum_ttest[5,i1]=tmpttest$p.value
  }
  for(i1 in 1:chr_num)
  {
    for(i2 in 1:if_bin_number)
    {
      tmpttest=t.test(random_cp[i1,i2,],mu=coverage_percent[i1,i2])
      ttest_result[i1,i2,1]=tmpttest$conf.int[1]
      ttest_result[i1,i2,2]=tmpttest$conf.int[2]
      ttest_result[i1,i2,3]=tmpttest$estimate
      ttest_result[i1,i2,4]=coverage_percent[i1,i2]
      ttest_result[i1,i2,5]=tmpttest$p.value

      tmpttest=t.test(random_bb[i1,i2,],mu=bb_ia_percent[i1,i2])
      bb_ttest_result[i1,i2,1]=tmpttest$conf.int[1]
      bb_ttest_result[i1,i2,2]=tmpttest$conf.int[2]
      bb_ttest_result[i1,i2,3]=tmpttest$estimate
      bb_ttest_result[i1,i2,4]=coverage_percent[i1,i2]
      bb_ttest_result[i1,i2,5]=tmpttest$p.value
    }
  }

  m_col=c("#e6f5c9","#ffffbf")

  if((slide_window==TRUE)||(slide_window=="TRUE")||(slide_window=="true"))
  {
    for(i1 in 1:chr_num)
    {
      chrname=chrom_info[i1,1]
      if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
      {
        pdf(paste(matrix_dir,"/",chrname,"_if_dis_slide.pdf",sep=""),width = 8,height = 8)

      }else
      {
        jpeg(paste(matrix_dir,"/",chrname,"_if_dis_slide.jpeg",sep=""),width=1000,height=1000,quality = 100)
      }
      boxplot(t(random_cp[i1,,]),col=m_col,xlab="interaction frequency distribution")
      #lines(ttest_result[i1,,4],lwd=3,col="red",lty=4)
      m_ly=numeric(if_bin_number)
      for(m_li in 1:(if_bin_number-2))
      {
        m_ly[m_li+1]=sum(ttest_result[i1,m_li:(m_li+2),4])/3
      }
      m_ly[1]=sum(ttest_result[i1,1:2,4])/2
      m_ly[if_bin_number]=sum(ttest_result[i1,(if_bin_number-1):if_bin_number,4])/2
      lines(m_ly,lwd=3,col="red",lty=4)
      dev.off()
    }
    if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
    {
      pdf(paste(matrix_dir,"/all_chrom_if_dis_slide.pdf",sep=""),width = 8,height = 8)

    }else
    {
      jpeg(paste(matrix_dir,"/all_chrom_if_dis_slide.jpeg",sep=""),width=1000,height=1000,quality = 100)
    }

    boxplot(random_sum_cp,col=m_col,xlab="interaction frequency distribution")
    #lines(sum_ttest[4,],lwd=3,col="red",lty=4)
    sum_ly=numeric(if_bin_number)
    for(sum_li in 1:(if_bin_number-2))
    {
      sum_ly[sum_li+1]=sum(sum_ttest[4,sum_li:(sum_li+2)])/3
    }
    sum_ly[1]=sum(sum_ttest[4,1:2])/2
    sum_ly[if_bin_number]=sum(sum_ttest[4,(if_bin_number-1):if_bin_number])/2
    lines(sum_ly,lwd=3,col="red",lty=4)
    dev.off()
  }else
  {
    for(i1 in 1:chr_num)
    {
      chrname=chrom_info[i1,1]
      if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
      {
        pdf(paste(matrix_dir,"/",chrname,"_if_dis.pdf",sep=""),width = 8,height = 8)

      }else
      {
        jpeg(paste(matrix_dir,"/",chrname,"_if_dis.jpeg",sep=""),width=1000,height=1000,quality = 100)
      }
      boxplot(t(random_cp[i1,,]),col=m_col,xlab="interaction frequency distribution")
      lines(ttest_result[i1,,4],lwd=3,col="red",lty=4)
      dev.off()
    }
    if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
    {
      pdf(paste(matrix_dir,"/all_chrom_if_dis.pdf",sep=""),width = 8,height = 8)

    }else
    {
      jpeg(paste(matrix_dir,"/all_chrom_if_dis.jpeg",sep=""),width=1000,height=1000,quality = 100)
    }

    boxplot(random_sum_cp,col=m_col,xlab="interaction frequency distribution")
    lines(sum_ttest[4,],lwd=3,col="red",lty=4)
    dev.off()
  }


}




network_analysis<-function(bedFile,matrix_dir="hg19",outputpdf=FALSE,chrom="all",chrstart=0,chrend=0,resolution=100,bedWindow=0,net_layout="layout.fruchterman.reingold",netplot=TRUE,NetClusterType="multileve",NetVertexSize=2,NetVertexChangeSize="degree",NetVertexLableDist=0.1,NetVertexColor="#7fbc41",NetVertexLabelCex=3,if_threshold=0)
{
  gc()
  #source("calbed.R")
  matrix_name_dir=list.files(path=matrix_dir,full.names=F,pattern=".matrix")
  matrix_full_dir=list.files(path=matrix_dir,full.names=T,pattern=".matrix")
  m_bed=load_bed(bedFile)
  chrNum=length(matrix_name_dir)
  if(chrom=="all")
  {
    for (i in 1:chrNum)
    {


      tmpNum=regexpr(".matrix",matrix_name_dir[i])
      chrName=substr(matrix_name_dir[i],1,tmpNum-1)
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
        chrBedToBedInter=find_bed_to_bed_interaction(chrCmap,chrBedMatrix,chrBedBin,chrBed,chrName,chrCmap,if_threshold)

        if((netplot==TRUE)||(netplot=="true")||(netplot=="TRUE"))
        {
          netgraph=data.frame("p1"=character(dim(chrBedToBedInter)[1]),"p2"=character(dim(chrBedToBedInter)[1]),"weight"=numeric(dim(chrBedToBedInter)[1]))
          netgraph[,1]=as.data.frame(paste(chrBedToBedInter[,2],":",(chrBedToBedInter[,3]+chrstart),"-",(chrBedToBedInter[,4]+chrstart),sep=""))
          netgraph[,2]=as.data.frame(paste(chrBedToBedInter[,8],":",(chrBedToBedInter[,9]+chrstart),"-",(chrBedToBedInter[,10]+chrstart),sep=""))
          netgraph[,3]=as.data.frame(chrBedToBedInter[,13])

          set.seed(1234)
          #pdf(paste(matrix_dir,"/",chrName,"_netplot.pdf",sep=""))

          if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
          {
            pdf(paste(matrix_dir,"/",chrName,"_netplot.pdf",sep=""),width = 8,height = 8)

          }else
          {
            jpeg(paste(matrix_dir,"/",chrName,"_netplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
          }
          #pdf("tmp_pdf.pdf",width=8,height=8)
          g = graph.data.frame(netgraph,directed = F)
          set.seed(1234)

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
          #colors<-c("#f7f4f9","#f7f4f9","#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f")
          colors<-c("#fff7f3","#fff7f3","#fde0dd","#fcc5c0","#fa9fb5","#f768a1","#dd3497","#b0017e","#7a0177","#49006a")
          #colors<-c("#fff7fb","#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#570b0","#45a8d","#23858")
          #colors<-c("#ffffcc","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#b00026")
          #colors<-c("#f7fbff","#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#8519c","#8306b")
          weight_range=range(E(g)$weight)
          E(g)$color=colors[1]
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])/10)]$color=colors[2]
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*2/10)]$color=colors[3]
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*3/10)]$color=colors[4]
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*4/10)]$color=colors[5]
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*5/10)]$color=colors[6]
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*6/10)]$color=colors[7]
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*7/10)]$color=colors[8]
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*8/10)]$color=colors[9]
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*9/10)]$color=colors[10]

          edge_width=0.05
          E(g)$width=edge_width
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])/10)]$width=2*edge_width
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*2/10)]$width=3*edge_width
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*3/10)]$width=4*edge_width
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*4/10)]$width=5*edge_width
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*5/10)]$width=6*edge_width
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*6/10)]$width=7*edge_width
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*7/10)]$width=8*edge_width
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*8/10)]$width=9*edge_width
          E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*9/10)]$width=10*edge_width

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

          if(net_layout=="layout.fruchterman.reingold")
          {

            if(NetClusterType=="NULL")
            {
              plot(g,layout=layout.fruchterman.reingold, vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              netcsv[1:(length(netdegree)),9]="NULL"
            }
            if(NetClusterType=="edgeBetweenness")
            {
              system.time(ec <- edge.betweenness.community(g))
              print(modularity(ec))
              netcsv[1:(length(netdegree)),9]=ec$membership
              plot(ec, g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
            }
            if(NetClusterType=="walktrap")
            {
              system.time(wc <- walktrap.community(g))
              netcsv[1:(length(netdegree)),9]=wc$membership

              print(modularity(wc))
              plot(wc , g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)

            }
            if(NetClusterType=="multileve")
            {
              system.time(mc <- multilevel.community(g, weights=NA))
              print(modularity(mc))
              plot(mc, g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              netcsv[1:(length(netdegree)),9]=mc$membership

            }
            if(NetClusterType=="labelPropagation")
            {
              system.time(lc <- label.propagation.community(g))
              print(modularity(lc))
              plot(lc , g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              netcsv[1:(length(netdegree)),9]=lc$membership

            }
          }else if(net_layout=="layout.circle")
          {
            if(NetClusterType=="NULL")
            {
              plot(g,layout=layout.circle, vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              netcsv[1:(length(netdegree)),9]="NULL"
            }
            if(NetClusterType=="edgeBetweenness")
            {
              system.time(ec <- edge.betweenness.community(g))
              print(modularity(ec))
              netcsv[1:(length(netdegree)),9]=ec$membership
              plot(ec, g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
            }
            if(NetClusterType=="walktrap")
            {
              system.time(wc <- walktrap.community(g))
              netcsv[1:(length(netdegree)),9]=wc$membership

              print(modularity(wc))
              plot(wc , g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)

            }
            if(NetClusterType=="multileve")
            {
              system.time(mc <- multilevel.community(g, weights=NA))
              print(modularity(mc))
              plot(mc, g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              netcsv[1:(length(netdegree)),9]=mc$membership

            }
            if(NetClusterType=="labelPropagation")
            {
              system.time(lc <- label.propagation.community(g))
              print(modularity(lc))
              plot(lc , g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              netcsv[1:(length(netdegree)),9]=lc$membership

            }
          }else
          {
            if(NetClusterType=="NULL")
            {
              plot(g,layout=layout.auto, vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              netcsv[1:(length(netdegree)),9]="NULL"
            }
            if(NetClusterType=="edgeBetweenness")
            {
              system.time(ec <- edge.betweenness.community(g))
              print(modularity(ec))
              netcsv[1:(length(netdegree)),9]=ec$membership
              plot(ec, g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
            }
            if(NetClusterType=="walktrap")
            {
              system.time(wc <- walktrap.community(g))
              netcsv[1:(length(netdegree)),9]=wc$membership

              print(modularity(wc))
              plot(wc , g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)

            }
            if(NetClusterType=="multileve")
            {
              system.time(mc <- multilevel.community(g, weights=NA))
              print(modularity(mc))
              plot(mc, g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              netcsv[1:(length(netdegree)),9]=mc$membership

            }
            if(NetClusterType=="labelPropagation")
            {
              system.time(lc <- label.propagation.community(g))
              print(modularity(lc))
              plot(lc , g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
              netcsv[1:(length(netdegree)),9]=lc$membership

            }
          }

          dev.off()

          write.csv(netcsv,paste(matrix_dir,"/",chrName,"_network.csv",sep=""),row.names = FALSE)
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

        if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
        {
          pdf(paste(matrix_dir,"/",chrName,"_bedplot.pdf",sep=""),width = 8,height = 8)

        }else
        {
          jpeg(paste(matrix_dir,"/",chrName,"_bedplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
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
        chrbeddensitydata=chrbeddensitydata+hmrange[1]
        colnames(chrbeddensitydata)="bed"
        chrbeddensity=ggplot(chrbeddensitydata)+geom_density(aes(x=bed))
        chrhmdensitydata=NULL
        for(iiii in 1:chrTotSize)
        {
          #pp[iiii,1]=length(which(chrCmap[,iiii]>0))
          chrhmdensitydata=c(chrhmdensitydata,which(chrCmap[,iiii]>0))
        }
        chrhmdensitydata=as.data.frame(chrhmdensitydata)
        chrhmdensitydata=chrhmdensitydata+hmrange[1]
        colnames(chrhmdensitydata)="chrom"
        chrcmapdensity=ggplot(chrhmdensitydata)+geom_density(aes(x=chrom))
        print(chrhm,vp=heatmapViewport)
        print(chrbedplot,vp=scatterViewport)
        print(chrbeddensity,vp=densityViewport)
        print(chrcmapdensity,vp=hmdensityViewport)

        dev.off()

        write.table(chrBedToBedInter,file=paste(matrix_dir,"/",chrName,"_BedToBedInter.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
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


      tmpNum=regexpr(".matrix",matrix_name_dir[i])
      chrName=substr(matrix_name_dir[i],1,tmpNum-1)
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
            chrBedToBedInter=find_bed_to_bed_interaction(chrCmap,chrBedMatrix,chrBedBin,chrBed,chrName,chrCmap,if_threshold)

            if((netplot==TRUE)||(netplot=="true")||(netplot=="TRUE"))
            {
              netgraph=data.frame("p1"=character(dim(chrBedToBedInter)[1]),"p2"=character(dim(chrBedToBedInter)[1]),"weight"=numeric(dim(chrBedToBedInter)[1]))
              netgraph[,1]=as.data.frame(paste(chrBedToBedInter[,2],":",(chrBedToBedInter[,3]+chrstart),"-",(chrBedToBedInter[,4]+chrstart),sep=""))
              netgraph[,2]=as.data.frame(paste(chrBedToBedInter[,8],":",(chrBedToBedInter[,9]+chrstart),"-",(chrBedToBedInter[,10]+chrstart),sep=""))
              netgraph[,3]=as.data.frame(chrBedToBedInter[,13])
              set.seed(1234)
              if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
              {
                pdf(paste(matrix_dir,"/",chrName,"_netplot.pdf",sep=""),width = 8,height = 8)

              }else
              {
                jpeg(paste(matrix_dir,"/",chrName,"_netplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
              }
              #pdf("tmp_pdf.pdf",width=8,height=8)
              g = graph.data.frame(netgraph,directed = F)
              set.seed(1234)

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
              #colors<-c("#f7f4f9","#f7f4f9","#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f")
              colors<-c("#fff7f3","#fff7f3","#fde0dd","#fcc5c0","#fa9fb5","#f768a1","#dd3497","#b0017e","#7a0177","#49006a")
              #colors<-c("#fff7fb","#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#570b0","#45a8d","#23858")
              #colors<-c("#ffffcc","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#b00026")
              #colors<-c("#f7fbff","#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#8519c","#8306b")
              weight_range=range(E(g)$weight)
              E(g)$color=colors[1]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])/10)]$color=colors[2]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*2/10)]$color=colors[3]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*3/10)]$color=colors[4]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*4/10)]$color=colors[5]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*5/10)]$color=colors[6]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*6/10)]$color=colors[7]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*7/10)]$color=colors[8]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*8/10)]$color=colors[9]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*9/10)]$color=colors[10]

              edge_width=0.05
              E(g)$width=edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])/10)]$width=2*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*2/10)]$width=3*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*3/10)]$width=4*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*4/10)]$width=5*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*5/10)]$width=6*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*6/10)]$width=7*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*7/10)]$width=8*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*8/10)]$width=9*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*9/10)]$width=10*edge_width

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


              if(net_layout=="layout.fruchterman.reingold")
              {

                if(NetClusterType=="NULL")
                {
                  plot(g,layout=layout.fruchterman.reingold, vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]="NULL"
                }
                if(NetClusterType=="edgeBetweenness")
                {
                  system.time(ec <- edge.betweenness.community(g))
                  print(modularity(ec))
                  netcsv[1:(length(netdegree)),9]=ec$membership
                  plot(ec, g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                }
                if(NetClusterType=="walktrap")
                {
                  system.time(wc <- walktrap.community(g))
                  netcsv[1:(length(netdegree)),9]=wc$membership

                  print(modularity(wc))
                  plot(wc , g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)

                }
                if(NetClusterType=="multileve")
                {
                  system.time(mc <- multilevel.community(g, weights=NA))
                  print(modularity(mc))
                  plot(mc, g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=mc$membership

                }
                if(NetClusterType=="labelPropagation")
                {
                  system.time(lc <- label.propagation.community(g))
                  print(modularity(lc))
                  plot(lc , g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=lc$membership

                }
              }else if(net_layout=="layout.circle")
              {
                if(NetClusterType=="NULL")
                {
                  plot(g,layout=layout.circle, vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]="NULL"
                }
                if(NetClusterType=="edgeBetweenness")
                {
                  system.time(ec <- edge.betweenness.community(g))
                  print(modularity(ec))
                  netcsv[1:(length(netdegree)),9]=ec$membership
                  plot(ec, g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                }
                if(NetClusterType=="walktrap")
                {
                  system.time(wc <- walktrap.community(g))
                  netcsv[1:(length(netdegree)),9]=wc$membership

                  print(modularity(wc))
                  plot(wc , g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)

                }
                if(NetClusterType=="multileve")
                {
                  system.time(mc <- multilevel.community(g, weights=NA))
                  print(modularity(mc))
                  plot(mc, g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=mc$membership

                }
                if(NetClusterType=="labelPropagation")
                {
                  system.time(lc <- label.propagation.community(g))
                  print(modularity(lc))
                  plot(lc , g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=lc$membership

                }
              }else
              {
                if(NetClusterType=="NULL")
                {
                  plot(g,layout=layout.auto, vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]="NULL"
                }
                if(NetClusterType=="edgeBetweenness")
                {
                  system.time(ec <- edge.betweenness.community(g))
                  print(modularity(ec))
                  netcsv[1:(length(netdegree)),9]=ec$membership
                  plot(ec, g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                }
                if(NetClusterType=="walktrap")
                {
                  system.time(wc <- walktrap.community(g))
                  netcsv[1:(length(netdegree)),9]=wc$membership

                  print(modularity(wc))
                  plot(wc , g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)

                }
                if(NetClusterType=="multileve")
                {
                  system.time(mc <- multilevel.community(g, weights=NA))
                  print(modularity(mc))
                  plot(mc, g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=mc$membership

                }
                if(NetClusterType=="labelPropagation")
                {
                  system.time(lc <- label.propagation.community(g))
                  print(modularity(lc))
                  plot(lc , g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=lc$membership

                }
              }
              dev.off()
              write.csv(netcsv,paste(matrix_dir,"/",chrName,"_network.csv",sep=""),row.names = FALSE)

            }
            bedIplot=cbind(rbind(chrBedToBedInter[,14],chrBedToBedInter[,15]),rbind(chrBedToBedInter[,15],chrBedToBedInter[,14]))
            write.table(chrBedToBedInter,file=paste(matrix_dir,"/",chrName,"_BedToBedInter.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

            hm_dim=dim(chrCmap)[1]
            chrCmap=as.matrix(chrCmap)
            hm_mean=mean(chrCmap)
            for(j in 1:hm_dim)
            {

              chrCmap[j,which(chrCmap[j,]>5*hm_mean)]=5*hm_mean
            }
            chrhmCmap=melt(chrCmap)
            print(paste("plot ",chrName,"bed picture",sep=""))
            if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
            {
              pdf(paste(matrix_dir,"/",chrName,"_bedplot.pdf",sep=""),width = 8,height = 8)

            }else
            {
              jpeg(paste(matrix_dir,"/",chrName,"_bedplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
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
            chrbeddensitydata=chrbeddensitydata+hmrange[1]

            colnames(chrbeddensitydata)="bed"
            chrbeddensity=ggplot(chrbeddensitydata)+geom_density(aes(x=bed))

            chrhmdensitydata=NULL
            for(iiii in 1:chrTotSize)
            {
              #pp[iiii,1]=length(which(chrCmap[,iiii]>0))
              chrhmdensitydata=c(chrhmdensitydata,which(chrCmap[,iiii]>0))
            }
            chrhmdensitydata=as.data.frame(chrhmdensitydata)
            chrhmdensitydata=chrhmdensitydata+hmrange[1]
            colnames(chrhmdensitydata)="chrom"


            chrcmapdensity=ggplot(chrhmdensitydata)+geom_density(aes(x=chrom))



            print(chrhm,vp=heatmapViewport)
            print(chrbedplot,vp=scatterViewport)
            print(chrbeddensity,vp=densityViewport)
            print(chrcmapdensity,vp=hmdensityViewport)

            dev.off()
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

            chrBedToBedInter=find_bed_to_bed_interaction(chrCmap,chrBedMatrix,chrBedBin,chrBed,chrName,chrCmap,if_threshold)

            if((netplot==TRUE)||(netplot=="true")||(netplot=="TRUE"))
            {
              netgraph=data.frame("p1"=character(dim(chrBedToBedInter)[1]),"p2"=character(dim(chrBedToBedInter)[1]),"weight"=numeric(dim(chrBedToBedInter)[1]))
              netgraph[,1]=as.data.frame(paste(chrBedToBedInter[,2],":",(chrBedToBedInter[,3]+chrstart),"-",(chrBedToBedInter[,4]+chrstart),sep=""))
              netgraph[,2]=as.data.frame(paste(chrBedToBedInter[,8],":",(chrBedToBedInter[,9]+chrstart),"-",(chrBedToBedInter[,10]+chrstart),sep=""))
              netgraph[,3]=as.data.frame(chrBedToBedInter[,13])

              set.seed(1234)

              if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
              {
                pdf(paste(matrix_dir,"/",chrName,"_netplot.pdf",sep=""),width = 8,height = 8)

              }else
              {
                jpeg(paste(matrix_dir,"/",chrName,"_netplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
              }
              #pdf("tmp_pdf.pdf",width=8,height=8)
              g = graph.data.frame(netgraph,directed = F)
              set.seed(1234)

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
              #colors<-c("#f7f4f9","#f7f4f9","#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f")
              colors<-c("#fff7f3","#fff7f3","#fde0dd","#fcc5c0","#fa9fb5","#f768a1","#dd3497","#b0017e","#7a0177","#49006a")
              #colors<-c("#fff7fb","#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#570b0","#45a8d","#23858")
              #colors<-c("#ffffcc","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#b00026")
              #colors<-c("#f7fbff","#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#8519c","#8306b")
              weight_range=range(E(g)$weight)
              E(g)$color=colors[1]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])/10)]$color=colors[2]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*2/10)]$color=colors[3]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*3/10)]$color=colors[4]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*4/10)]$color=colors[5]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*5/10)]$color=colors[6]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*6/10)]$color=colors[7]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*7/10)]$color=colors[8]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*8/10)]$color=colors[9]
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*9/10)]$color=colors[10]

              edge_width=0.05
              E(g)$width=edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])/10)]$width=2*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*2/10)]$width=3*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*3/10)]$width=4*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*4/10)]$width=5*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*5/10)]$width=6*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*6/10)]$width=7*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*7/10)]$width=8*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*8/10)]$width=9*edge_width
              E(g)[weight>=(weight_range[1]+(weight_range[2]-weight_range[1])*9/10)]$width=10*edge_width

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

              if(net_layout=="layout.fruchterman.reingold")
              {

                if(NetClusterType=="NULL")
                {
                  plot(g,layout=layout.fruchterman.reingold, vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]="NULL"
                }
                if(NetClusterType=="edgeBetweenness")
                {
                  system.time(ec <- edge.betweenness.community(g))
                  print(modularity(ec))
                  netcsv[1:(length(netdegree)),9]=ec$membership
                  plot(ec, g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                }
                if(NetClusterType=="walktrap")
                {
                  system.time(wc <- walktrap.community(g))
                  netcsv[1:(length(netdegree)),9]=wc$membership

                  print(modularity(wc))
                  plot(wc , g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)

                }
                if(NetClusterType=="multileve")
                {
                  system.time(mc <- multilevel.community(g, weights=NA))
                  print(modularity(mc))
                  plot(mc, g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=mc$membership

                }
                if(NetClusterType=="labelPropagation")
                {
                  system.time(lc <- label.propagation.community(g))
                  print(modularity(lc))
                  plot(lc , g,layout=layout.fruchterman.reingold,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=lc$membership

                }
              }else if(net_layout=="layout.circle")
              {
                if(NetClusterType=="NULL")
                {
                  plot(g,layout=layout.circle, vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]="NULL"
                }
                if(NetClusterType=="edgeBetweenness")
                {
                  system.time(ec <- edge.betweenness.community(g))
                  print(modularity(ec))
                  netcsv[1:(length(netdegree)),9]=ec$membership
                  plot(ec, g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                }
                if(NetClusterType=="walktrap")
                {
                  system.time(wc <- walktrap.community(g))
                  netcsv[1:(length(netdegree)),9]=wc$membership

                  print(modularity(wc))
                  plot(wc , g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)

                }
                if(NetClusterType=="multileve")
                {
                  system.time(mc <- multilevel.community(g, weights=NA))
                  print(modularity(mc))
                  plot(mc, g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=mc$membership

                }
                if(NetClusterType=="labelPropagation")
                {
                  system.time(lc <- label.propagation.community(g))
                  print(modularity(lc))
                  plot(lc , g,layout=layout.circle,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=lc$membership

                }
              }else
              {
                if(NetClusterType=="NULL")
                {
                  plot(g,layout=layout.auto, vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]="NULL"
                }
                if(NetClusterType=="edgeBetweenness")
                {
                  system.time(ec <- edge.betweenness.community(g))
                  print(modularity(ec))
                  netcsv[1:(length(netdegree)),9]=ec$membership
                  plot(ec, g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                }
                if(NetClusterType=="walktrap")
                {
                  system.time(wc <- walktrap.community(g))
                  netcsv[1:(length(netdegree)),9]=wc$membership

                  print(modularity(wc))
                  plot(wc , g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)

                }
                if(NetClusterType=="multileve")
                {
                  system.time(mc <- multilevel.community(g, weights=NA))
                  print(modularity(mc))
                  plot(mc, g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=mc$membership

                }
                if(NetClusterType=="labelPropagation")
                {
                  system.time(lc <- label.propagation.community(g))
                  print(modularity(lc))
                  plot(lc , g,layout=layout.auto,vertex.label.dist=NetVertexLableDist, vertex.color=NetVertexColor, edge.arrow.size=0.05,vertex.label.cex=NetVertexLabelCex)
                  netcsv[1:(length(netdegree)),9]=lc$membership

                }
              }
              dev.off()
              write.csv(netcsv,paste(matrix_dir,"/",chrName,"_network.csv",sep=""),row.names = FALSE)

            }


            bedIplot=cbind(rbind(chrBedToBedInter[,14],chrBedToBedInter[,15]),rbind(chrBedToBedInter[,15],chrBedToBedInter[,14]))
            write.table(chrBedToBedInter,file=paste(matrix_dir,"/",chrName,"_BedToBedInter.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

            hm_dim=dim(chrCmap)[1]
            chrCmap=as.matrix(chrCmap)
            hm_mean=mean(chrCmap)
            for(j in 1:hm_dim)
            {

              chrCmap[j,which(chrCmap[j,]>5*hm_mean)]=5*hm_mean
            }
            chrhmCmap=melt(chrCmap)
            print(paste("plot ",chrName,"bed picture",sep=""))
            if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
            {
              pdf(paste(matrix_dir,"/",chrName,"_bedplot.pdf",sep=""),width = 8,height = 8)

            }else
            {
              jpeg(paste(matrix_dir,"/",chrName,"_bedplot.jpeg",sep=""),width=1000,height=1000,quality = 100)
            }

            grid.newpage()
            heatmapViewport <- viewport(height=0.5, width=0.5, x=0.25,y=0.5)
            scatterViewport <- viewport(height=0.5, width=0.5, x=0.75,y=0.5)
            densityViewport <- viewport(height=0.25,width=0.5, x=0.75,y=0.125)
            hmdensityViewport <- viewport(height=0.25,width=0.5,x=0.25,y=0.125)
            jit=position_jitter(width=0.5)
            hmrange=range(chrhmCmap[,1])
            #print(bedIplot)
            #print(hmrange)
            bedIplot=bedIplot+hmrange[1]
            chrhm = ggplot(chrhmCmap, aes(x=Var1, y=Var2, fill=value))+scale_y_discrete(breaks=seq(0, 10, 5))+xlab('chrom')+ylab("chrom")+scale_fill_gradient(low='white', high='red')+geom_tile()+guides(fill=FALSE)
            chrbedplot=qplot(bedIplot[1,],bedIplot[2,],alpha=I(1/10),size=I(1))+xlab('chrom')+ylab("chrom")+geom_jitter(position=jit,colour="black",alpha=1/100)

            chrbeddensitydata=as.data.frame(c(chrBedToBedInter[,14],chrBedToBedInter[,15]))
            chrbeddensitydata=chrbeddensitydata+hmrange[1]
            colnames(chrbeddensitydata)="bed"
            chrbeddensity=ggplot(chrbeddensitydata)+geom_density(aes(x=bed))

            chrhmdensitydata=NULL
            for(iiii in 1:chrTotSize)
            {
              #pp[iiii,1]=length(which(chrCmap[,iiii]>0))
              chrhmdensitydata=c(chrhmdensitydata,which(chrCmap[,iiii]>0))
            }
            chrhmdensitydata=as.data.frame(chrhmdensitydata)
            chrhmdensitydata=chrhmdensitydata+hmrange[1]
            colnames(chrhmdensitydata)="chrom"


            chrcmapdensity=ggplot(chrhmdensitydata)+geom_density(aes(x=chrom))



            print(chrhm,vp=heatmapViewport)
            print(chrbedplot,vp=scatterViewport)
            print(chrbeddensity,vp=densityViewport)
            print(chrcmapdensity,vp=hmdensityViewport)

            dev.off()


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


circos_plot<-function(bedFile,wig_dir="wig",matrix_dir="hg19",bedWindow=0,outputpdf=TRUE,chrom="all",chrstart=0,chrend=0,resolution=100,circosLineWidth=0.01,if_threshold=0,circosLinecolor="ReadCounts",circosTrackWidth=40)
{
  #source("calbed.R")
  gc()
  matrix_name_dir=list.files(path=matrix_dir,full.names=F,pattern=".matrix")
  matrix_full_dir=list.files(path=matrix_dir,full.names=T,pattern=".matrix")
  all_wig_file=list.files(path=wig_dir,full.names=T,pattern=".wig")
  all_wig_name=list.files(path=wig_dir,full.names = FALSE,pattern = ".wig")
  all_wig_num=length(all_wig_file)
  m_bed=load_bed(bedFile)
  chrNum=length(matrix_name_dir)
  options(stringsAsFactors = FALSE);

  genomeFrame=data.frame("seg.name"=character(0),"seg.start"=numeric(0),"seg.end"=numeric(0),"the.v"=character(0),"NO"=character(0),stringsAsFactors=FALSE)
  genomeFrameNum=1
  print("generate genome frame")
  for (i in 1:chrNum)
  {
    chrCmap=read.table(file=matrix_full_dir[i], fill=TRUE, stringsAsFactors=FALSE)
    chrTotSize=dim(chrCmap)[1]
    tmpNum=regexpr(".matrix",matrix_name_dir[i])
    chrName=substr(matrix_name_dir[i],1,tmpNum-1)

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
    print(paste(chrName," frame generate finished!",sep = ""))


    rm(chrCmap)
    gc()
  }

  if(chrom=="all")
  {
    for (i in 1:chrNum)
    {

      tmpNum=regexpr(".matrix",matrix_name_dir[i])
      chrName=substr(matrix_name_dir[i],1,tmpNum-1)
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
        print(paste("plot ",chrName,"circos picture",sep=""))

        chrCircosMapping=calculate_omiccircos_data(chrCmap,chrBedMatrix,chrBedBin,chrBed,chrName,chrCmap,if_threshold)

        chrCircosMapping[,1]=chrNo
        chrCircosMapping[,4]=chrNo
        chrCircosDb=segAnglePo(genomeFrame, seg=chrName)
        seg.num<-length(unique(genomeFrame[,1]))
        colors<-rainbow(seg.num, alpha=0.5)


        if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
        {
          pdf(paste(matrix_dir,"/",chrName,"_circos.pdf",sep=""),width = 8,height = 8)

        }else
        {
          jpeg(paste(matrix_dir,"/",chrName,"_circos.jpeg",sep=""),width=1000,height=1000,quality = 100)
        }
        #pdf("tmp_circos.pdf",width = 8,height = 8)
        options(stringsAsFactors = FALSE);

        par(mar=c(2, 2, 2, 2));
        plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
        if(circosLinecolor=="rainbow")
        {
          colors<-rainbow(seg.num, alpha=0.5)
          circosLinecolor=colors
          circos(R=340, cir=chrCircosDb, W=40,  mapping=chrCircosMapping, type="link2", lwd=circosLineWidth,col=colors);
        }else if(circosLinecolor=="ReadCounts")
        {
          #colors<-c("#f7f4f9","#f7f4f9","#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f")
          colors<-c("#fff7f3","#fff7f3","#fde0dd","#fcc5c0","#fa9fb5","#f768a1","#dd3497","#b0017e","#7a0177","#49006a")
          #colors<-c("#fff7fb","#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")
          #colors<-c("#ffffcc","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#b00026")
          #colors<-c("#f7fbff","#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b")
          chrCircosMapping[,7]=as.numeric(chrCircosMapping[,7])
          chrCircosMapping<-chrCircosMapping[with(chrCircosMapping, order(read_count,decreasing = FALSE)), ]
          ccm_dim=dim(chrCircosMapping)[1]
          for(ccm in 1:10)
          {
            ccm_start=1+ceiling((ccm_dim/10)*(ccm-1))
            ccm_end=ceiling((ccm_dim/10)*ccm)
            circos(R=340, cir=chrCircosDb, W=40,  mapping=chrCircosMapping[ccm_start:ccm_end,], type="link2", lwd=circosLineWidth*ccm,col=colors[ccm]);
          }
        }else
        {
          circos(R=340, cir=chrCircosDb, W=40,  mapping=chrCircosMapping, type="link2", lwd=circosLineWidth,col=circosLinecolor);
        }

        all_wig_hm=NULL
        all_bed_wig_hm=NULL

        if(all_wig_num>0)
        {
          for(ww in 1:all_wig_num)
          {
            m_wig=load_wig(all_wig_file[ww],resolution*1000,chrName,chrTotSize,chrstart,chrend)


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

              if(is.null(all_bed_wig_hm))
              {
                all_bed_wig_hm=m_wig_bed
              }else
              {
                all_bed_wig_hm=cbind(all_bed_wig_hm,m_wig_bed[,3])
              }
              if(is.null(all_wig_hm))
              {
                all_wig_hm=m_wig
                all_wig_hm[,3]=(all_wig_hm[,3]-range(all_wig_hm[,3])[1])/(range(all_wig_hm[,3])[2]-range(all_wig_hm[,3])[1])
                tmpname=colnames(all_wig_hm)
                colnames(all_wig_hm)=c(tmpname[1:2],all_wig_name[ww])
              }else
              {
                tmpname=colnames(all_wig_hm)
                all_wig_hm=cbind(all_wig_hm,((m_wig[,3]-range(m_wig[,3])[1])/(range(m_wig[,3])[2]-range(m_wig[,3])[1])))
                colnames(all_wig_hm)=c(tmpname,all_wig_name[ww])

              }

            }
          }
        }


        circos(R=360, cir=chrCircosDb, W=circosTrackWidth, mapping=all_wig_hm,  col.v=3,  type="heatmap2",  col.bar=TRUE, lwd=0.1, col="blue")
        circos(R=(360+circosTrackWidth), cir=chrCircosDb, W=1,   type="chr", print.chr.lab=FALSE, scale=TRUE)
        text(730,820,"track name Outside-to-inside",cex=0.55,family="mono")
        for(namei in 1:all_wig_num)
        {
          text(740,820-namei*20,all_wig_name[namei],cex=0.55,family="mono")
        }
        dev.off() # text(400,820,family="mono",wigFile3,cex=0.7)
        print(paste(chrName," circos plot finish",sep = ""))


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

      tmpNum=regexpr(".matrix",matrix_name_dir[i])
      chrName=substr(matrix_name_dir[i],1,tmpNum-1)
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



          chrCircosMapping=calculate_omiccircos_data(chrCmap,chrBedMatrix,chrBedBin,chrBed,chrName,chrCmap,if_threshold)
          chrCircosMapping[,1]=chrNo
          chrCircosMapping[,4]=chrNo
          chrCircosDb=segAnglePo(genomeFrame, seg=chrName)
          seg.num<-length(unique(genomeFrame[,1]))


          if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
          {
            pdf(paste(matrix_dir,"/",chrName,"_circos.pdf",sep=""),width = 8,height = 8)

          }else
          {
            jpeg(paste(matrix_dir,"/",chrName,"_circos.jpeg",sep=""),width=1000,height=1000,quality = 100)
          }
          #pdf("tmp_circos.pdf",width = 8,height = 8)
          options(stringsAsFactors = FALSE);

          par(mar=c(2, 2, 2, 2));
          plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
          if(circosLinecolor=="rainbow")
          {
            colors<-rainbow(seg.num, alpha=0.5)
            circosLinecolor=colors
            circos(R=340, cir=chrCircosDb, W=40,  mapping=chrCircosMapping, type="link2", lwd=circosLineWidth,col=colors);
          }else if(circosLinecolor=="ReadCounts")
          {
            #colors<-c("#f7f4f9","#f7f4f9","#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f")
            colors<-c("#fff7f3","#fff7f3","#fde0dd","#fcc5c0","#fa9fb5","#f768a1","#dd3497","#b0017e","#7a0177","#49006a")
            #colors<-c("#fff7fb","#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")
            #colors<-c("#ffffcc","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#b00026")
            #colors<-c("#f7fbff","#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b")
            chrCircosMapping[,7]=as.numeric(chrCircosMapping[,7])
            chrCircosMapping<-chrCircosMapping[with(chrCircosMapping, order(read_count,decreasing = FALSE)), ]
            ccm_dim=dim(chrCircosMapping)[1]
            for(ccm in 1:10)
            {
              ccm_start=1+ceiling((ccm_dim/10)*(ccm-1))
              ccm_end=ceiling((ccm_dim/10)*ccm)
              circos(R=340, cir=chrCircosDb, W=40,  mapping=chrCircosMapping[ccm_start:ccm_end,], type="link2", lwd=circosLineWidth*ccm,col=colors[ccm]);
            }
          }else
          {
            circos(R=340, cir=chrCircosDb, W=40,  mapping=chrCircosMapping, type="link2", lwd=circosLineWidth,col=circosLinecolor);
          }

          all_wig_hm=NULL
          all_bed_wig_hm=NULL
          print(paste("plot ",chrName,"circos picture",sep=""))

          if(all_wig_num>0)
          {
            for(ww in 1:all_wig_num)
            {
              m_wig=load_wig(all_wig_file[ww],resolution*1000,chrName,chrTotSize,chrstart,chrend)


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

                if(is.null(all_bed_wig_hm))
                {
                  all_bed_wig_hm=m_wig_bed
                }else
                {
                  all_bed_wig_hm=cbind(all_bed_wig_hm,m_wig_bed[,3])
                }
                if(is.null(all_wig_hm))
                {
                  all_wig_hm=m_wig
                  all_wig_hm[,3]=(all_wig_hm[,3]-range(all_wig_hm[,3])[1])/(range(all_wig_hm[,3])[2]-range(all_wig_hm[,3])[1])
                  tmpname=colnames(all_wig_hm)
                  colnames(all_wig_hm)=c(tmpname[1:2],all_wig_name[ww])
                }else
                {
                  tmpname=colnames(all_wig_hm)
                  all_wig_hm=cbind(all_wig_hm,((m_wig[,3]-range(m_wig[,3])[1])/(range(m_wig[,3])[2]-range(m_wig[,3])[1])))
                  colnames(all_wig_hm)=c(tmpname,all_wig_name[ww])

                }

              }
            }
          }


          circos(R=360, cir=chrCircosDb, W=circosTrackWidth, mapping=all_wig_hm,  col.v=3,  type="heatmap2",  col.bar=TRUE, lwd=0.1, col="blue")
          circos(R=(360+circosTrackWidth), cir=chrCircosDb, W=1,   type="chr", print.chr.lab=FALSE, scale=TRUE)
          text(730,820,"track name Outside-to-inside",cex=0.55,family="mono")
          for(namei in 1:all_wig_num)
          {
            text(740,820-namei*20,all_wig_name[namei],cex=0.55,family="mono")
          }
          dev.off() # text(400,820,family="mono",wigFile3,cex=0.7)
          print(paste(chrName," circos plot finish",sep = ""))


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

statistical_analysis<-function(bedFile,wig_dir="wig",bedWindow=0,matrix_dir="hg19",outputpdf=TRUE,chrom="all",chrstart=0,chrend=0,resolution=100,groupNum=100,dist_method="euclidean",clust_method="complete",clust_label=TRUE,clust_k=5,threshold=0,hm_trace=TRUE)
{
  #source("calbed.R")
  matrix_name_dir=list.files(path=matrix_dir,full.names=F,pattern=".matrix")
  matrix_full_dir=list.files(path=matrix_dir,full.names=T,pattern=".matrix")
  m_bed=load_bed(bedFile)
  chrNum=length(matrix_name_dir)
  all_wig_file=list.files(path=wig_dir,full.names=T,pattern=".wig")
  all_wig_name=list.files(path=wig_dir,full.names = FALSE,pattern = ".wig")
  all_wig_num=length(all_wig_file)

  if(chrom=="all")
  {
    for (i in 1:chrNum)
    {
      tmpNum=regexpr(".matrix",matrix_name_dir[i])
      chrName=substr(matrix_name_dir[i],1,tmpNum-1)
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

        logfile = paste(paste(matrix_dir,"/",chrName,"_statistic.txt",sep=""))
        if (file.exists(logfile) ==TRUE){file.remove(logfile)}
        starttime = paste("Analysis start time:" , as.character(Sys.time()))
        write(starttime,file=logfile,append=TRUE)
        print(starttime)



        if(all_wig_num>0)
        {
          clust_wig=NULL
          all_wig_info=NULL
          all_wig_info_loc=matrix(data = 0,nrow = all_wig_num,ncol = 2)
          print(paste(Sys.time()," start make clusters"))

          for(ww in 1:all_wig_num)
          {
            chrBedWig=load_bed_wig(all_wig_file[ww],chrBed,chrName,0,0,0)
            clust_wig=cbind(clust_wig,chrBedWig[,4])
            if(is.null(all_wig_info))
            {
              wig_info_start=1
            }else
            {
              wig_info_start=dim(all_wig_info)[1]+1
            }
            one_wig_info=load_all_wig(all_wig_file[ww])
            one_wig_info=one_wig_info[one_wig_info$chr==chrName,]
            one_wig_info[,2]=as.numeric(one_wig_info[,2])
            one_wig_info[,3]=as.numeric(one_wig_info[,3])
            one_wig_info[,4]=as.numeric(one_wig_info[,4])
            if(chrend>0)
            {
              one_wig_info=one_wig_info[which(one_wig_info[,3]<chrend),]
            }
            if(chrstart>0)
            {
              one_wig_info[,2]=one_wig_info[,2]-chrstart
              one_wig_info[,3]=one_wig_info[,3]-chrstart
              one_wig_info=one_wig_info[which(one_wig_info[,2]>=0),]
            }
            wig_info_end=wig_info_start+dim(one_wig_info)[1]-1
            all_wig_info_loc[ww,1]=wig_info_start
            all_wig_info_loc[ww,2]=wig_info_end
            all_wig_info=rbind(all_wig_info,one_wig_info)
          }
          if((is.null(clust_wig))==FALSE)
          {
            for(cii in 1:all_wig_num)
            {
              clust_wig[is.na(clust_wig[,cii]),cii]=0
              clust_wig[is.nan(clust_wig[,cii]),cii]=0
            }
          }



          if((is.null(clust_wig))==FALSE)
          {
            tt=c(1:(dim(clust_wig)[1]))
            clust_name=paste(tt,"_",chrBed[,1],":",(chrBed[,2]+chrstart),"-",(chrBed[,3]+chrstart),sep = "")
            row.names(clust_wig)=clust_name

            colnames(clust_wig)=all_wig_name
            mydata=clust_wig

            for(cii in 1:all_wig_num)
            {
              mydata[is.na(mydata[,cii]),cii]=0
              mydata[is.nan(mydata[,cii]),cii]=0


            }


            suppressPackageStartupMessages(library("lattice"))
            suppressPackageStartupMessages(library("flexclust"))
            if(all_wig_num>1)
            {
              bcl <- bootFlexclust(mydata, k=2:7, nboot=50, FUN=cclust, multicore=FALSE)
              if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
              {
                pdf(paste(matrix_dir,"/",chrName,"_cluster_k_density.pdf",sep=""),width = 8,height = 8)

              }else
              {
                jpeg(paste(matrix_dir,"/",chrName,"_cluster_k_density.jpeg",sep=""),width=1000,height=1000,quality = 100)
              }
              plot(bcl)
              densityplot(bcl, from=0)

              dev.off()
            }

            out.dist=dist(mydata,method=dist_method) #manhattan,euclidean,minkowski,chebyshev,mahalanobis,canberra
            out.hclust=hclust(out.dist,method=clust_method) #average,centroid,median,complete,single,ward.D,density
            if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
            {
              pdf(paste(matrix_dir,"/",chrName,"_cluster_tree.pdf",sep=""),width = 8,height = 8)

            }else
            {
              jpeg(paste(matrix_dir,"/",chrName,"_cluster_tree.jpeg",sep=""),width=1000,height=1000,quality = 100)
            }

            if((clust_label==TRUE)||(clust_label=="TRUE")||(clust_label=="true"))
            {
              plot(out.hclust)

            }else
            {
              ll_length=length(out.hclust$labels)
              tmp_label=out.hclust
              for(lli in 1:ll_length)
              {
                tmp_label$labels[lli]=""
              }

              plot(tmp_label)
              cc_list=rect.hclust(tmp_label,clust_k)

            }

            cc_list=rect.hclust(out.hclust,clust_k)
            cluster.id=cutree(out.hclust,clust_k)
            dev.off()
            row.names(chrBed)=clust_name
            m_ttest_result=NULL
            for(ww in 1:all_wig_num)
            {
              tmp_ttest_result=data.frame("wig1_pvalue"=numeric(dim(clust_wig)[1]),"wig1_difference"=character(dim(clust_wig)[1]),stringsAsFactors=FALSE)
              tmp_wig_ttest=NULL
              write("",file=logfile,append=TRUE)
              write("",file=logfile,append=TRUE)
              write("",file=logfile,append=TRUE)

              for(iii in 1:(dim(clust_wig)[1]))
              {
                if(is.na(clust_wig[iii,ww])==FALSE)
                {
                  tmp_wig_ttest=rbind(tmp_wig_ttest,t.test(all_wig_info[(all_wig_info_loc[ww,1]:all_wig_info_loc[ww,2]),4],mu=clust_wig[iii,ww]))
                }else
                {
                  tmp_wig_ttest=rbind(tmp_wig_ttest,t.test(all_wig_info[(all_wig_info_loc[ww,1]:all_wig_info_loc[ww,2]),4],mu=0))
                }
              }
              m_wig_equal=which(tmp_wig_ttest[,3]>0.05)
              m_wig_mean=mean(all_wig_info[(all_wig_info_loc[ww,1]:all_wig_info_loc[ww,2]),4])
              m_wig_nequal=which(tmp_wig_ttest[,3]<=0.05)
              m_wig_more=m_wig_nequal[which(clust_wig[m_wig_nequal[],ww]>m_wig_mean)]
              m_wig_less=m_wig_nequal[which(clust_wig[m_wig_nequal[],ww]<=m_wig_mean)]
              tmp_ttest_result[,1]=as.data.frame(as.matrix(tmp_wig_ttest[,3]))
              tmp_ttest_result[,1]=as.numeric(tmp_ttest_result[,1])
              tmp_ttest_result[m_wig_more[],2]="more"
              tmp_ttest_result[m_wig_less[],2]="less"
              tmp_ttest_result[m_wig_equal[],2]="equal"
              if(is.null(m_ttest_result))
              {
                m_ttest_result=tmp_ttest_result
              }else
              {
                m_ttest_result=cbind(m_ttest_result,tmp_ttest_result)

              }
              wig_test=rbind(cbind(na.omit(clust_wig[,ww]),1),cbind(na.omit(all_wig_info[(all_wig_info_loc[ww,1]:all_wig_info_loc[ww,2]),4]),2))
              bed_wig_mean=mean(na.omit(clust_wig[,ww]))
              all_wig_mean=mean(na.omit(all_wig_info[(all_wig_info_loc[ww,1]:all_wig_info_loc[ww,2]),4]))
              wig_test=as.data.frame(wig_test)
              colnames(wig_test)=c("wig_value","group")
              rownames(wig_test)=c(1:(dim(wig_test)[1]))
              wig_test$group=as.factor(wig_test$group)
              wig_kruskal=kruskal.test(wig_value~group, data=wig_test)
              wig_kruskalmc=kruskalmc(wig_value~group, data=wig_test, probs=0.05)
              wig_mult <- oneway_test(wig_value~group, data=wig_test,
                                      ytrafo = function(data) trafo(data, numeric_trafo = rank),
                                      xtrafo = function(data) trafo(data, factor_trafo = function(x)
                                        model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                                      teststat = "max", distribution = approximate(B = 90000))
              wig_pvalue=pvalue(wig_mult, method = "single-step")
              write(paste("the statistic test between BED WIG and  WIG : ",all_wig_name[ww],sep=""),file=logfile,append=TRUE)
              write("",file=logfile,append=TRUE)
              write("test name : Kruskal-Wallis rank sum test",file=logfile,append=TRUE)
              write(paste("Kruskal-Wallis chi-squared : ",wig_kruskal$statistic,sep = ""),file=logfile,append=TRUE)
              write(paste("Kruskal-Wallis df : ",wig_kruskal$parameter,sep = ""),file=logfile,append=TRUE)
              write(paste("Kruskal-Wallis p value : ",wig_kruskal$p.value,sep = ""),file=logfile,append=TRUE)
              write("",file=logfile,append=TRUE)
              write("test name : Multiple comparison test after Kruskal-Wallis",file=logfile,append=TRUE)
              write(paste("significance level : ",wig_kruskalmc$signif.level,sep = ""),file=logfile,append=TRUE)
              write(paste("observed difference  : ",wig_kruskalmc$dif.com$obs.dif,sep = ""),file=logfile,append=TRUE)
              write(paste("critical difference  : ",wig_kruskalmc$dif.com$critical.dif,sep = ""),file=logfile,append=TRUE)
              write(paste("exist difference  : ",wig_kruskalmc$dif.com$difference,sep = ""),file=logfile,append=TRUE)
              write(paste("BED WIG mean  : ",bed_wig_mean,sep = ""),file=logfile,append=TRUE)
              write(paste("ALL WIG mean  : ",all_wig_mean,sep = ""),file=logfile,append=TRUE)
              write("",file=logfile,append=TRUE)
              write("",file=logfile,append=TRUE)
              write("",file=logfile,append=TRUE)
            }

            tt_name=NULL
            for(ww in 1:all_wig_num)
            {
              tt_name=c(tt_name,paste(all_wig_name[ww],"_pvalue",sep = ""),paste(all_wig_name[ww],"_difference",sep=""))
            }
            colnames(m_ttest_result)=tt_name
            clust_group=cbind(na.omit(cbind(chrBed,clust_wig,m_ttest_result)),cluster.id)

            write.csv(clust_group,paste(matrix_dir,"/",chrName,"_cluster.csv",sep=""),row.names = FALSE)
            clust_heatmap=NULL


            clust_order_num=NULL
            tmp_order_num=0

            for(ii in 1:clust_k)
            {
              clust_order_num=rbind(clust_order_num,length(cc_list[[ii]]))
            }
            clust_order_num2=clust_order_num
            for(ii in 1:clust_k)
            {
              clust_order_num[ii]=sum(clust_order_num2[1:ii,])
            }
            clust_order_num[clust_k]=clust_order_num[clust_k]+1

            for( ii in 1:clust_k)
            {
              clust_heatmap=rbind(clust_heatmap,clust_group[which(clust_group[,"cluster.id"]==ii),])

            }

            rownames(clust_heatmap)=NULL
            print(paste(Sys.time()," print cluster heatmap"))

            if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
            {
              pdf(paste(matrix_dir,"/",chrName,"_cluster_heatmap.pdf",sep=""),width = 8,height = 8)

            }else
            {
              jpeg(paste(matrix_dir,"/",chrName,"_cluster_heatmap.jpeg",sep=""),width=1000,height=1000,quality = 100)
            }
            hm_data=as.matrix(clust_group[,5:(5+all_wig_num-1)])
            # for(hmi in 1:dim(hm_data)[2])
            #{
            #   hm_data[,hmi]=((hm_data[,hmi]-range(hm_data[,hmi])[1])/(range(hm_data[,hmi])[2]-range(hm_data[,hmi])[1])-0.5)*2
            #
            # }
            for(hmi in 1:dim(hm_data)[2])
            {
              hm_data[,hmi]=(hm_data[,hmi]-range(hm_data[,hmi])[1])/(range(hm_data[,hmi])[2]-range(hm_data[,hmi])[1])

            }
            clust_col_num=c(1:dim(hm_data)[2])

            if(hm_trace==TRUE)
            {
              heatmap.2(hm_data[out.hclust$order,],rowsep = clust_order_num,sepcolor="black",sepwidth = c(0.1,0.1),srtCol = 25,adjCol = c(0.6,0.8),cexCol = 0.7,col=whitered,dendrogram = "none",Rowv=FALSE,Colv=FALSE,adjRow = c(-500,-500),
                        breaks=256,
                        key.title=NA,
                        key.xlab=NA,
                        key.par=list(mgp=c(1.5, 0.5, 0),
                                     mar=c(1, 2.5, 1, 0)),
                        key.xtickfun=function() {
                          cex <- par("cex")*par("cex.axis")
                          side <- 1
                          line <- 0
                          col <- par("col.axis")
                          font <- par("font.axis")
                          mtext("low", side=side, at=0, adj=0,
                                line=line, cex=cex, col=col, font=font)
                          mtext("high", side=side, at=1, adj=1,
                                line=line, cex=cex, col=col, font=font)
                          return(list(labels=FALSE, tick=FALSE))
                        })
            }else
            {
              heatmap.2(hm_data[out.hclust$order,],rowsep = clust_order_num,sepcolor="black",sepwidth = c(0.1,0.1),srtCol = 25,adjCol = c(0.6,0.8),cexCol = 0.7,col=whitered,dendrogram = "none",Rowv=FALSE,Colv=FALSE,adjRow = c(-500,-500),
                        breaks=256,
                        trace = "none",
                        key.title=NA,
                        key.xlab=NA,
                        key.par=list(mgp=c(1.5, 0.5, 0),
                                     mar=c(1, 2.5, 1, 0)),
                        key.xtickfun=function() {
                          cex <- par("cex")*par("cex.axis")
                          side <- 1
                          line <- 0
                          col <- par("col.axis")
                          font <- par("font.axis")
                          mtext("low", side=side, at=0, adj=0,
                                line=line, cex=cex, col=col, font=font)
                          mtext("high", side=side, at=1, adj=1,
                                line=line, cex=cex, col=col, font=font)
                          return(list(labels=FALSE, tick=FALSE))
                        })
            }

            #heatmap(as.matrix(clust_heatmap[,5:(5+all_wig_num-1)]),Rowv=NA,Colv=NA,cexCol = 1,labCol = "")
            dev.off()

          }
        }



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
        write(paste("b2b frequency mean  : ",mean(bb_info[,3]),sep = ""),file=logfile,append=TRUE)
        write(paste("b2o frequency mean  : ",mean(bedTOall_info[,3]),sep = ""),file=logfile,append=TRUE)
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
      tmpNum=regexpr(".matrix",matrix_name_dir[i])
      chrName=substr(matrix_name_dir[i],1,tmpNum-1)
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




          logfile = paste(paste(matrix_dir,"/",chrName,"_statistic.txt",sep=""))
          if (file.exists(logfile) ==TRUE){file.remove(logfile)}
          starttime = paste("Analysis start time:" , as.character(Sys.time()))
          write(starttime,file=logfile,append=TRUE)

          print(starttime)
          if(all_wig_num>0)
          {
            clust_wig=NULL
            all_wig_info=NULL
            all_wig_info_loc=matrix(data = 0,nrow = all_wig_num,ncol = 2)
            print(paste(Sys.time()," start make clusters"))

            for(ww in 1:all_wig_num)
            {
              chrBedWig=load_bed_wig(all_wig_file[ww],chrBed,chrName,chrstart,chrend,0)
              clust_wig=cbind(clust_wig,chrBedWig[,4])
              if(is.null(all_wig_info))
              {
                wig_info_start=1
              }else
              {
                wig_info_start=dim(all_wig_info)[1]+1
              }
              one_wig_info=load_all_wig(all_wig_file[ww])
              one_wig_info=one_wig_info[one_wig_info$chr==chrom,]
              one_wig_info[,2]=as.numeric(one_wig_info[,2])
              one_wig_info[,3]=as.numeric(one_wig_info[,3])
              one_wig_info[,4]=as.numeric(one_wig_info[,4])
              if(chrend>0)
              {
                one_wig_info=one_wig_info[which(one_wig_info[,3]<chrend),]
              }
              if(chrstart>0)
              {
                one_wig_info[,2]=one_wig_info[,2]-chrstart
                one_wig_info[,3]=one_wig_info[,3]-chrstart
                one_wig_info=one_wig_info[which(one_wig_info[,2]>=0),]
              }
              wig_info_end=wig_info_start+dim(one_wig_info)[1]-1
              all_wig_info_loc[ww,1]=wig_info_start
              all_wig_info_loc[ww,2]=wig_info_end
              all_wig_info=rbind(all_wig_info,one_wig_info)
            }

            if((is.null(clust_wig))==FALSE)
            {
              for(cii in 1:all_wig_num)
              {
                clust_wig[is.na(clust_wig[,cii]),cii]=0
                clust_wig[is.nan(clust_wig[,cii]),cii]=0
              }
            }



            if((is.null(clust_wig))==FALSE)
            {
              tt=c(1:(dim(clust_wig)[1]))
              clust_name=paste(tt,"_",chrBed[,1],":",(chrBed[,2]+chrstart),"-",(chrBed[,3]+chrstart),sep = "")
              row.names(clust_wig)=clust_name

              colnames(clust_wig)=all_wig_name
              mydata=clust_wig

              for(cii in 1:all_wig_num)
              {
                mydata[is.na(mydata[,cii]),cii]=0
                mydata[is.nan(mydata[,cii]),cii]=0


              }



              suppressPackageStartupMessages(library("lattice"))
              suppressPackageStartupMessages(library("flexclust"))
              if(all_wig_num>1)
              {
                bcl <- bootFlexclust(mydata, k=2:7, nboot=50, FUN=cclust, multicore=FALSE)
                if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
                {
                  pdf(paste(matrix_dir,"/",chrName,"_cluster_k_density.pdf",sep=""),width = 8,height = 8)

                }else
                {
                  jpeg(paste(matrix_dir,"/",chrName,"_cluster_k_density.jpeg",sep=""),width=1000,height=1000,quality = 100)
                }
                plot(bcl)
                densityplot(bcl, from=0)

                dev.off()
              }

              out.dist=dist(mydata,method=dist_method) #manhattan,euclidean,minkowski,chebyshev,mahalanobis,canberra
              out.hclust=hclust(out.dist,method=clust_method) #average,centroid,median,complete,single,ward.D,density
              if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
              {
                pdf(paste(matrix_dir,"/",chrName,"_cluster_tree.pdf",sep=""),width = 8,height = 8)

              }else
              {
                jpeg(paste(matrix_dir,"/",chrName,"_cluster_tree.jpeg",sep=""),width=1000,height=1000,quality = 100)
              }

              if((clust_label==TRUE)||(clust_label=="TRUE")||(clust_label=="true"))
              {
                plot(out.hclust)

              }else
              {
                ll_length=length(out.hclust$labels)
                tmp_label=out.hclust
                for(lli in 1:ll_length)
                {
                  tmp_label$labels[lli]=""
                }

                plot(tmp_label)
                cc_list=rect.hclust(tmp_label,clust_k)

              }

              cc_list=rect.hclust(out.hclust,clust_k)
              cluster.id=cutree(out.hclust,clust_k)
              dev.off()
              row.names(chrBed)=clust_name
              m_ttest_result=NULL
              for(ww in 1:all_wig_num)
              {
                tmp_ttest_result=data.frame("wig1_pvalue"=numeric(dim(clust_wig)[1]),"wig1_difference"=character(dim(clust_wig)[1]),stringsAsFactors=FALSE)
                tmp_wig_ttest=NULL
                write("",file=logfile,append=TRUE)
                write("",file=logfile,append=TRUE)
                write("",file=logfile,append=TRUE)

                for(iii in 1:(dim(clust_wig)[1]))
                {
                  if(is.na(clust_wig[iii,ww])==FALSE)
                  {
                    tmp_wig_ttest=rbind(tmp_wig_ttest,t.test(all_wig_info[(all_wig_info_loc[ww,1]:all_wig_info_loc[ww,2]),4],mu=clust_wig[iii,ww]))
                  }else
                  {
                    tmp_wig_ttest=rbind(tmp_wig_ttest,t.test(all_wig_info[(all_wig_info_loc[ww,1]:all_wig_info_loc[ww,2]),4],mu=0))
                  }
                }
                m_wig_equal=which(tmp_wig_ttest[,3]>0.05)
                m_wig_mean=mean(all_wig_info[(all_wig_info_loc[ww,1]:all_wig_info_loc[ww,2]),4])
                m_wig_nequal=which(tmp_wig_ttest[,3]<=0.05)
                m_wig_more=m_wig_nequal[which(clust_wig[m_wig_nequal[],ww]>m_wig_mean)]
                m_wig_less=m_wig_nequal[which(clust_wig[m_wig_nequal[],ww]<=m_wig_mean)]
                tmp_ttest_result[,1]=as.data.frame(as.matrix(tmp_wig_ttest[,3]))
                tmp_ttest_result[,1]=as.numeric(tmp_ttest_result[,1])
                tmp_ttest_result[m_wig_more[],2]="more"
                tmp_ttest_result[m_wig_less[],2]="less"
                tmp_ttest_result[m_wig_equal[],2]="equal"
                if(is.null(m_ttest_result))
                {
                  m_ttest_result=tmp_ttest_result
                }else
                {
                  m_ttest_result=cbind(m_ttest_result,tmp_ttest_result)

                }
                wig_test=rbind(cbind(na.omit(clust_wig[,ww]),1),cbind(na.omit(all_wig_info[(all_wig_info_loc[ww,1]:all_wig_info_loc[ww,2]),4]),2))
                bed_wig_mean=mean(na.omit(clust_wig[,ww]))
                all_wig_mean=mean(na.omit(all_wig_info[(all_wig_info_loc[ww,1]:all_wig_info_loc[ww,2]),4]))
                wig_test=as.data.frame(wig_test)
                colnames(wig_test)=c("wig_value","group")
                rownames(wig_test)=c(1:(dim(wig_test)[1]))
                wig_test$group=as.factor(wig_test$group)
                wig_kruskal=kruskal.test(wig_value~group, data=wig_test)
                wig_kruskalmc=kruskalmc(wig_value~group, data=wig_test, probs=0.05)
                wig_mult <- oneway_test(wig_value~group, data=wig_test,
                                        ytrafo = function(data) trafo(data, numeric_trafo = rank),
                                        xtrafo = function(data) trafo(data, factor_trafo = function(x)
                                          model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                                        teststat = "max", distribution = approximate(B = 90000))
                wig_pvalue=pvalue(wig_mult, method = "single-step")
                write(paste("the statistic test between BED WIG and  WIG : ",all_wig_name[ww],sep=""),file=logfile,append=TRUE)
                write("",file=logfile,append=TRUE)
                write("test name : Kruskal-Wallis rank sum test",file=logfile,append=TRUE)
                write(paste("Kruskal-Wallis chi-squared : ",wig_kruskal$statistic,sep = ""),file=logfile,append=TRUE)
                write(paste("Kruskal-Wallis df : ",wig_kruskal$parameter,sep = ""),file=logfile,append=TRUE)
                write(paste("Kruskal-Wallis p value : ",wig_kruskal$p.value,sep = ""),file=logfile,append=TRUE)
                write("",file=logfile,append=TRUE)
                write("test name : Multiple comparison test after Kruskal-Wallis",file=logfile,append=TRUE)
                write(paste("significance level : ",wig_kruskalmc$signif.level,sep = ""),file=logfile,append=TRUE)
                write(paste("observed difference  : ",wig_kruskalmc$dif.com$obs.dif,sep = ""),file=logfile,append=TRUE)
                write(paste("critical difference  : ",wig_kruskalmc$dif.com$critical.dif,sep = ""),file=logfile,append=TRUE)
                write(paste("exist difference  : ",wig_kruskalmc$dif.com$difference,sep = ""),file=logfile,append=TRUE)
                write(paste("BED WIG mean  : ",bed_wig_mean,sep = ""),file=logfile,append=TRUE)
                write(paste("ALL WIG mean  : ",all_wig_mean,sep = ""),file=logfile,append=TRUE)
                write("",file=logfile,append=TRUE)
                write("",file=logfile,append=TRUE)
                write("",file=logfile,append=TRUE)
              }
              tt_name=NULL
              for(ww in 1:all_wig_num)
              {
                tt_name=c(tt_name,paste(all_wig_name[ww],"_pvalue",sep = ""),paste(all_wig_name[ww],"_difference",sep=""))
              }
              colnames(m_ttest_result)=tt_name
              clust_group=cbind(na.omit(cbind(chrBed,clust_wig,m_ttest_result)),cluster.id)

              write.csv(clust_group,paste(matrix_dir,"/",chrName,"_cluster.csv",sep=""),row.names = FALSE)
              clust_heatmap=NULL


              clust_order_num=NULL
              tmp_order_num=0

              for(ii in 1:clust_k)
              {
                clust_order_num=rbind(clust_order_num,length(cc_list[[ii]]))
              }
              clust_order_num2=clust_order_num
              for(ii in 1:clust_k)
              {
                clust_order_num[ii]=sum(clust_order_num2[1:ii,])
              }
              clust_order_num[clust_k]=clust_order_num[clust_k]+1


              for( ii in 1:clust_k)
              {
                clust_heatmap=rbind(clust_heatmap,clust_group[which(clust_group[,"cluster.id"]==ii),])

              }
              print(paste(Sys.time()," print cluster heatmap"))

              rownames(clust_heatmap)=NULL
              if((outputpdf==TRUE)||(outputpdf=="TRUE")||(outputpdf=="true"))
              {
                pdf(paste(matrix_dir,"/",chrName,"_cluster_heatmap.pdf",sep=""),width = 8,height = 8)

              }else
              {
                jpeg(paste(matrix_dir,"/",chrName,"_cluster_heatmap.jpeg",sep=""),width=1000,height=1000,quality = 100)
              }
              hm_data=as.matrix(clust_group[,5:(5+all_wig_num-1)])
              # for(hmi in 1:dim(hm_data)[2])
              #{
              #   hm_data[,hmi]=((hm_data[,hmi]-range(hm_data[,hmi])[1])/(range(hm_data[,hmi])[2]-range(hm_data[,hmi])[1])-0.5)*2
              #
              # }
              for(hmi in 1:dim(hm_data)[2])
              {
                hm_data[,hmi]=(hm_data[,hmi]-range(hm_data[,hmi])[1])/(range(hm_data[,hmi])[2]-range(hm_data[,hmi])[1])

              }
              clust_col_num=c(1:dim(hm_data)[2])


              if(hm_trace==TRUE)
              {
                heatmap.2(hm_data[out.hclust$order,],rowsep = clust_order_num,sepcolor="black",sepwidth = c(0.1,0.1),srtCol = 25,adjCol = c(0.6,0.8),cexCol = 0.7,col=whitered,dendrogram = "none",Rowv=FALSE,Colv=FALSE,adjRow = c(-500,-500),
                          breaks=256,
                          key.title=NA,
                          key.xlab=NA,
                          key.par=list(mgp=c(1.5, 0.5, 0),
                                       mar=c(1, 2.5, 1, 0)),
                          key.xtickfun=function() {
                            cex <- par("cex")*par("cex.axis")
                            side <- 1
                            line <- 0
                            col <- par("col.axis")
                            font <- par("font.axis")
                            mtext("low", side=side, at=0, adj=0,
                                  line=line, cex=cex, col=col, font=font)
                            mtext("high", side=side, at=1, adj=1,
                                  line=line, cex=cex, col=col, font=font)
                            return(list(labels=FALSE, tick=FALSE))
                          })
              }else
              {
                heatmap.2(hm_data[out.hclust$order,],rowsep = clust_order_num,sepcolor="black",sepwidth = c(0.1,0.1),srtCol = 25,adjCol = c(0.6,0.8),cexCol = 0.7,col=whitered,dendrogram = "none",Rowv=FALSE,Colv=FALSE,adjRow = c(-500,-500),
                          breaks=256,
                          trace = "none",
                          key.title=NA,
                          key.xlab=NA,
                          key.par=list(mgp=c(1.5, 0.5, 0),
                                       mar=c(1, 2.5, 1, 0)),
                          key.xtickfun=function() {
                            cex <- par("cex")*par("cex.axis")
                            side <- 1
                            line <- 0
                            col <- par("col.axis")
                            font <- par("font.axis")
                            mtext("low", side=side, at=0, adj=0,
                                  line=line, cex=cex, col=col, font=font)
                            mtext("high", side=side, at=1, adj=1,
                                  line=line, cex=cex, col=col, font=font)
                            return(list(labels=FALSE, tick=FALSE))
                          })
              }
              #heatmap(as.matrix(clust_heatmap[,5:(5+all_wig_num-1)]),Rowv=NA,Colv=NA,cexCol = 1,labCol = "")
              dev.off()

            }
          }

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
          write(paste("b2b frequency mean  : ",mean(bb_info[,3]),sep = ""),file=logfile,append=TRUE)
          write(paste("b2o frequency mean  : ",mean(bedTOall_info[,3]),sep = ""),file=logfile,append=TRUE)
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
