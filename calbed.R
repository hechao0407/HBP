#load bed files
load_bed<-function(bedfile)  
{
  
  #Reading all file and selecting the useful columns
  uniques<-read.table(file=bedfile, fill=TRUE, stringsAsFactors=FALSE,skip=1)
  uniques<-uniques[,c(1,2,3,4)]

  
  #Ordering and naming
  colnames(uniques)<-c("chr","start","end","name")
  uniques<-na.omit(uniques)
  ordereduniques<-uniques[with(uniques, order(chr,as.numeric(start))), ]
  if (dim(ordereduniques)[1]==0){
    stop('No data loaded! \n\n') 
  }
  else{
    return(ordereduniques)
  }
}



choose_chr_bed<-function(m_bed,chrom_list="chr2L")
{
  
  sec=NULL
  for (j in 1:length(chrom_list)){
    sec=rbind(sec,m_bed[m_bed$chr==chrom_list[j],])
  }
  return(sec)
}


#merge bed files to bin
check_bed_bin<-function(m_bed,bin=2000)
{
  ooo=matrix(data=0, nrow = dim(m_bed)[1], ncol = 4)
  for (i in 1:dim(m_bed)[1])
  {
    
    st=(ceiling(m_bed[i,2]/bin))
    ed=(ceiling(m_bed[i,3]/bin))
    st1=(ceiling((m_bed[i,2]-500)/bin))
    ed1=(ceiling((m_bed[i,3]+500)/bin))
    
    if(st!=ed)
    {
      ooo[i,2]=st
      ooo[i,3]=ed
      if(st1!=st)
      {
        ooo[i,1]=st1
      }
      if(ed1!=ed)
      {
        ooo[i,4]=ed1
      }
    }
    else
    {
      ooo[i,2]=st
      
      if(st1!=st)
      {
        ooo[i,1]=st1
      }
      if(ed1!=ed)
      {
        ooo[i,4]=ed1
      }
    }
    
  }
  return(ooo)
}


convert_bed_to_matrix<-function(bedfile,bin=2000,chrom_list="chr2L",tot_size=0,bed_window=2000){
  
  
  #Get only some chromosome
  if(tot_size==0)
  {
  	print("please input tot_size")
  }
  else
  {
  	cmap=matrix(data=0, nrow = tot_size, ncol = 1)
  	sec=NULL
	  for (j in 1:length(chrom_list)){
	    sec=rbind(sec,bedfile[bedfile$chr==chrom_list[j],])
	  }
	  ordereduniques=sec
	  
	   #Data insertion
	  count=0
	  for (i in 1:dim(ordereduniques)[1]){
      st=(ceiling(ordereduniques[i,2]/bin))
      ed=(ceiling(ordereduniques[i,3]/bin))
      st1=(ceiling((ordereduniques[i,2]-bed_window)/bin))
      ed1=(ceiling((ordereduniques[i,3]+bed_window)/bin))
      
      if(st<1)
      {
      	st=1
      }
      if(st1<1)
      {
      	st1=1
      }
      if(ed<1)
      {
      	ed=1
      }
      if(ed1<1)
      {
      	ed1=1
      }
      
      if(st>tot_size)
      {
      	st=tot_size
      }
      if(st1>tot_size)
      {
      	st1=tot_size
      }
      if(ed>tot_size)
      {
      	ed=tot_size
      }
      if(ed1>tot_size)
      {
      	ed1=tot_size
      }
      
      if(st!=ed)
      {
        cmap[st,1]=cmap[st,1]+1
        count=count+1
        cmap[ed,1]=cmap[ed,1]+1
        count=count+1
        if(st1!=st)
        {
          cmap[st1,1]=cmap[st1,1]+1
          count=count+1
        }
        if(ed1!=ed)
        {
          cmap[ed1,1]=cmap[ed1,1]+1
          count=count+1
        }
      }
      else
      {
        cmap[ed,1]=cmap[ed,1]+1
        count=count+1
        if(st1!=st)
        {
          cmap[st1,1]=cmap[st1,1]+1
          count=count+1
        }
        if(ed1!=ed)
        {
          cmap[ed1,1]=cmap[ed1,1]+1
          count=count+1
        }
      }
    }
      
	  
  }
  
  cat(sprintf("Mapped Fragments: %s\n",count))
  return(cmap)
}


mirny_cutoff<-function(cmap,th=1)
{
  ooo=matrix(data=0, nrow = dim(cmap)[1], ncol = dim(cmap)[1])
  rownames(ooo)<-rownames(cmap)
  colnames(ooo)<-colnames(cmap)
  n=dim(cmap)[1]
  
  
  for(i in 1:n)
  {
    
    for(j in i:n)
    {
      
      if(cmap[i,j]<th)
      {
        ooo[i,j]=0
      }
      else
      {
        ooo[i,j]=cmap[i,j]
      }
    }
  }
  
  return (ooo)
}


dif_analysis<-function(cmap,bed_file,group_num=100)
{
	lll=bed_file
  n=dim(lll)[1]
  n_count=0
  for (i in 1:n)
  {
    if(lll[i]!=0)
    {
      n_count=n_count+1
    }
  }
  random_group=matrix(data=0, nrow = group_num, ncol = dim(cmap)[1])
  random_result=matrix(data=0, nrow = group_num , ncol = 3)
  
  for (i in 1:group_num)
  {
    tmp_site=sample(1:n,size=n_count)
    random_group[i,tmp_site]= 2
  }
  
  
  for(t in 1:group_num)
  {
    for (i in 1:dim(lll)[1])
    {
      for (j in i:dim(lll)[1])
      {
        if(cmap[i,j]!=0)
        {
          if(random_group[t,i]!=0)
          {
            if(random_group[t,j]!=0)
            {
              random_result[t,1]=random_result[t,1]+1
            }
            random_result[t,2]=random_result[t,2]+1
          }
          random_result[t,3]=random_result[t,3]+1
          
        }
      }
    #  cat(sprintf("%d<------>%d,%d  \n",group_num,t,i))
    }
   
  }
  
  bed_count=0
  all_count=0
  mm_count=0
  
  for (i in 1:dim(lll)[1])
  {
    for (j in i:dim(lll)[1])
    {
      if(cmap[i,j]!=0)
      {
        if(lll[i]!=0)
        {
          if(lll[j]!=0)
          {
            mm_count=mm_count+1
          }
          bed_count=bed_count+1
        }
        all_count=all_count+1
        
      }
    }
  }
  
  cat(sprintf("random_mean %f,random_median %f,random_var %f,random_max %f,random_min %f \n",mean(random_result[,1]),median(random_result[,1]),var(random_result[,1]),max(random_result[,1]),min(random_result[,1])))
  cat(sprintf("%d<------>%d<------>%d  \n",all_count,bed_count,mm_count))
  #return (random_result)
  return(c(mean(random_result[,1]),median(random_result[,1]),var(random_result[,1]),max(random_result[,1]),min(random_result[,1]),mean(random_result[,2]),median(random_result[,2]),var(random_result[,2]),max(random_result[,2]),min(random_result[,2]),mean(random_result[,3]),median(random_result[,3]),var(random_result[,3]),max(random_result[,3]),min(random_result[,3]),mm_count,bed_count,all_count))
}


#output all beds that take part in the interaction
find_bed_info<-function(cmap,bed_matrix,bed_bin,bedfile,chr="chr4",st_cmap)
{
  n=dim(cmap)[1]
  bed_count=0
  nm_count=0
  all_count=0
  cccc=0;
  print(dim(bed_matrix)[1])
  print(dim(bed_bin)[1])
  ooo<-data.frame("bed_id"=numeric(0),"chr_nm"=character(0),"n_start"=numeric(0),"n_end"=numeric(0),"n_location"=numeric(0),"bed_seq"=character(0),"read_count"=numeric(0),stringsAsFactors=FALSE)

  for (i in 1:dim(bed_matrix)[1])
  {
    for (j in i:dim(bed_matrix)[1])
    {
      if(cmap[i,j]!=0)
      {
        if(bed_matrix[i]!=0)
        {
          for(m in 1:dim(bed_bin)[1])
          {
            for(k in 1:4)
            {
              if(bed_bin[m,k]!=0)
              {
                if(bed_bin[m,k]==i)
                {
                  if(bedfile[m,1]==chr)
                  {
                    cccc=cccc+1
                    ooo[cccc,1]=m
                    ooo[cccc,2]=bedfile[m,1]
                    ooo[cccc,3]=bedfile[m,2]
                    ooo[cccc,4]=bedfile[m,3]
                    ooo[cccc,5]=k
                    ooo[cccc,6]=bedfile[m,4]
                    ooo[cccc,7]=st_cmap[i,j]
                    #cccc=cccc+1 
                    #print(paste(cccc,"$$$$",cmap[i,j],"$$$$",bed_matrix[i],"$$$$",bed_bin[m,k],"$$$$",i,"$$$$",chr,sep=""))
                  }
                }
                
                
              }
            }
          }
        }        
      }
    }
  }
  
  
  return (ooo)
}


#output the bed that interact with our bedfile
find_bed_to_bed<-function(cmap,bed_matrix,bed_bin,bedfile,chr="chr4",st_cmap)
{
  n=dim(cmap)[1]
  bed_count=0
  nm_count=0
  all_count=0
  cccc=0;
  ooo<-data.frame("bed_id"=numeric(0),"chr_nm"=character(0),"n_start"=numeric(0),"n_end"=numeric(0),"n_location"=numeric(0),"bed_seq"=character(0),"read_count"=numeric(0),stringsAsFactors=FALSE)

  for (i in 1:dim(bed_matrix)[1])
  {
    for (j in i:dim(bed_matrix)[1])
    {
      if(cmap[i,j]!=0)
      {
        if(bed_matrix[i]!=0)
        {
          if(bed_matrix[j]!=0)
          {
            for(m in 1:dim(bed_bin)[1])
            {
              for(k in 1:4)
              {
                if(bed_bin[m,k]!=0)
                {
                  if(bed_bin[m,k]==i)
                  {
                    if(bedfile[m,1]==chr)
                    {
                      cccc=cccc+1
                      ooo[cccc,1]=m
                      ooo[cccc,2]=bedfile[m,1]
                      ooo[cccc,3]=bedfile[m,2]
                      ooo[cccc,4]=bedfile[m,3]
                      ooo[cccc,5]=k
                      ooo[cccc,6]=bedfile[m,4]
                      ooo[cccc,7]=st_cmap[i,j]
                      #print(cccc)
                     
                    }
                  }
                  
                  
                }
              }
            }
          }
        }        
      }
    }
  }
 
  
  return (ooo)
}



#output the interaction that our bed interact with ourbed
find_bed_to_bed_interaction2<-function(cmap,bed_matrix,bed_bin,bedfile,chr="chr4",st_cmap,m_threshold=0)
{
  n=dim(cmap)[1]
  bed_count=0
  nm_count=0
  all_count=0
  cccc=0;

  ooo<-data.frame("bed_id1"=numeric(0),"chr_nm1"=character(0),"n_start1"=numeric(0),"n_end1"=numeric(0),"n_location1"=numeric(0),"bed_seq1"=character(0),"bed_id2"=numeric(0),"chr_nm2"=character(0),"n_start2"=numeric(0),"n_end2"=numeric(0),"n_location2"=numeric(0),"bed_seq2"=character(0),"read_count"=numeric(0),"bed_bin1"=numeric(0),"bed_bin2"=numeric(0),stringsAsFactors=FALSE)
 
  for (i in 1:(dim(bed_matrix)[1]-1))
  {
    for (j in (i+1):dim(bed_matrix)[1])
    {
      if(cmap[i,j]>m_threshold)
      {
        if(bed_matrix[i]!=0)
        {
          if(bed_matrix[j]!=0)
          {
            for(m in 1:dim(bed_bin)[1])
            {
              for(k in 1:4)
              {
                if(bed_bin[m,k]!=0)
                {
                  if(bed_bin[m,k]==i)
                  {
                    if(bedfile[m,1]==chr)
                    {
                      for(mm in 1:dim(bed_bin)[1])
                      {
                        for(kk in 1:4)
                        {
                          if(bed_bin[mm,kk]!=0)
                          {
                            if(bed_bin[mm,kk]==j)
                            {
                              if(bedfile[mm,1]==chr)
                              {
                              	cccc=cccc+1
                                ooo[cccc,1]=m
                                ooo[cccc,2]=bedfile[m,1]
                                ooo[cccc,3]=bedfile[m,2]
                                ooo[cccc,4]=bedfile[m,3]
                                ooo[cccc,5]=k
                                ooo[cccc,6]=bedfile[m,4]
                                ooo[cccc,7]=mm
                                ooo[cccc,8]=bedfile[mm,1]
                                ooo[cccc,9]=bedfile[mm,2]
                                ooo[cccc,10]=bedfile[mm,3]
                                ooo[cccc,11]=kk
                                ooo[cccc,12]=bedfile[mm,4]
                                ooo[cccc,13]=st_cmap[i,j]
                                ooo[cccc,14]=i
                                ooo[cccc,15]=j
                                #cccc=cccc+1 
                                
                              }
                            }
                          }
                        }
                      }
                      
                    }
                  }
                }
              }
            }
          }
        }        
      }
    }
  }
  
  
  return (ooo)
}


#output bed that interact with node without our bed
find_bed_to_other<-function(cmap,bed_matrix,bed_bin,bedfile,chr="chr4",st_cmap)
{
  n=dim(cmap)[1]
  bed_count=0
  nm_count=0
  all_count=0
  cccc=0;

  ooo<-data.frame("bed_id"=numeric(0),"chr_nm"=character(0),"n_start"=numeric(0),"n_end"=numeric(0),"n_location"=numeric(0),"bed_seq"=character(0),"reads_count"=numeric(0),stringsAsFactors=FALSE)
  
  
  for (i in 1:dim(bed_matrix)[1])
  {
    for (j in i:dim(bed_matrix)[1])
    {
      if(cmap[i,j]!=0)
      {
        if(bed_matrix[i]!=0)
        {
          if(bed_matrix[j]==0)
          {
            for(m in 1:dim(bed_bin)[1])
            {
              for(k in 1:4)
              {
                if(bed_bin[m,k]!=0)
                {
                  if(bed_bin[m,k]==i)
                  {
                    if(bedfile[m,1]==chr)
                    {
                    	cccc=cccc+1
                      ooo[cccc,1]=m
                      ooo[cccc,2]=bedfile[m,1]
                      ooo[cccc,3]=bedfile[m,2]
                      ooo[cccc,4]=bedfile[m,3]
                      ooo[cccc,5]=k
                      ooo[cccc,6]=bedfile[m,4]
                      ooo[cccc,7]=st_cmap[i,j]
                      #cccc=cccc+1 
                    }
                  }
                }
              }
            }
          }
        }        
      }
    }
  }
  
  return (ooo)
}


#output the interaction that our bed interact with ourbed
calculate_omiccircos_data<-function(cmap,bed_matrix,bed_bin,bedfile,chr="chr4",st_cmap,m_th=0)
{
  n=dim(cmap)[1]
  bed_count=0
  nm_count=0
  all_count=0
  cccc=0;
  
  ooo<-data.frame("chr_nm1"=character(0),"n_start1"=numeric(0),"bed_seq1"=character(0),"chr_nm2"=character(0),"n_start2"=numeric(0),"bed_seq2"=character(0),"read_count"=numeric(0),stringsAsFactors=FALSE)
  
  for (i in 1:(dim(bed_matrix)[1]-1))
  {
    for (j in (i+1):dim(bed_matrix)[1])
    {
      if(cmap[i,j]>m_th)
      {
        if(bed_matrix[i]!=0)
        {
          if(bed_matrix[j]!=0)
          {
            for(m in 1:dim(bed_bin)[1])
            {
              for(k in 1:4)
              {
                if(bed_bin[m,k]!=0)
                {
                  if(bed_bin[m,k]==i)
                  {
                    if(bedfile[m,1]==chr)
                    {
                      for(mm in 1:dim(bed_bin)[1])
                      {
                        for(kk in 1:4)
                        {
                          if(bed_bin[mm,kk]!=0)
                          {
                            if(bed_bin[mm,kk]==j)
                            {
                              if(bedfile[mm,1]==chr)
                              {
                                cccc=cccc+1
                                ooo[cccc,1]=bedfile[m,1]
                                ooo[cccc,2]=bedfile[m,2]
                                ooo[cccc,3]=bedfile[m,4]
                                ooo[cccc,4]=bedfile[mm,1]
                                ooo[cccc,5]=bedfile[mm,2]
                                ooo[cccc,6]=bedfile[mm,4]
                                ooo[cccc,7]=st_cmap[i,j]
                                #cccc=cccc+1 
                                
                              }
                            }
                          }
                        }
                      }
                      
                    }
                  }
                }
              }
            }
          }
        }        
      }
    }
  }
  
  
  return (ooo)
}


load_wig<-function(wigfile,resolution=2000,chrom="chr2L",chrTotSize,chrstart,chrend)
{
  
  
  tmp_wig<-read.table(file=wigfile, fill=TRUE, stringsAsFactors=FALSE,skip=1)
  
  ooo<-data.frame("seg.name"=character(chrTotSize),"seg.po"=numeric(chrTotSize),"value"=numeric(chrTotSize),stringsAsFactors=FALSE)
  ooo_num<-data.frame("wig_num"=numeric(chrTotSize))
  sec=NULL
  
  sec=rbind(sec,tmp_wig[tmp_wig[,1]==chrom,])
  sec[,2]=as.numeric(sec[,2])
  sec[,3]=as.numeric(sec[,3])
  sec[,4]=as.numeric(sec[,4])
  if(chrend>0)
  {
    sec=sec[which(sec[,3]<chrend),]
    sec[,2]=sec[,2]-chrstart
    sec[,3]=sec[,3]-chrstart
    sec=sec[which(sec[,2]>0),]
  }
  
  
  
  secdim=dim(sec)[1]
  if(is.null(secdim))
  {
    return(NULL)
  }
  else if(secdim<2)
  {
    return(NULL)
  }
  else
  {
    ooo[,1]=chrom
    ooo[,2:3]=0
    ooo_num[,1]=0
    wig_bin=1
    sec_num=dim(sec)[1]
    
    
    for(i in 1:sec_num)
    {
      wig_bin=ceiling((sec[i,2]+sec[i,3])/(2*resolution))
      if(wig_bin>chrTotSize)
      {
        wig_bin=chrTotSize
      }
      ooo[wig_bin,3]=ooo[wig_bin,3]+sec[i,4]
      ooo_num[wig_bin,1]=ooo_num[wig_bin,1]+1
      ooo[wig_bin,2]=wig_bin
    }
    for(i in 1:chrTotSize)
    {
      if(ooo[i,2]==0)
      {
        ooo[i,2]=i
      }
      if(ooo_num[i,1]>0)
      {
        ooo[i,3]=(ooo[i,3])/(ooo_num[i,1])
        
      }
      else
      {
        ooo[i,3]=0
      }
    }
    return (ooo)
  }
  #tmp_wig=sec
  
}


convertToColors <- function(mat) {
  # Produce 'normalized' version of matrix, with values ranging from 0 to 1
  rng <- range(mat, na.rm = TRUE)
  m <- (mat - rng[1])/diff(rng)
  # Convert to a matrix of sRGB color strings
  m2 <- m; class(m2) <- "character"
  m2[!is.na(m2)] <- rgb(colorRamp(heat.colors(10))(m[!is.na(m)]), max = 255)
  m2[is.na(m2)] <- "transparent"
  return(m2)
}


convertToBedColors <- function(mat) {
  # Produce 'normalized' version of matrix, with values ranging from 0 to 1
  matlength=dim(mat)[1]
  pointwindow=5
  m0=mat
  for(i in 1:matlength)
  {
    a=which(mat[i,]==1)
    aLength=length(a)
    if(aLength>0)
    {
      if((i>pointwindow)&&(i<(matlength-pointwindow)))
      {
        for(j in 1:aLength)
        {
          if((a[j]>pointwindow)&&(a[j]<(matlength-pointwindow)))
          {
            m0[(i-pointwindow):(i+pointwindow),(a[j]-pointwindow):(a[j]+pointwindow)]=1
            
          }
          else if((a[j]<=pointwindow))
          {
            m0[(i-pointwindow):(i+pointwindow),1:(a[j]+pointwindow)]=1
          }
          else if(a[j]>=(matlength-pointwindow))
          {
            m0[(i-pointwindow):(i+pointwindow),(a[j]-pointwindow):matlength]=1
          }
        }
      }
      else if((i<=pointwindow))
      {
        for(j in 1:aLength)
        {
          if((a[j]>pointwindow)&&(a[j]<(matlength-pointwindow)))
          {
            m0[1:(i+pointwindow),(a[j]-pointwindow):(a[j]+pointwindow)]=1
            
          }
          else if((a[j]<=pointwindow))
          {
            m0[1:(i+pointwindow),1:(a[j]+pointwindow)]=1
          }
          else if(a[j]>=(matlength-pointwindow))
          {
            m0[1:(i+pointwindow),(a[j]-pointwindow):matlength]=1
          }
        }
      }
      else if((i>=(matlength-pointwindow)))
      {
        for(j in 1:aLength)
        {
          if((a[j]>pointwindow)&&(a[j]<(matlength-pointwindow)))
          {
            m0[(i-pointwindow):matlength,(a[j]-pointwindow):(a[j]+pointwindow)]=1
            
          }
          else if((a[j]<=pointwindow))
          {
            m0[(i-pointwindow):matlength,1:(a[j]+pointwindow)]=1
          }
          else if(a[j]>=(matlength-pointwindow))
          {
            m0[(i-pointwindow):matlength,(a[j]-pointwindow):matlength]=1
          }
        }
      }
    }
  }
  m=1-m0
  # Convert to a matrix of sRGB color strings
  m2 <- m; class(m2) <- "character"
  m2[!is.na(m2)] <- rgb(colorRamp(topo.colors(10))(m[!is.na(m)]), max = 255)
  m2[lower.tri(m2)] <- NA
  m2[is.na(m2)] <- "transparent"
  return(m2)
}



load_all_wig<-function(wigfile)  
{
  
  #Reading all file and selecting the useful columns
  uniques<-read.table(file=wigfile, fill=TRUE, stringsAsFactors=FALSE,skip=1)
  uniques<-uniques[,c(1,2,3,4)]
  
  
  #Ordering and naming
  colnames(uniques)<-c("chr","start","end","value")
  #unique[which(is.na(uniques[,4])),4]=0
  #uniques<-na.omit(uniques)
  #ordereduniques<-uniques[with(uniques, order(chr,as.numeric(start))), ]
  if (dim(uniques)[1]==0){
    stop('No data loaded! \n\n') 
  }
  else{
    return(uniques)
  }
}


load_bed_wig<-function(wigfile,chrBed,chrom,chrstart,chrend,m_win)
{
  all_wig=load_all_wig(wigfile)
  all_wig=all_wig[all_wig$chr==chrom,]
  all_wig[,2]=as.numeric(all_wig[,2])
  all_wig[,3]=as.numeric(all_wig[,3])
  m_result=data.frame("chrom"=character(0),"start"=numeric(0),"end"=numeric(0),"wig_value"=numeric(0),stringsAsFactors=FALSE)
  dim_bed=dim(chrBed)[1]
  m_result[1:dim_bed,1:3]=chrBed[,1:3]
  if(chrend>0)
  {
    all_wig=all_wig[which(all_wig[,3]<chrend),]
  }
  if(chrstart>0)
  {
    all_wig[,2]=all_wig[,2]-chrstart
    all_wig[,3]=all_wig[,3]-chrstart
    all_wig=all_wig[which(all_wig[,2]>=0),]
  }
  
  for(i in 1:dim_bed)
  {
    tmp_bed_wig=which(all_wig[,2]>=chrBed[i,2]-m_win)
    tmp_bed_wig=tmp_bed_wig[which(all_wig[tmp_bed_wig,3]<=chrBed[i,3]+m_win)]
    m_result[i,4]=mean(as.numeric(all_wig[tmp_bed_wig,4]))
  }
  return(m_result)
}

dm3_bed_convert<-function(bedfile,outbed)
{
  uniques<-read.table(file=bedfile, fill=TRUE, stringsAsFactors=FALSE)
  uniques<-uniques[,c(1,2,3,4)]
  colnames(uniques)<-c("chr","start","end","name")
  uniques<-na.omit(uniques)
  n_bed<-uniques[with(uniques, order(chr,as.numeric(start))), ]
  n_bed[which(n_bed[,1]=="2L"),1]="chr1"
  n_bed[which(n_bed[,1]=="2LHet"),1]="chr2"
  n_bed[which(n_bed[,1]=="2R"),1]="chr3"
  n_bed[which(n_bed[,1]=="2RHet"),1]="chr4"
  n_bed[which(n_bed[,1]=="3L"),1]="chr5"
  n_bed[which(n_bed[,1]=="3LHet"),1]="chr6"
  n_bed[which(n_bed[,1]=="3R"),1]="chr7"
  n_bed[which(n_bed[,1]=="3RHet"),1]="chr8"
  n_bed[which(n_bed[,1]=="4"),1]="chr9"
  n_bed[which(n_bed[,1]=="U"),1]="chr10"
  n_bed[which(n_bed[,1]=="Uextra"),1]="chr11"
  n_bed[which(n_bed[,1]=="X"),1]="chr12"
  n_bed[which(n_bed[,1]=="XHet"),1]="chr13"
  n_bed[which(n_bed[,1]=="YHet"),1]="chr14"
  n_bed[which(n_bed[,1]=="M"),1]="chr15"
  write.table(n_bed, file = outbed,quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE)
}
dm3_bed_convert2<-function(bedfile,outbed)
{
  uniques<-read.table(file=bedfile, fill=TRUE, stringsAsFactors=FALSE)
  uniques<-uniques[,c(1,2,3,4)]
  colnames(uniques)<-c("chr","start","end","name")
  uniques<-na.omit(uniques)
  n_bed<-uniques[with(uniques, order(chr,as.numeric(start))), ]
  n_bed[which(n_bed[,1]=="chr2L"),1]="chr1"
  n_bed[which(n_bed[,1]=="chr2LHet"),1]="chr2"
  n_bed[which(n_bed[,1]=="chr2R"),1]="chr3"
  n_bed[which(n_bed[,1]=="chr2RHet"),1]="chr4"
  n_bed[which(n_bed[,1]=="chr3L"),1]="chr5"
  n_bed[which(n_bed[,1]=="chr3LHet"),1]="chr6"
  n_bed[which(n_bed[,1]=="chr3R"),1]="chr7"
  n_bed[which(n_bed[,1]=="chr3RHet"),1]="chr8"
  n_bed[which(n_bed[,1]=="chr4"),1]="chr9"
  n_bed[which(n_bed[,1]=="chrU"),1]="chr10"
  n_bed[which(n_bed[,1]=="chrUextra"),1]="chr11"
  n_bed[which(n_bed[,1]=="chrX"),1]="chr12"
  n_bed[which(n_bed[,1]=="chrXHet"),1]="chr13"
  n_bed[which(n_bed[,1]=="chrYHet"),1]="chr14"
  n_bed[which(n_bed[,1]=="chrM"),1]="chr15"
  write.table(n_bed, file = outbed,quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE)
}


find_bed_to_bed_interaction<-function(cmap,bed_matrix,bed_bin,bedfile,chr="chr4",st_cmap,m_threshold=0)
{
  n=dim(cmap)[1]
  bed_count=0
  nm_count=0
  all_count=0
  cccc=0;
  
  ooo<-data.frame("bed_id1"=numeric(0),"chr_nm1"=character(0),"n_start1"=numeric(0),"n_end1"=numeric(0),"n_location1"=numeric(0),"bed_seq1"=character(0),"bed_id2"=numeric(0),"chr_nm2"=character(0),"n_start2"=numeric(0),"n_end2"=numeric(0),"n_location2"=numeric(0),"bed_seq2"=character(0),"read_count"=numeric(0),"bed_bin1"=numeric(0),"bed_bin2"=numeric(0),stringsAsFactors=FALSE)
  
  for(i in 1:(dim(bed_matrix)[1]-1))
  {
    jj=which(cmap[i,(i+1):(dim(bed_matrix)[1])]>m_threshold)
    jjj=jj+i
    if(bed_matrix[i]!=0)
    {
      jjjj=which(bed_matrix[jjj]!=0)
      jjjjj=jjj[jjjj]
      j_length=length(jjjjj)
      if(j_length>=1)
      {
        for(k in 1:4)
        {
          m_all=which(bed_bin[,k]==i)
          m_all_length=length(m_all)
          if(m_all_length>0)
          {
            for(mn in 1:m_all_length)
            {
              m=m_all[mn]
              if(bedfile[m,1]==chr)
              {
                for(hh in 1:j_length)
                {
                  mj=jjjjj[hh]
                  for(kk in 1:4)
                  {
                    mm_all=which(bed_bin[,kk]==mj)
                    mm_all_length=length(mm_all)
                    if(mm_all_length>0)
                    {
                      for(mmn in 1:mm_all_length)
                      {
                        mm=mm_all[mmn]
                        if(bedfile[mm,1]==chr)
                        {
                          cccc=cccc+1
                          ooo[cccc,1]=m
                          ooo[cccc,2]=bedfile[m,1]
                          ooo[cccc,3]=bedfile[m,2]
                          ooo[cccc,4]=bedfile[m,3]
                          ooo[cccc,5]=k
                          ooo[cccc,6]=bedfile[m,4]
                          ooo[cccc,7]=mm
                          ooo[cccc,8]=bedfile[mm,1]
                          ooo[cccc,9]=bedfile[mm,2]
                          ooo[cccc,10]=bedfile[mm,3]
                          ooo[cccc,11]=kk
                          ooo[cccc,12]=bedfile[mm,4]
                          ooo[cccc,13]=st_cmap[i,mj]
                          ooo[cccc,14]=i
                          ooo[cccc,15]=mj
                          #cccc=cccc+1 
                          
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return (ooo)
}
