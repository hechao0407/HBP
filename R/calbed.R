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


#output the interaction that our bed interact with ourbed

calculate_omiccircos_data<-function(cmap,bed_matrix,bed_bin,bedfile,chr="chr4",st_cmap,m_threshold=0)
{
  n=dim(cmap)[1]
  bed_count=0
  nm_count=0
  all_count=0
  cccc=0;

  ooo<-data.frame("chr_nm1"=character(0),"n_start1"=numeric(0),"bed_seq1"=character(0),"chr_nm2"=character(0),"n_start2"=numeric(0),"bed_seq2"=character(0),"read_count"=numeric(0),stringsAsFactors=FALSE)

  ooo=lapply(rbind(1:(dim(bed_matrix)[1]-1)), calculate_omiccircos_data_solo,cmap=cmap,bed_matrix=bed_matrix,bed_bin=bed_bin,bedfile=bedfile,chr=chr,st_cmap=st_cmap,m_threshold=m_threshold)
  final_ooo<-data.frame("chr_nm1"=character(0),"n_start1"=numeric(0),"bed_seq1"=character(0),"chr_nm2"=character(0),"n_start2"=numeric(0),"bed_seq2"=character(0),"read_count"=numeric(0),stringsAsFactors=FALSE)
  for(tt in 1:(dim(bed_matrix)[1]-1))
  {
    if(!is.null(ooo[tt][[1]]))
    {
      final_ooo=rbind(final_ooo,as.data.frame(ooo[tt][[1]]))
    }
  }
  colnames(final_ooo)=c("chr_nm1","n_start1","bed_seq1","chr_nm2","n_start2","bed_seq2","read_count")
  if(is.factor(final_ooo[,2]))
  {
    final_ooo[,2]=as.numeric(levels(final_ooo[,2])[final_ooo[,2]])
    final_ooo[,5]=as.numeric(levels(final_ooo[,5])[final_ooo[,5]])
    final_ooo[,7]=as.numeric(levels(final_ooo[,7])[final_ooo[,7]])
  }else
  {
    final_ooo[,2]=as.numeric(final_ooo[,2])
    final_ooo[,5]=as.numeric(final_ooo[,5])
    final_ooo[,7]=as.numeric(final_ooo[,7])
  }

 # final_ooo[,2]=as.numeric(levels(final_ooo[,2])[final_ooo[,2]])
 # final_ooo[,5]=as.numeric(levels(final_ooo[,5])[final_ooo[,5]])
 # final_ooo[,7]=as.numeric(levels(final_ooo[,7])[final_ooo[,7]])
#
 # final_ooo[,1]=as.character(levels(final_ooo[,1])[final_ooo[,1]])
 # final_ooo[,3]=as.character(levels(final_ooo[,3])[final_ooo[,3]])
 # final_ooo[,4]=as.character(levels(final_ooo[,4])[final_ooo[,4]])
 # final_ooo[,6]=as.character(levels(final_ooo[,6])[final_ooo[,6]])

  return (final_ooo)
}
calculate_omiccircos_data_solo<-function(i,cmap,bed_matrix,bed_bin,bedfile,chr="chr4",st_cmap,m_threshold=0)
{
  jj=which(cmap[i,(i+1):(dim(bed_matrix)[1])]>m_threshold)
  jjj=jj+i
  ooo=NULL
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
                        ooo=rbind(ooo,cbind(bedfile[m,1],bedfile[m,2],bedfile[m,4],bedfile[mm,1],bedfile[mm,2],bedfile[mm,4],st_cmap[i,mj]))
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

  tmp_ys=i%%100
  if(tmp_ys==0)
  {
    print(paste("circos bin ",i," finish",sep = ""))
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


find_bed_to_bed_interaction<-function(cmap,bed_matrix,bed_bin,bedfile,chr="chr4",st_cmap,m_threshold=0)
{
  n=dim(cmap)[1]
  bed_count=0
  nm_count=0
  all_count=0
  cccc=0;

  ooo<-data.frame("bed_id1"=numeric(0),"chr_nm1"=character(0),"n_start1"=numeric(0),"n_end1"=numeric(0),"n_location1"=numeric(0),"bed_seq1"=character(0),"bed_id2"=numeric(0),"chr_nm2"=character(0),"n_start2"=numeric(0),"n_end2"=numeric(0),"n_location2"=numeric(0),"bed_seq2"=character(0),"read_count"=numeric(0),"bed_bin1"=numeric(0),"bed_bin2"=numeric(0),stringsAsFactors=FALSE)

  ooo=lapply(rbind(1:(dim(bed_matrix)[1]-1)), find_bed_to_bed_interaction_solo,cmap=cmap,bed_matrix=bed_matrix,bed_bin=bed_bin,bedfile=bedfile,chr=chr,st_cmap=st_cmap,m_threshold=m_threshold)
  final_ooo<-data.frame("bed_id1"=numeric(0),"chr_nm1"=character(0),"n_start1"=numeric(0),"n_end1"=numeric(0),"n_location1"=numeric(0),"bed_seq1"=character(0),"bed_id2"=numeric(0),"chr_nm2"=character(0),"n_start2"=numeric(0),"n_end2"=numeric(0),"n_location2"=numeric(0),"bed_seq2"=character(0),"read_count"=numeric(0),"bed_bin1"=numeric(0),"bed_bin2"=numeric(0),stringsAsFactors=FALSE)
  for(tt in 1:(dim(bed_matrix)[1]-1))
  {
    if(!is.null(ooo[tt][[1]]))
    {
      final_ooo=rbind(final_ooo,as.data.frame(ooo[tt][[1]]))
    }
  }
  colnames(final_ooo)=c("bed_id1","chr_nm1","n_start1","n_end1","n_location1","bed_seq1","bed_id2","chr_nm2","n_start2","n_end2","n_location2","bed_seq2","read_count","bed_bin1","bed_bin2")
  if(is.factor(final_ooo[,3]))
  {
    final_ooo[,1]=as.numeric(levels(final_ooo[,1])[final_ooo[,1]])
    final_ooo[,3]=as.numeric(levels(final_ooo[,3])[final_ooo[,3]])
    final_ooo[,4]=as.numeric(levels(final_ooo[,4])[final_ooo[,4]])
    final_ooo[,5]=as.numeric(levels(final_ooo[,5])[final_ooo[,5]])
    final_ooo[,7]=as.numeric(levels(final_ooo[,7])[final_ooo[,7]])
    final_ooo[,10]=as.numeric(levels(final_ooo[,10])[final_ooo[,10]])
    final_ooo[,9]=as.numeric(levels(final_ooo[,9])[final_ooo[,9]])
    final_ooo[,11]=as.numeric(levels(final_ooo[,11])[final_ooo[,11]])
    final_ooo[,13]=as.numeric(levels(final_ooo[,13])[final_ooo[,13]])
    final_ooo[,14]=as.numeric(levels(final_ooo[,14])[final_ooo[,14]])
    final_ooo[,15]=as.numeric(levels(final_ooo[,15])[final_ooo[,15]])
    final_ooo[,2]=as.character(levels(final_ooo[,2])[final_ooo[,2]])
    final_ooo[,5]=as.character(levels(final_ooo[,5])[final_ooo[,5]])
    final_ooo[,6]=as.character(levels(final_ooo[,6])[final_ooo[,6]])
    final_ooo[,8]=as.character(levels(final_ooo[,8])[final_ooo[,8]])
    final_ooo[,12]=as.character(levels(final_ooo[,12])[final_ooo[,12]])
  }else
  {
    final_ooo[,3]=as.numeric(final_ooo[,3])
    final_ooo[,4]=as.numeric(final_ooo[,4])
    final_ooo[,9]=as.numeric(final_ooo[,9])
    final_ooo[,10]=as.numeric(final_ooo[,10])

    final_ooo[,13]=as.numeric(final_ooo[,13])
    final_ooo[,14]=as.numeric(final_ooo[,14])

    final_ooo[,15]=as.numeric(final_ooo[,15])


  }




  return (final_ooo)
}
find_bed_to_bed_interaction_solo<-function(i,cmap,bed_matrix,bed_bin,bedfile,chr="chr4",st_cmap,m_threshold=0)
{
  jj=which(cmap[i,(i+1):(dim(bed_matrix)[1])]>m_threshold)
  jjj=jj+i
  ooo=NULL
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
                        ooo=rbind(ooo,cbind(m,bedfile[m,1],bedfile[m,2],bedfile[m,3],k,bedfile[m,4],mm,bedfile[mm,1],bedfile[mm,2],bedfile[mm,3],kk,bedfile[mm,4],st_cmap[i,mj],i,mj))

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


whitered<-function (n)
{
  colorpanel(n, "white", "red")
}
