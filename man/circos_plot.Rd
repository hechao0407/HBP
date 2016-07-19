\name{circos_plot}
\alias{circos_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot the circos picture with tracks
}
\description{
plot the circos picture with tracks
}
\usage{
circos_plot(bedFile, wig_dir = "wig", matrix_dir = "hg19", bedWindow = 0, outputpdf = "TRUE", chrom = "all", chrstart = 0, chrend = 0, resolution = 100, circosHmSize = 0.1, circosCmpSize = 5, circosLineWidth = 0.01, circosLinecolor = "ReadCounts")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bedFile}{
%%     ~~Describe \code{bedFile} here~~
}
  \item{wig_dir}{
%%     ~~Describe \code{wig_dir} here~~
}
  \item{matrix_dir}{
%%     ~~Describe \code{matrix_dir} here~~
}
  \item{bedWindow}{
%%     ~~Describe \code{bedWindow} here~~
}
  \item{outputpdf}{
%%     ~~Describe \code{outputpdf} here~~
}
  \item{chrom}{
%%     ~~Describe \code{chrom} here~~
}
  \item{chrstart}{
%%     ~~Describe \code{chrstart} here~~
}
  \item{chrend}{
%%     ~~Describe \code{chrend} here~~
}
  \item{resolution}{
%%     ~~Describe \code{resolution} here~~
}
  \item{circosHmSize}{
%%     ~~Describe \code{circosHmSize} here~~
}
  \item{circosCmpSize}{
%%     ~~Describe \code{circosCmpSize} here~~
}
  \item{circosLineWidth}{
%%     ~~Describe \code{circosLineWidth} here~~
}
  \item{circosLinecolor}{
%%     ~~Describe \code{circosLinecolor} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (bedFile, wig_dir = "wig", matrix_dir = "hg19", bedWindow = 0, 
    outputpdf = "TRUE", chrom = "all", chrstart = 0, chrend = 0, 
    resolution = 100, circosHmSize = 0.1, circosCmpSize = 5, 
    circosLineWidth = 0.01, circosLinecolor = "ReadCounts") 
{
    source("calbed.R")
    gc()
    matrix_name_dir = list.files(path = matrix_dir, full.names = F, 
        pattern = ".matrix")
    matrix_full_dir = list.files(path = matrix_dir, full.names = T, 
        pattern = ".matrix")
    all_wig_file = list.files(path = wig_dir, full.names = T, 
        pattern = ".wig")
    all_wig_name = list.files(path = wig_dir, full.names = FALSE, 
        pattern = ".wig")
    all_wig_num = length(all_wig_file)
    m_bed = load_bed(bedFile)
    chrNum = length(matrix_name_dir)
    options(stringsAsFactors = FALSE)
    genomeFrame = data.frame(seg.name = character(0), seg.start = numeric(0), 
        seg.end = numeric(0), the.v = character(0), NO = character(0), 
        stringsAsFactors = FALSE)
    genomeFrameNum = 1
    for (i in 1:chrNum) {
        chrCmap = read.table(file = matrix_full_dir[i], fill = TRUE, 
            stringsAsFactors = FALSE)
        chrTotSize = dim(chrCmap)[1]
        tmpNum = regexpr(".matrix", matrix_name_dir[i])
        chrName = substr(matrix_name_dir[i], 1, tmpNum - 1)
        if (chrName == chrom) {
            if (chrend > 0) {
                tmpCmapStart = abs(ceiling((chrstart/(resolution * 
                  1000))))
                tmpCmapEnd = abs(ceiling((chrend/(resolution * 
                  1000))))
                if (tmpCmapEnd < tmpCmapStart) {
                  print("please input correct start and end number")
                  break
                }
                if (tmpCmapEnd > chrTotSize) {
                  tmpCmapEnd = chrTotSize
                }
                if (tmpCmapStart > chrTotSize) {
                  tmpCmapStart = chrTotSize
                }
                chrCmap = chrCmap[tmpCmapStart:tmpCmapEnd, tmpCmapStart:tmpCmapEnd]
                chrTotSize = dim(chrCmap)[1]
            }
        }
        for (ii in 1:chrTotSize) {
            genomeFrame[genomeFrameNum, 1] = chrName
            genomeFrame[genomeFrameNum, 2] = resolution * 1000 * 
                (ii - 1)
            genomeFrame[genomeFrameNum, 3] = resolution * 1000 * 
                ii
            genomeFrame[genomeFrameNum, 4] = NA
            genomeFrame[genomeFrameNum, 5] = NA
            genomeFrameNum = genomeFrameNum + 1
        }
        rm(chrCmap)
        gc()
    }
    if (chrom == "all") {
        for (i in 1:chrNum) {
            tmpNum = regexpr(".matrix", matrix_name_dir[i])
            chrName = substr(matrix_name_dir[i], 1, tmpNum - 
                1)
            chrBed = choose_chr_bed(m_bed, chrName)
            chrBedNum = dim(chrBed)[1]
            if (chrBedNum > 0) {
                chrCmap = read.table(file = matrix_full_dir[i], 
                  fill = TRUE, stringsAsFactors = FALSE)
                chrTotSize = dim(chrCmap)[1]
                tmpNum2 = regexpr("r", chrName)
                chrNo = substr(chrName, tmpNum2 + 1, nchar(chrName))
                chrBedBin = check_bed_bin(chrBed, resolution * 
                  1000)
                chrBedMatrix = convert_bed_to_matrix(chrBed, 
                  resolution * 1000, chrName, chrTotSize, bedWindow)
                chrCircosMapping = calculate_omiccircos_data(chrCmap, 
                  chrBedMatrix, chrBedBin, chrBed, chrName, chrCmap)
                chrCircosMapping[, 1] = chrNo
                chrCircosMapping[, 4] = chrNo
                chrCircosDb = segAnglePo(genomeFrame, seg = chrName)
                seg.num <- length(unique(genomeFrame[, 1]))
                colors <- rainbow(seg.num, alpha = 0.5)
                if ((outputpdf == "TRUE") || (outputpdf == "true")) {
                  pdf(paste(matrix_dir, "/", chrName, "_circos.pdf", 
                    sep = ""), width = 8, height = 8)
                }
                else {
                  jpeg(paste(matrix_dir, "/", chrName, "_circos.jpeg", 
                    sep = ""), width = 1000, height = 1000, quality = 100)
                }
                options(stringsAsFactors = FALSE)
                par(mar = c(2, 2, 2, 2))
                plot(c(1, 800), c(1, 800), type = "n", axes = FALSE, 
                  xlab = "", ylab = "", main = "")
                if (circosLinecolor == "rainbow") {
                  colors <- rainbow(seg.num, alpha = 0.5)
                  circosLinecolor = colors
                  circos(R = 340, cir = chrCircosDb, W = 40, 
                    mapping = chrCircosMapping, type = "link2", 
                    lwd = circosLineWidth, col = colors)
                }
                else if (circosLinecolor == "ReadCounts") {
                  colors <- c("#f7fbff", "#f7fbff", "#deebf7", 
                    "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", 
                    "#2171b5", "#08519c", "#08306b")
                  chrCircosMapping[, 7] = as.numeric(chrCircosMapping[, 
                    7])
                  chrCircosMapping <- chrCircosMapping[with(chrCircosMapping, 
                    order(read_count, decreasing = FALSE)), ]
                  ccm_dim = dim(chrCircosMapping)[1]
                  for (ccm in 1:10) {
                    ccm_start = 1 + ceiling((ccm_dim/10) * (ccm - 
                      1))
                    ccm_end = ceiling((ccm_dim/10) * ccm)
                    circos(R = 340, cir = chrCircosDb, W = 40, 
                      mapping = chrCircosMapping[ccm_start:ccm_end, 
                        ], type = "link2", lwd = circosLineWidth * 
                        ccm, col = colors[ccm])
                  }
                }
                else {
                  circos(R = 340, cir = chrCircosDb, W = 40, 
                    mapping = chrCircosMapping, type = "link2", 
                    lwd = circosLineWidth, col = circosLinecolor)
                }
                all_wig_hm = NULL
                all_bed_wig_hm = NULL
                print(paste("plot ", chrName, "circos picture", 
                  sep = ""))
                if (all_wig_num > 0) {
                  for (ww in 1:all_wig_num) {
                    m_wig = load_wig(all_wig_file[ww], resolution * 
                      1000, chrName, chrTotSize, chrstart, chrend)
                    chrBedWig = m_wig[which(chrBedMatrix[, 1] == 
                      1), ]
                    chrBedWig[, 2] = chrBedWig[, 2] * resolution * 
                      1000
                    m_wig[, 2] = m_wig[, 2] * resolution * 1000
                    m_wig_ttest = NULL
                    if ((length(chrBedWig)) > 0) {
                      for (iii in 1:(dim(chrBedWig)[1])) {
                        m_wig_ttest = rbind(m_wig_ttest, t.test(m_wig[, 
                          3], mu = chrBedWig[iii, 3]))
                      }
                      m_wig_equal = which(m_wig_ttest[, 3] > 
                        0.05)
                      m_wig_mean = mean(m_wig[, 3])
                      m_wig_nequal = which(m_wig_ttest[, 3] <= 
                        0.05)
                      m_wig_more = m_wig_nequal[which(chrBedWig[m_wig_nequal[], 
                        3] > m_wig_mean)]
                      m_wig_less = m_wig_nequal[which(chrBedWig[m_wig_nequal[], 
                        3] <= m_wig_mean)]
                      m_wig_bed = chrBedWig
                      m_wig_bed[, 3] = 0
                      m_wig_bed[m_wig_more, 3] = 1
                      m_wig_bed[m_wig_less, 3] = -1
                      if (is.null(all_bed_wig_hm)) {
                        all_bed_wig_hm = m_wig_bed
                      }
                      else {
                        all_bed_wig_hm = cbind(all_bed_wig_hm, 
                          m_wig_bed[, 3])
                      }
                      if (is.null(all_wig_hm)) {
                        all_wig_hm = m_wig
                        all_wig_hm[, 3] = (all_wig_hm[, 3] - 
                          range(all_wig_hm[, 3])[1])/(range(all_wig_hm[, 
                          3])[2] - range(all_wig_hm[, 3])[1])
                        tmpname = colnames(all_wig_hm)
                        colnames(all_wig_hm) = c(tmpname[1:2], 
                          all_wig_name[ww])
                      }
                      else {
                        tmpname = colnames(all_wig_hm)
                        all_wig_hm = cbind(all_wig_hm, ((m_wig[, 
                          3] - range(m_wig[, 3])[1])/(range(m_wig[, 
                          3])[2] - range(m_wig[, 3])[1])))
                        colnames(all_wig_hm) = c(tmpname, all_wig_name[ww])
                      }
                    }
                  }
                }
                circos(R = 360, cir = chrCircosDb, W = 40, mapping = all_wig_hm, 
                  col.v = 3, type = "heatmap2", col.bar = TRUE, 
                  lwd = 0.1, col = "blue")
                circos(R = 400, cir = chrCircosDb, W = 1, type = "chr", 
                  print.chr.lab = FALSE, scale = TRUE)
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
    else {
        for (i in 1:chrNum) {
            tmpNum = regexpr(".matrix", matrix_name_dir[i])
            chrName = substr(matrix_name_dir[i], 1, tmpNum - 
                1)
            if (chrName == chrom) {
                chrBed = choose_chr_bed(m_bed, chrName)
                if (chrend > 0) {
                  chrBed = chrBed[which(chrBed[, 3] < chrend), 
                    ]
                  chrBed[, 2] = chrBed[, 2] - chrstart
                  chrBed[, 3] = chrBed[, 3] - chrstart
                  chrBed = chrBed[which(chrBed[, 2] > 0), ]
                }
                chrBedNum = dim(chrBed)[1]
                if (chrBedNum > 0) {
                  chrCmap = read.table(file = matrix_full_dir[i], 
                    fill = TRUE, stringsAsFactors = FALSE)
                  chrTotSize = dim(chrCmap)[1]
                  if (chrend > 0) {
                    tmpCmapStart = abs(ceiling((chrstart/(resolution * 
                      1000))))
                    tmpCmapEnd = abs(ceiling((chrend/(resolution * 
                      1000))))
                    if (tmpCmapEnd < tmpCmapStart) {
                      print("please input correct start and end number")
                      break
                    }
                    if (tmpCmapEnd > chrTotSize) {
                      tmpCmapEnd = chrTotSize
                    }
                    if (tmpCmapStart > chrTotSize) {
                      tmpCmapStart = chrTotSize
                    }
                    chrCmap = chrCmap[tmpCmapStart:tmpCmapEnd, 
                      tmpCmapStart:tmpCmapEnd]
                    chrTotSize = dim(chrCmap)[1]
                  }
                  print(paste("plot ", chrName, "circos picture", 
                    sep = ""))
                  tmpNum2 = regexpr("r", chrName)
                  chrNo = substr(chrName, tmpNum2 + 1, nchar(chrName))
                  chrBedBin = check_bed_bin(chrBed, resolution * 
                    1000)
                  chrBedMatrix = convert_bed_to_matrix(chrBed, 
                    resolution * 1000, chrName, chrTotSize, bedWindow)
                  chrCircosMapping = calculate_omiccircos_data(chrCmap, 
                    chrBedMatrix, chrBedBin, chrBed, chrName, 
                    chrCmap)
                  chrCircosMapping[, 1] = chrNo
                  chrCircosMapping[, 4] = chrNo
                  chrCircosDb = segAnglePo(genomeFrame, seg = chrName)
                  seg.num <- length(unique(genomeFrame[, 1]))
                  if ((outputpdf == "TRUE") || (outputpdf == 
                    "true")) {
                    pdf(paste(matrix_dir, "/", chrName, "_circos.pdf", 
                      sep = ""), width = 8, height = 8)
                  }
                  else {
                    jpeg(paste(matrix_dir, "/", chrName, "_circos.jpeg", 
                      sep = ""), width = 1000, height = 1000, 
                      quality = 100)
                  }
                  options(stringsAsFactors = FALSE)
                  par(mar = c(2, 2, 2, 2))
                  plot(c(1, 800), c(1, 800), type = "n", axes = FALSE, 
                    xlab = "", ylab = "", main = "")
                  if (circosLinecolor == "rainbow") {
                    colors <- rainbow(seg.num, alpha = 0.5)
                    circosLinecolor = colors
                    circos(R = 340, cir = chrCircosDb, W = 40, 
                      mapping = chrCircosMapping, type = "link2", 
                      lwd = circosLineWidth, col = colors)
                  }
                  else if (circosLinecolor == "ReadCounts") {
                    colors <- c("#f7fbff", "#f7fbff", "#deebf7", 
                      "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", 
                      "#2171b5", "#08519c", "#08306b")
                    chrCircosMapping[, 7] = as.numeric(chrCircosMapping[, 
                      7])
                    chrCircosMapping <- chrCircosMapping[with(chrCircosMapping, 
                      order(read_count, decreasing = FALSE)), 
                      ]
                    ccm_dim = dim(chrCircosMapping)[1]
                    for (ccm in 1:10) {
                      ccm_start = 1 + ceiling((ccm_dim/10) * 
                        (ccm - 1))
                      ccm_end = ceiling((ccm_dim/10) * ccm)
                      circos(R = 340, cir = chrCircosDb, W = 40, 
                        mapping = chrCircosMapping[ccm_start:ccm_end, 
                          ], type = "link2", lwd = circosLineWidth * 
                          ccm, col = colors[ccm])
                    }
                  }
                  else {
                    circos(R = 340, cir = chrCircosDb, W = 40, 
                      mapping = chrCircosMapping, type = "link2", 
                      lwd = circosLineWidth, col = circosLinecolor)
                  }
                  all_wig_hm = NULL
                  all_bed_wig_hm = NULL
                  print(paste("plot ", chrName, "circos picture", 
                    sep = ""))
                  if (all_wig_num > 0) {
                    for (ww in 1:all_wig_num) {
                      m_wig = load_wig(all_wig_file[ww], resolution * 
                        1000, chrName, chrTotSize, chrstart, 
                        chrend)
                      chrBedWig = m_wig[which(chrBedMatrix[, 
                        1] == 1), ]
                      chrBedWig[, 2] = chrBedWig[, 2] * resolution * 
                        1000
                      m_wig[, 2] = m_wig[, 2] * resolution * 
                        1000
                      m_wig_ttest = NULL
                      if ((length(chrBedWig)) > 0) {
                        for (iii in 1:(dim(chrBedWig)[1])) {
                          m_wig_ttest = rbind(m_wig_ttest, t.test(m_wig[, 
                            3], mu = chrBedWig[iii, 3]))
                        }
                        m_wig_equal = which(m_wig_ttest[, 3] > 
                          0.05)
                        m_wig_mean = mean(m_wig[, 3])
                        m_wig_nequal = which(m_wig_ttest[, 3] <= 
                          0.05)
                        m_wig_more = m_wig_nequal[which(chrBedWig[m_wig_nequal[], 
                          3] > m_wig_mean)]
                        m_wig_less = m_wig_nequal[which(chrBedWig[m_wig_nequal[], 
                          3] <= m_wig_mean)]
                        m_wig_bed = chrBedWig
                        m_wig_bed[, 3] = 0
                        m_wig_bed[m_wig_more, 3] = 1
                        m_wig_bed[m_wig_less, 3] = -1
                        if (is.null(all_bed_wig_hm)) {
                          all_bed_wig_hm = m_wig_bed
                        }
                        else {
                          all_bed_wig_hm = cbind(all_bed_wig_hm, 
                            m_wig_bed[, 3])
                        }
                        if (is.null(all_wig_hm)) {
                          all_wig_hm = m_wig
                          all_wig_hm[, 3] = (all_wig_hm[, 3] - 
                            range(all_wig_hm[, 3])[1])/(range(all_wig_hm[, 
                            3])[2] - range(all_wig_hm[, 3])[1])
                          tmpname = colnames(all_wig_hm)
                          colnames(all_wig_hm) = c(tmpname[1:2], 
                            all_wig_name[ww])
                        }
                        else {
                          tmpname = colnames(all_wig_hm)
                          all_wig_hm = cbind(all_wig_hm, ((m_wig[, 
                            3] - range(m_wig[, 3])[1])/(range(m_wig[, 
                            3])[2] - range(m_wig[, 3])[1])))
                          colnames(all_wig_hm) = c(tmpname, all_wig_name[ww])
                        }
                      }
                    }
                  }
                  circos(R = 360, cir = chrCircosDb, W = 40, 
                    mapping = all_wig_hm, col.v = 3, type = "heatmap2", 
                    col.bar = TRUE, lwd = 0.1, col = "blue")
                  circos(R = 400, cir = chrCircosDb, W = 1, type = "chr", 
                    print.chr.lab = FALSE, scale = TRUE)
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
