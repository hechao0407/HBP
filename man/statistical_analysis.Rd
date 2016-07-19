\name{statistical_analysis}
\alias{statistical_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
statistical analysis and make clusters
}
\description{
statistical analysis and make clusters
}
\usage{
statistical_analysis(bedFile, wig_dir = "wig", bedWindow = 0, matrix_dir = "hg19", outputpdf = "TRUE", chrom = "all", chrstart = 0, chrend = 0, resolution = 100, groupNum = 100, dist_method = "euclidean", clust_method = "complete", clust_k = 5, threshold = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bedFile}{
%%     ~~Describe \code{bedFile} here~~
}
  \item{wig_dir}{
%%     ~~Describe \code{wig_dir} here~~
}
  \item{bedWindow}{
%%     ~~Describe \code{bedWindow} here~~
}
  \item{matrix_dir}{
%%     ~~Describe \code{matrix_dir} here~~
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
  \item{groupNum}{
%%     ~~Describe \code{groupNum} here~~
}
  \item{dist_method}{
%%     ~~Describe \code{dist_method} here~~
}
  \item{clust_method}{
%%     ~~Describe \code{clust_method} here~~
}
  \item{clust_k}{
%%     ~~Describe \code{clust_k} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
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
function (bedFile, wig_dir = "wig", bedWindow = 0, matrix_dir = "hg19", 
    outputpdf = "TRUE", chrom = "all", chrstart = 0, chrend = 0, 
    resolution = 100, groupNum = 100, dist_method = "euclidean", 
    clust_method = "complete", clust_k = 5, threshold = 0) 
{
    source("calbed.R")
    matrix_name_dir = list.files(path = matrix_dir, full.names = F, 
        pattern = ".matrix")
    matrix_full_dir = list.files(path = matrix_dir, full.names = T, 
        pattern = ".matrix")
    m_bed = load_bed(bedFile)
    chrNum = length(matrix_name_dir)
    all_wig_file = list.files(path = wig_dir, full.names = T, 
        pattern = ".wig")
    all_wig_name = list.files(path = wig_dir, full.names = FALSE, 
        pattern = ".wig")
    all_wig_num = length(all_wig_file)
    if (chrom == "all") {
        for (i in 1:chrNum) {
            tmpNum = regexpr(".matrix", matrix_name_dir[i])
            chrName = substr(matrix_name_dir[i], 1, tmpNum - 
                1)
            chrBed = choose_chr_bed(m_bed, chrName)
            chrBedNum = dim(chrBed)[1]
            print(matrix_full_dir[i])
            if (chrBedNum > 0) {
                print(chrName)
                chrCmap = read.table(file = matrix_full_dir[i], 
                  fill = TRUE, stringsAsFactors = FALSE)
                chrTotSize = dim(chrCmap)[1]
                chrBedBin = check_bed_bin(chrBed, resolution * 
                  1000)
                chrBedMatrix = convert_bed_to_matrix(chrBed, 
                  resolution * 1000, chrName, chrTotSize, bedWindow)
                logfile = paste(paste(matrix_dir, "/", chrName, 
                  "_statistic.txt", sep = ""))
                if (file.exists(logfile) == TRUE) {
                  file.remove(logfile)
                }
                starttime = paste("Analysis start time:", as.character(Sys.time()))
                write(starttime, file = logfile, append = TRUE)
                if (all_wig_num > 0) {
                  clust_wig = NULL
                  all_wig_info = NULL
                  all_wig_info_loc = matrix(data = 0, nrow = all_wig_num, 
                    ncol = 2)
                  for (ww in 1:all_wig_num) {
                    chrBedWig = load_bed_wig(all_wig_file[ww], 
                      chrBed, chrName, 0, 0, 0)
                    clust_wig = cbind(clust_wig, chrBedWig[, 
                      4])
                    if (is.null(all_wig_info)) {
                      wig_info_start = 1
                    }
                    else {
                      wig_info_start = dim(all_wig_info)[1] + 
                        1
                    }
                    one_wig_info = load_all_wig(all_wig_file[ww])
                    one_wig_info = one_wig_info[one_wig_info$chr == 
                      chrName, ]
                    one_wig_info[, 2] = as.numeric(one_wig_info[, 
                      2])
                    one_wig_info[, 3] = as.numeric(one_wig_info[, 
                      3])
                    one_wig_info[, 4] = as.numeric(one_wig_info[, 
                      4])
                    if (chrend > 0) {
                      one_wig_info = one_wig_info[which(one_wig_info[, 
                        3] < chrend), ]
                    }
                    if (chrstart > 0) {
                      one_wig_info[, 2] = one_wig_info[, 2] - 
                        chrstart
                      one_wig_info[, 3] = one_wig_info[, 3] - 
                        chrstart
                      one_wig_info = one_wig_info[which(one_wig_info[, 
                        2] >= 0), ]
                    }
                    wig_info_end = wig_info_start + dim(one_wig_info)[1] - 
                      1
                    all_wig_info_loc[ww, 1] = wig_info_start
                    all_wig_info_loc[ww, 2] = wig_info_end
                    all_wig_info = rbind(all_wig_info, one_wig_info)
                  }
                  if ((is.null(clust_wig)) == FALSE) {
                    tt = c(1:(dim(clust_wig)[1]))
                    clust_name = paste(tt, "_", chrBed[, 1], 
                      ":", (chrBed[, 2] + chrstart), "-", (chrBed[, 
                        3] + chrstart), sep = "")
                    row.names(clust_wig) = clust_name
                    colnames(clust_wig) = all_wig_name
                    mydata = na.omit(clust_wig)
                    suppressPackageStartupMessages(library("lattice"))
                    suppressPackageStartupMessages(library("flexclust"))
                    if (all_wig_num > 1) {
                      bcl <- bootFlexclust(mydata, k = 2:7, nboot = 50, 
                        FUN = cclust, multicore = FALSE)
                      if ((outputpdf == "TRUE") || (outputpdf == 
                        "true")) {
                        pdf(paste(matrix_dir, "/", chrName, "_cluster_k_density.pdf", 
                          sep = ""), width = 8, height = 8)
                      }
                      else {
                        jpeg(paste(matrix_dir, "/", chrName, 
                          "_cluster_k_density.jpeg", sep = ""), 
                          width = 1000, height = 1000, quality = 100)
                      }
                      plot(bcl)
                      densityplot(bcl, from = 0)
                      dev.off()
                    }
                    out.dist = dist(mydata, method = dist_method)
                    out.hclust = hclust(out.dist, method = clust_method)
                    if ((outputpdf == "TRUE") || (outputpdf == 
                      "true")) {
                      pdf(paste(matrix_dir, "/", chrName, "_cluster_tree.pdf", 
                        sep = ""), width = 8, height = 8)
                    }
                    else {
                      jpeg(paste(matrix_dir, "/", chrName, "_cluster_tree.jpeg", 
                        sep = ""), width = 1000, height = 1000, 
                        quality = 100)
                    }
                    plot(out.hclust)
                    rect.hclust(out.hclust, clust_k)
                    cluster.id = cutree(out.hclust, clust_k)
                    dev.off()
                    row.names(chrBed) = clust_name
                    m_ttest_result = NULL
                    for (ww in 1:all_wig_num) {
                      tmp_ttest_result = data.frame(wig1_pvalue = numeric(dim(clust_wig)[1]), 
                        wig1_difference = character(dim(clust_wig)[1]), 
                        stringsAsFactors = FALSE)
                      tmp_wig_ttest = NULL
                      write("", file = logfile, append = TRUE)
                      write("", file = logfile, append = TRUE)
                      write("", file = logfile, append = TRUE)
                      for (iii in 1:(dim(clust_wig)[1])) {
                        if (is.na(clust_wig[iii, 1]) == FALSE) {
                          tmp_wig_ttest = rbind(tmp_wig_ttest, 
                            t.test(all_wig_info[(all_wig_info_loc[ww, 
                              1]:all_wig_info_loc[ww, 2]), 4], 
                              mu = clust_wig[iii, 1]))
                        }
                        else {
                          tmp_wig_ttest = rbind(tmp_wig_ttest, 
                            t.test(all_wig_info[(all_wig_info_loc[ww, 
                              1]:all_wig_info_loc[ww, 2]), 4], 
                              mu = 0))
                        }
                      }
                      m_wig_equal = which(tmp_wig_ttest[, 3] > 
                        0.05)
                      m_wig_mean = mean(all_wig_info[(all_wig_info_loc[ww, 
                        1]:all_wig_info_loc[ww, 2]), 4])
                      m_wig_nequal = which(tmp_wig_ttest[, 3] <= 
                        0.05)
                      m_wig_more = m_wig_nequal[which(clust_wig[m_wig_nequal[], 
                        1] > m_wig_mean)]
                      m_wig_less = m_wig_nequal[which(clust_wig[m_wig_nequal[], 
                        1] <= m_wig_mean)]
                      tmp_ttest_result[, 1] = as.data.frame(as.matrix(tmp_wig_ttest[, 
                        3]))
                      tmp_ttest_result[, 1] = as.numeric(tmp_ttest_result[, 
                        1])
                      tmp_ttest_result[m_wig_more[], 2] = "more"
                      tmp_ttest_result[m_wig_less[], 2] = "less"
                      tmp_ttest_result[m_wig_equal[], 2] = "equal"
                      if (is.null(m_ttest_result)) {
                        m_ttest_result = tmp_ttest_result
                      }
                      else {
                        m_ttest_result = cbind(m_ttest_result, 
                          tmp_ttest_result)
                      }
                      wig_test = rbind(cbind(na.omit(clust_wig[, 
                        1]), 1), cbind(na.omit(all_wig_info[(all_wig_info_loc[ww, 
                        1]:all_wig_info_loc[ww, 2]), 4]), 2))
                      wig_test = as.data.frame(wig_test)
                      colnames(wig_test) = c("wig_value", "group")
                      rownames(wig_test) = c(1:(dim(wig_test)[1]))
                      wig_test$group = as.factor(wig_test$group)
                      wig_kruskal = kruskal.test(wig_value ~ 
                        group, data = wig_test)
                      wig_kruskalmc = kruskalmc(wig_value ~ group, 
                        data = wig_test, probs = 0.05)
                      wig_mult <- oneway_test(wig_value ~ group, 
                        data = wig_test, ytrafo = function(data) trafo(data, 
                          numeric_trafo = rank), xtrafo = function(data) trafo(data, 
                          factor_trafo = function(x) model.matrix(~x - 
                            1) \%*\% t(contrMat(table(x), "Tukey"))), 
                        teststat = "max", distribution = approximate(B = 90000))
                      wig_pvalue = pvalue(wig_mult, method = "single-step")
                      write(paste("the statistic test between BED WIG and  WIG : ", 
                        all_wig_name[ww], sep = ""), file = logfile, 
                        append = TRUE)
                      write("", file = logfile, append = TRUE)
                      write("test name : Kruskal-Wallis rank sum test", 
                        file = logfile, append = TRUE)
                      write(paste("Kruskal-Wallis chi-squared : ", 
                        wig_kruskal$statistic, sep = ""), file = logfile, 
                        append = TRUE)
                      write(paste("Kruskal-Wallis df : ", wig_kruskal$parameter, 
                        sep = ""), file = logfile, append = TRUE)
                      write(paste("Kruskal-Wallis p value : ", 
                        wig_kruskal$p.value, sep = ""), file = logfile, 
                        append = TRUE)
                      write("", file = logfile, append = TRUE)
                      write("test name : Multiple comparison test after Kruskal-Wallis", 
                        file = logfile, append = TRUE)
                      write(paste("significance level : ", wig_kruskalmc$signif.level, 
                        sep = ""), file = logfile, append = TRUE)
                      write(paste("observed difference  : ", 
                        wig_kruskalmc$dif.com$obs.dif, sep = ""), 
                        file = logfile, append = TRUE)
                      write(paste("critical difference  : ", 
                        wig_kruskalmc$dif.com$critical.dif, sep = ""), 
                        file = logfile, append = TRUE)
                      write(paste("exist difference  : ", wig_kruskalmc$dif.com$difference, 
                        sep = ""), file = logfile, append = TRUE)
                      write("", file = logfile, append = TRUE)
                      write("", file = logfile, append = TRUE)
                      write("", file = logfile, append = TRUE)
                    }
                    tt_name = NULL
                    for (ww in 1:all_wig_num) {
                      tt_name = c(tt_name, paste(all_wig_name[ww], 
                        "_pvalue", sep = ""), paste(all_wig_name[ww], 
                        "_difference", sep = ""))
                    }
                    colnames(m_ttest_result) = tt_name
                    clust_group = cbind(na.omit(cbind(chrBed, 
                      clust_wig, m_ttest_result)), cluster.id)
                    write.csv(clust_group, paste(matrix_dir, 
                      "/", chrName, "_cluster.csv", sep = ""), 
                      row.names = FALSE)
                    clust_heatmap = NULL
                    clust_order_num = NULL
                    tmp_order_num = 0
                    for (ii in 1:clust_k) {
                      clust_heatmap = rbind(clust_heatmap, clust_group[which(clust_group[, 
                        "cluster.id"] == ii), ])
                      if (ii != clust_k) {
                        clust_order_num = rbind(clust_order_num, 
                          dim(clust_heatmap)[1])
                      }
                    }
                    rownames(clust_heatmap) = NULL
                    if ((outputpdf == "TRUE") || (outputpdf == 
                      "true")) {
                      pdf(paste(matrix_dir, "/", chrName, "_cluster_heatmap.pdf", 
                        sep = ""), width = 8, height = 8)
                    }
                    else {
                      jpeg(paste(matrix_dir, "/", chrName, "_cluster_heatmap.jpeg", 
                        sep = ""), width = 1000, height = 1000, 
                        quality = 100)
                    }
                    hm_data = as.matrix(clust_heatmap[, 5:(5 + 
                      all_wig_num - 1)])
                    for (hmi in 1:dim(hm_data)[2]) {
                      hm_data[, hmi] = (hm_data[, hmi] - range(hm_data[, 
                        hmi])[1])/(range(hm_data[, hmi])[2] - 
                        range(hm_data[, hmi])[1])
                    }
                    clust_col_num = c(1:dim(hm_data)[2])
                    heatmap.2(hm_data, rowsep = clust_order_num, 
                      sepcolor = "black", sepwidth = c(0.1, 0.1), 
                      srtCol = 25, adjCol = c(0.6, 0.8), cexCol = 0.7, 
                      col = whitered, dendrogram = "none", Rowv = FALSE, 
                      Colv = FALSE, adjRow = c(-500, -500), breaks = 256, 
                      key.title = NA, key.xlab = NA, key.par = list(mgp = c(1.5, 
                        0.5, 0), mar = c(1, 2.5, 1, 0)), key.xtickfun = function() {
                        cex <- par("cex") * par("cex.axis")
                        side <- 1
                        line <- 0
                        col <- par("col.axis")
                        font <- par("font.axis")
                        mtext("low", side = side, at = 0, adj = 0, 
                          line = line, cex = cex, col = col, 
                          font = font)
                        mtext("high", side = side, at = 1, adj = 1, 
                          line = line, cex = cex, col = col, 
                          font = font)
                        return(list(labels = FALSE, tick = FALSE))
                      })
                    dev.off()
                  }
                }
                n_count = 0
                for (ii in 1:chrTotSize) {
                  if (chrBedMatrix[ii] != 0) {
                    n_count = n_count + 1
                  }
                }
                random_group = matrix(data = 0, nrow = groupNum, 
                  ncol = chrTotSize)
                random_result = matrix(data = 0, nrow = groupNum, 
                  ncol = 3)
                for (ii in 1:groupNum) {
                  tmp_site = sample(1:chrTotSize, size = n_count)
                  random_group[ii, tmp_site] = 2
                }
                tmpCmap = chrCmap
                tmpCmap[lower.tri(tmpCmap)] = NA
                for (t in 1:groupNum) {
                  random_result[t, 1] = 0
                  random_result[t, 2] = 0
                  for (ii in 1:chrTotSize) {
                    if (random_group[t, ii] > 0) {
                      b = which(tmpCmap[ii, ] > 0)
                      c = which((random_group[t, b]) > 0)
                      random_result[t, 1] = random_result[t, 
                        1] + length(c)
                      random_result[t, 2] = random_result[t, 
                        2] + length(b)
                    }
                  }
                }
                bed_count = 0
                all_count = 0
                bb_count = 0
                bb_info = NULL
                bedTOall_info = NULL
                for (ii in 1:chrTotSize) {
                  if (chrBedMatrix[ii] > 0) {
                    b = which(tmpCmap[ii, ] > 0)
                    c = which((chrBedMatrix[b]) > 0)
                    if (length(b) > 0) {
                      bedTOall_info = rbind(bedTOall_info, cbind(ii, 
                        b, t(tmpCmap[ii, b])))
                    }
                    if (length(c) > 0) {
                      bb_info = rbind(bb_info, cbind(ii, b[c], 
                        t(tmpCmap[ii, b[c]])))
                    }
                    bb_count = bb_count + length(c)
                    bed_count = bed_count + length(b)
                  }
                }
                if_test = rbind(cbind(bb_info[, 3], 1), cbind(bedTOall_info[, 
                  3], 2))
                if_test = as.data.frame(if_test)
                colnames(if_test) = c("if_value", "group")
                rownames(if_test) = c(1:(dim(if_test)[1]))
                if_test$group = as.factor(if_test$group)
                if_kruskal = kruskal.test(if_value ~ group, data = if_test)
                if_kruskalmc = kruskalmc(if_value ~ group, data = if_test, 
                  probs = 0.05)
                mult <- oneway_test(if_value ~ group, data = if_test, 
                  ytrafo = function(data) trafo(data, numeric_trafo = rank), 
                  xtrafo = function(data) trafo(data, factor_trafo = function(x) model.matrix(~x - 
                    1) \%*\% t(contrMat(table(x), "Tukey"))), teststat = "max", 
                  distribution = approximate(B = 90000))
                if_pvalue = pvalue(mult, method = "single-step")
                dif_result = t.test(random_result[, 1], mu = bb_count)
                write("", file = logfile, append = TRUE)
                write("", file = logfile, append = TRUE)
                write("", file = logfile, append = TRUE)
                write("the statistic test of interaction frequency between b2b and b2o :", 
                  file = logfile, append = TRUE)
                write("", file = logfile, append = TRUE)
                write("test name : Kruskal-Wallis rank sum test", 
                  file = logfile, append = TRUE)
                write(paste("Kruskal-Wallis chi-squared : ", 
                  if_kruskal$statistic, sep = ""), file = logfile, 
                  append = TRUE)
                write(paste("Kruskal-Wallis df : ", if_kruskal$parameter, 
                  sep = ""), file = logfile, append = TRUE)
                write(paste("Kruskal-Wallis p value : ", if_kruskal$p.value, 
                  sep = ""), file = logfile, append = TRUE)
                write("", file = logfile, append = TRUE)
                write("test name : Multiple comparison test after Kruskal-Wallis", 
                  file = logfile, append = TRUE)
                write(paste("significance level : ", if_kruskalmc$signif.level, 
                  sep = ""), file = logfile, append = TRUE)
                write(paste("observed difference  : ", if_kruskalmc$dif.com$obs.dif, 
                  sep = ""), file = logfile, append = TRUE)
                write(paste("critical difference  : ", if_kruskalmc$dif.com$critical.dif, 
                  sep = ""), file = logfile, append = TRUE)
                write(paste("exist difference  : ", if_kruskalmc$dif.com$difference, 
                  sep = ""), file = logfile, append = TRUE)
                write("", file = logfile, append = TRUE)
                write("", file = logfile, append = TRUE)
                write("", file = logfile, append = TRUE)
                write("", file = logfile, append = TRUE)
                write("the statistic test of interaction number between b2b and o2o :", 
                  file = logfile, append = TRUE)
                write("", file = logfile, append = TRUE)
                write("test name : t-test", file = logfile, append = TRUE)
                write(paste("numbers of random group : ", groupNum, 
                  sep = ""), file = logfile, append = TRUE)
                write(paste("95 percent confidence interval of random group : ", 
                  dif_result$conf.int[1], " ~ ", dif_result$conf.int[2], 
                  sep = ""), file = logfile, append = TRUE)
                write(paste("numbers of b2b : ", bb_count, sep = ""), 
                  file = logfile, append = TRUE)
                write(paste("t test p value : ", dif_result$p.value, 
                  sep = ""), file = logfile, append = TRUE)
                if (dif_result$p.value < 0.05) {
                  write("exist difference  : TRUE", file = logfile, 
                    append = TRUE)
                }
                else {
                  write("exist difference  : FALSE", file = logfile, 
                    append = TRUE)
                }
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
                print(matrix_full_dir[i])
                if (chrBedNum > 0) {
                  print(chrName)
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
                  chrBedBin = check_bed_bin(chrBed, resolution * 
                    1000)
                  chrBedMatrix = convert_bed_to_matrix(chrBed, 
                    resolution * 1000, chrName, chrTotSize, bedWindow)
                  logfile = paste(paste(matrix_dir, "/", chrName, 
                    "_statistic.txt", sep = ""))
                  if (file.exists(logfile) == TRUE) {
                    file.remove(logfile)
                  }
                  starttime = paste("Analysis start time:", as.character(Sys.time()))
                  write(starttime, file = logfile, append = TRUE)
                  if (all_wig_num > 0) {
                    clust_wig = NULL
                    all_wig_info = NULL
                    all_wig_info_loc = matrix(data = 0, nrow = all_wig_num, 
                      ncol = 2)
                    for (ww in 1:all_wig_num) {
                      chrBedWig = load_bed_wig(all_wig_file[ww], 
                        chrBed, chrName, chrstart, chrend, 0)
                      clust_wig = cbind(clust_wig, chrBedWig[, 
                        4])
                      if (is.null(all_wig_info)) {
                        wig_info_start = 1
                      }
                      else {
                        wig_info_start = dim(all_wig_info)[1] + 
                          1
                      }
                      one_wig_info = load_all_wig(all_wig_file[ww])
                      one_wig_info = one_wig_info[one_wig_info$chr == 
                        chrom, ]
                      one_wig_info[, 2] = as.numeric(one_wig_info[, 
                        2])
                      one_wig_info[, 3] = as.numeric(one_wig_info[, 
                        3])
                      one_wig_info[, 4] = as.numeric(one_wig_info[, 
                        4])
                      if (chrend > 0) {
                        one_wig_info = one_wig_info[which(one_wig_info[, 
                          3] < chrend), ]
                      }
                      if (chrstart > 0) {
                        one_wig_info[, 2] = one_wig_info[, 2] - 
                          chrstart
                        one_wig_info[, 3] = one_wig_info[, 3] - 
                          chrstart
                        one_wig_info = one_wig_info[which(one_wig_info[, 
                          2] >= 0), ]
                      }
                      wig_info_end = wig_info_start + dim(one_wig_info)[1] - 
                        1
                      all_wig_info_loc[ww, 1] = wig_info_start
                      all_wig_info_loc[ww, 2] = wig_info_end
                      all_wig_info = rbind(all_wig_info, one_wig_info)
                    }
                    if ((is.null(clust_wig)) == FALSE) {
                      tt = c(1:(dim(clust_wig)[1]))
                      clust_name = paste(tt, "_", chrBed[, 1], 
                        ":", (chrBed[, 2] + chrstart), "-", (chrBed[, 
                          3] + chrstart), sep = "")
                      row.names(clust_wig) = clust_name
                      colnames(clust_wig) = all_wig_name
                      mydata = na.omit(clust_wig)
                      suppressPackageStartupMessages(library("lattice"))
                      suppressPackageStartupMessages(library("flexclust"))
                      if (all_wig_num > 1) {
                        bcl <- bootFlexclust(mydata, k = 2:7, 
                          nboot = 50, FUN = cclust, multicore = FALSE)
                        if ((outputpdf == "TRUE") || (outputpdf == 
                          "true")) {
                          pdf(paste(matrix_dir, "/", chrName, 
                            "_cluster_k_density.pdf", sep = ""), 
                            width = 8, height = 8)
                        }
                        else {
                          jpeg(paste(matrix_dir, "/", chrName, 
                            "_cluster_k_density.jpeg", sep = ""), 
                            width = 1000, height = 1000, quality = 100)
                        }
                        plot(bcl)
                        densityplot(bcl, from = 0)
                        dev.off()
                      }
                      out.dist = dist(mydata, method = dist_method)
                      out.hclust = hclust(out.dist, method = clust_method)
                      if ((outputpdf == "TRUE") || (outputpdf == 
                        "true")) {
                        pdf(paste(matrix_dir, "/", chrName, "_cluster_tree.pdf", 
                          sep = ""), width = 8, height = 8)
                      }
                      else {
                        jpeg(paste(matrix_dir, "/", chrName, 
                          "_cluster_tree.jpeg", sep = ""), width = 1000, 
                          height = 1000, quality = 100)
                      }
                      plot(out.hclust)
                      rect.hclust(out.hclust, clust_k)
                      cluster.id = cutree(out.hclust, clust_k)
                      dev.off()
                      row.names(chrBed) = clust_name
                      m_ttest_result = NULL
                      for (ww in 1:all_wig_num) {
                        tmp_ttest_result = data.frame(wig1_pvalue = numeric(dim(clust_wig)[1]), 
                          wig1_difference = character(dim(clust_wig)[1]), 
                          stringsAsFactors = FALSE)
                        tmp_wig_ttest = NULL
                        write("", file = logfile, append = TRUE)
                        write("", file = logfile, append = TRUE)
                        write("", file = logfile, append = TRUE)
                        for (iii in 1:(dim(clust_wig)[1])) {
                          if (is.na(clust_wig[iii, 1]) == FALSE) {
                            tmp_wig_ttest = rbind(tmp_wig_ttest, 
                              t.test(all_wig_info[(all_wig_info_loc[ww, 
                                1]:all_wig_info_loc[ww, 2]), 
                                4], mu = clust_wig[iii, 1]))
                          }
                          else {
                            tmp_wig_ttest = rbind(tmp_wig_ttest, 
                              t.test(all_wig_info[(all_wig_info_loc[ww, 
                                1]:all_wig_info_loc[ww, 2]), 
                                4], mu = 0))
                          }
                        }
                        m_wig_equal = which(tmp_wig_ttest[, 3] > 
                          0.05)
                        m_wig_mean = mean(all_wig_info[(all_wig_info_loc[ww, 
                          1]:all_wig_info_loc[ww, 2]), 4])
                        m_wig_nequal = which(tmp_wig_ttest[, 
                          3] <= 0.05)
                        m_wig_more = m_wig_nequal[which(clust_wig[m_wig_nequal[], 
                          1] > m_wig_mean)]
                        m_wig_less = m_wig_nequal[which(clust_wig[m_wig_nequal[], 
                          1] <= m_wig_mean)]
                        tmp_ttest_result[, 1] = as.data.frame(as.matrix(tmp_wig_ttest[, 
                          3]))
                        tmp_ttest_result[, 1] = as.numeric(tmp_ttest_result[, 
                          1])
                        tmp_ttest_result[m_wig_more[], 2] = "more"
                        tmp_ttest_result[m_wig_less[], 2] = "less"
                        tmp_ttest_result[m_wig_equal[], 2] = "equal"
                        if (is.null(m_ttest_result)) {
                          m_ttest_result = tmp_ttest_result
                        }
                        else {
                          m_ttest_result = cbind(m_ttest_result, 
                            tmp_ttest_result)
                        }
                        wig_test = rbind(cbind(na.omit(clust_wig[, 
                          1]), 1), cbind(na.omit(all_wig_info[(all_wig_info_loc[ww, 
                          1]:all_wig_info_loc[ww, 2]), 4]), 2))
                        wig_test = as.data.frame(wig_test)
                        colnames(wig_test) = c("wig_value", "group")
                        rownames(wig_test) = c(1:(dim(wig_test)[1]))
                        wig_test$group = as.factor(wig_test$group)
                        wig_kruskal = kruskal.test(wig_value ~ 
                          group, data = wig_test)
                        wig_kruskalmc = kruskalmc(wig_value ~ 
                          group, data = wig_test, probs = 0.05)
                        wig_mult <- oneway_test(wig_value ~ group, 
                          data = wig_test, ytrafo = function(data) trafo(data, 
                            numeric_trafo = rank), xtrafo = function(data) trafo(data, 
                            factor_trafo = function(x) model.matrix(~x - 
                              1) \%*\% t(contrMat(table(x), "Tukey"))), 
                          teststat = "max", distribution = approximate(B = 90000))
                        wig_pvalue = pvalue(wig_mult, method = "single-step")
                        write(paste("the statistic test between BED WIG and  WIG : ", 
                          all_wig_name[ww], sep = ""), file = logfile, 
                          append = TRUE)
                        write("", file = logfile, append = TRUE)
                        write("test name : Kruskal-Wallis rank sum test", 
                          file = logfile, append = TRUE)
                        write(paste("Kruskal-Wallis chi-squared : ", 
                          wig_kruskal$statistic, sep = ""), file = logfile, 
                          append = TRUE)
                        write(paste("Kruskal-Wallis df : ", wig_kruskal$parameter, 
                          sep = ""), file = logfile, append = TRUE)
                        write(paste("Kruskal-Wallis p value : ", 
                          wig_kruskal$p.value, sep = ""), file = logfile, 
                          append = TRUE)
                        write("", file = logfile, append = TRUE)
                        write("test name : Multiple comparison test after Kruskal-Wallis", 
                          file = logfile, append = TRUE)
                        write(paste("significance level : ", 
                          wig_kruskalmc$signif.level, sep = ""), 
                          file = logfile, append = TRUE)
                        write(paste("observed difference  : ", 
                          wig_kruskalmc$dif.com$obs.dif, sep = ""), 
                          file = logfile, append = TRUE)
                        write(paste("critical difference  : ", 
                          wig_kruskalmc$dif.com$critical.dif, 
                          sep = ""), file = logfile, append = TRUE)
                        write(paste("exist difference  : ", wig_kruskalmc$dif.com$difference, 
                          sep = ""), file = logfile, append = TRUE)
                        write("", file = logfile, append = TRUE)
                        write("", file = logfile, append = TRUE)
                        write("", file = logfile, append = TRUE)
                      }
                      tt_name = NULL
                      for (ww in 1:all_wig_num) {
                        tt_name = c(tt_name, paste(all_wig_name[ww], 
                          "_pvalue", sep = ""), paste(all_wig_name[ww], 
                          "_difference", sep = ""))
                      }
                      colnames(m_ttest_result) = tt_name
                      clust_group = cbind(na.omit(cbind(chrBed, 
                        clust_wig, m_ttest_result)), cluster.id)
                      write.csv(clust_group, paste(matrix_dir, 
                        "/", chrName, "_cluster.csv", sep = ""), 
                        row.names = FALSE)
                      clust_heatmap = NULL
                      clust_order_num = NULL
                      tmp_order_num = 0
                      for (ii in 1:clust_k) {
                        clust_heatmap = rbind(clust_heatmap, 
                          clust_group[which(clust_group[, "cluster.id"] == 
                            ii), ])
                        if (ii != clust_k) {
                          clust_order_num = rbind(clust_order_num, 
                            dim(clust_heatmap)[1])
                        }
                      }
                      rownames(clust_heatmap) = NULL
                      if ((outputpdf == "TRUE") || (outputpdf == 
                        "true")) {
                        pdf(paste(matrix_dir, "/", chrName, "_cluster_heatmap.pdf", 
                          sep = ""), width = 8, height = 8)
                      }
                      else {
                        jpeg(paste(matrix_dir, "/", chrName, 
                          "_cluster_heatmap.jpeg", sep = ""), 
                          width = 1000, height = 1000, quality = 100)
                      }
                      hm_data = as.matrix(clust_heatmap[, 5:(5 + 
                        all_wig_num - 1)])
                      for (hmi in 1:dim(hm_data)[2]) {
                        hm_data[, hmi] = (hm_data[, hmi] - range(hm_data[, 
                          hmi])[1])/(range(hm_data[, hmi])[2] - 
                          range(hm_data[, hmi])[1])
                      }
                      clust_col_num = c(1:dim(hm_data)[2])
                      heatmap.2(hm_data, rowsep = clust_order_num, 
                        sepcolor = "black", sepwidth = c(0.1, 
                          0.1), srtCol = 25, adjCol = c(0.6, 
                          0.8), cexCol = 0.7, col = whitered, 
                        dendrogram = "none", Rowv = FALSE, Colv = FALSE, 
                        adjRow = c(-500, -500), breaks = 256, 
                        key.title = NA, key.xlab = NA, key.par = list(mgp = c(1.5, 
                          0.5, 0), mar = c(1, 2.5, 1, 0)), key.xtickfun = function() {
                          cex <- par("cex") * par("cex.axis")
                          side <- 1
                          line <- 0
                          col <- par("col.axis")
                          font <- par("font.axis")
                          mtext("low", side = side, at = 0, adj = 0, 
                            line = line, cex = cex, col = col, 
                            font = font)
                          mtext("high", side = side, at = 1, 
                            adj = 1, line = line, cex = cex, 
                            col = col, font = font)
                          return(list(labels = FALSE, tick = FALSE))
                        })
                      dev.off()
                    }
                  }
                  n_count = 0
                  for (ii in 1:chrTotSize) {
                    if (chrBedMatrix[ii] != 0) {
                      n_count = n_count + 1
                    }
                  }
                  random_group = matrix(data = 0, nrow = groupNum, 
                    ncol = chrTotSize)
                  random_result = matrix(data = 0, nrow = groupNum, 
                    ncol = 3)
                  for (ii in 1:groupNum) {
                    tmp_site = sample(1:chrTotSize, size = n_count)
                    random_group[ii, tmp_site] = 2
                  }
                  tmpCmap = chrCmap
                  tmpCmap[lower.tri(tmpCmap)] = NA
                  for (t in 1:groupNum) {
                    random_result[t, 1] = 0
                    random_result[t, 2] = 0
                    for (ii in 1:chrTotSize) {
                      if (random_group[t, ii] > 0) {
                        b = which(tmpCmap[ii, ] > 0)
                        c = which((random_group[t, b]) > 0)
                        random_result[t, 1] = random_result[t, 
                          1] + length(c)
                        random_result[t, 2] = random_result[t, 
                          2] + length(b)
                      }
                    }
                  }
                  bed_count = 0
                  all_count = 0
                  bb_count = 0
                  bb_info = NULL
                  bedTOall_info = NULL
                  for (ii in 1:chrTotSize) {
                    if (chrBedMatrix[ii] > 0) {
                      b = which(tmpCmap[ii, ] > 0)
                      c = which((chrBedMatrix[b]) > 0)
                      if (length(b) > 0) {
                        bedTOall_info = rbind(bedTOall_info, 
                          cbind(ii, b, t(tmpCmap[ii, b])))
                      }
                      if (length(c) > 0) {
                        bb_info = rbind(bb_info, cbind(ii, b[c], 
                          t(tmpCmap[ii, b[c]])))
                      }
                      bb_count = bb_count + length(c)
                      bed_count = bed_count + length(b)
                    }
                  }
                  if_test = rbind(cbind(bb_info[, 3], 1), cbind(bedTOall_info[, 
                    3], 2))
                  if_test = as.data.frame(if_test)
                  colnames(if_test) = c("if_value", "group")
                  rownames(if_test) = c(1:(dim(if_test)[1]))
                  if_test$group = as.factor(if_test$group)
                  if_kruskal = kruskal.test(if_value ~ group, 
                    data = if_test)
                  if_kruskalmc = kruskalmc(if_value ~ group, 
                    data = if_test, probs = 0.05)
                  mult <- oneway_test(if_value ~ group, data = if_test, 
                    ytrafo = function(data) trafo(data, numeric_trafo = rank), 
                    xtrafo = function(data) trafo(data, factor_trafo = function(x) model.matrix(~x - 
                      1) \%*\% t(contrMat(table(x), "Tukey"))), 
                    teststat = "max", distribution = approximate(B = 90000))
                  if_pvalue = pvalue(mult, method = "single-step")
                  dif_result = t.test(random_result[, 1], mu = bb_count)
                  write("", file = logfile, append = TRUE)
                  write("", file = logfile, append = TRUE)
                  write("", file = logfile, append = TRUE)
                  write("the statistic test of interaction frequency between b2b and b2o :", 
                    file = logfile, append = TRUE)
                  write("", file = logfile, append = TRUE)
                  write("test name : Kruskal-Wallis rank sum test", 
                    file = logfile, append = TRUE)
                  write(paste("Kruskal-Wallis chi-squared : ", 
                    if_kruskal$statistic, sep = ""), file = logfile, 
                    append = TRUE)
                  write(paste("Kruskal-Wallis df : ", if_kruskal$parameter, 
                    sep = ""), file = logfile, append = TRUE)
                  write(paste("Kruskal-Wallis p value : ", if_kruskal$p.value, 
                    sep = ""), file = logfile, append = TRUE)
                  write("", file = logfile, append = TRUE)
                  write("test name : Multiple comparison test after Kruskal-Wallis", 
                    file = logfile, append = TRUE)
                  write(paste("significance level : ", if_kruskalmc$signif.level, 
                    sep = ""), file = logfile, append = TRUE)
                  write(paste("observed difference  : ", if_kruskalmc$dif.com$obs.dif, 
                    sep = ""), file = logfile, append = TRUE)
                  write(paste("critical difference  : ", if_kruskalmc$dif.com$critical.dif, 
                    sep = ""), file = logfile, append = TRUE)
                  write(paste("exist difference  : ", if_kruskalmc$dif.com$difference, 
                    sep = ""), file = logfile, append = TRUE)
                  write("", file = logfile, append = TRUE)
                  write("", file = logfile, append = TRUE)
                  write("", file = logfile, append = TRUE)
                  write("", file = logfile, append = TRUE)
                  write("the statistic test of interaction number between b2b and o2o :", 
                    file = logfile, append = TRUE)
                  write("", file = logfile, append = TRUE)
                  write("test name : t-test", file = logfile, 
                    append = TRUE)
                  write(paste("numbers of random group : ", groupNum, 
                    sep = ""), file = logfile, append = TRUE)
                  write(paste("95 percent confidence interval of random group : ", 
                    dif_result$conf.int[1], " ~ ", dif_result$conf.int[2], 
                    sep = ""), file = logfile, append = TRUE)
                  write(paste("numbers of b2b : ", bb_count, 
                    sep = ""), file = logfile, append = TRUE)
                  write(paste("t test p value : ", dif_result$p.value, 
                    sep = ""), file = logfile, append = TRUE)
                  if (dif_result$p.value < 0.05) {
                    write("exist difference  : TRUE", file = logfile, 
                      append = TRUE)
                  }
                  else {
                    write("exist difference  : FALSE", file = logfile, 
                      append = TRUE)
                  }
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
