\name{network_analysis}
\alias{network_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
to plot interaction network picture and make clusters
}
\description{
to plot interaction network picture and make clusters
}
\usage{
network_analysis(bedFile, matrix_dir = "hg19", outputpdf = "FALSE", chrom = "all", chrstart = 0, chrend = 0, resolution = 100, bedWindow = 0, net_layout = "layout.fruchterman.reingold", netplot = "TRUE", NetClusterType = "multileve", NetVertexSize = 2, NetVertexChangeSize = "degree", NetVertexLableDist = 0.1, NetVertexColor = "#7fbc41", NetVertexLabelCex = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bedFile}{
%%     ~~Describe \code{bedFile} here~~
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
  \item{bedWindow}{
%%     ~~Describe \code{bedWindow} here~~
}
  \item{net_layout}{
%%     ~~Describe \code{net_layout} here~~
}
  \item{netplot}{
%%     ~~Describe \code{netplot} here~~
}
  \item{NetClusterType}{
%%     ~~Describe \code{NetClusterType} here~~
}
  \item{NetVertexSize}{
%%     ~~Describe \code{NetVertexSize} here~~
}
  \item{NetVertexChangeSize}{
%%     ~~Describe \code{NetVertexChangeSize} here~~
}
  \item{NetVertexLableDist}{
%%     ~~Describe \code{NetVertexLableDist} here~~
}
  \item{NetVertexColor}{
%%     ~~Describe \code{NetVertexColor} here~~
}
  \item{NetVertexLabelCex}{
%%     ~~Describe \code{NetVertexLabelCex} here~~
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
function (bedFile, matrix_dir = "hg19", outputpdf = "FALSE", 
    chrom = "all", chrstart = 0, chrend = 0, resolution = 100, 
    bedWindow = 0, net_layout = "layout.fruchterman.reingold", 
    netplot = "TRUE", NetClusterType = "multileve", NetVertexSize = 2, 
    NetVertexChangeSize = "degree", NetVertexLableDist = 0.1, 
    NetVertexColor = "#7fbc41", NetVertexLabelCex = 3) 
{
    gc()
    source("calbed.R")
    matrix_name_dir = list.files(path = matrix_dir, full.names = F, 
        pattern = ".matrix")
    matrix_full_dir = list.files(path = matrix_dir, full.names = T, 
        pattern = ".matrix")
    m_bed = load_bed(bedFile)
    chrNum = length(matrix_name_dir)
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
                chrBedToBedInter = find_bed_to_bed_interaction(chrCmap, 
                  chrBedMatrix, chrBedBin, chrBed, chrName, chrCmap)
                if ((netplot == "true") || (netplot == "TRUE")) {
                  netgraph = data.frame(p1 = character(dim(chrBedToBedInter)[1]), 
                    p2 = character(dim(chrBedToBedInter)[1]), 
                    weight = numeric(dim(chrBedToBedInter)[1]))
                  netgraph[, 1] = as.data.frame(paste(chrBedToBedInter[, 
                    2], ":", chrBedToBedInter[, 3], "-", chrBedToBedInter[, 
                    4], sep = ""))
                  netgraph[, 2] = as.data.frame(paste(chrBedToBedInter[, 
                    8], ":", chrBedToBedInter[, 9], "-", chrBedToBedInter[, 
                    10], sep = ""))
                  netgraph[, 3] = as.data.frame(chrBedToBedInter[, 
                    13])
                  set.seed(1234)
                  if ((outputpdf == "TRUE") || (outputpdf == 
                    "true")) {
                    pdf(paste(matrix_dir, "/", chrName, "_netplot.pdf", 
                      sep = ""), width = 8, height = 8)
                  }
                  else {
                    jpeg(paste(matrix_dir, "/", chrName, "_netplot.jpeg", 
                      sep = ""), width = 1000, height = 1000, 
                      quality = 100)
                  }
                  g = graph.data.frame(netgraph, directed = F)
                  set.seed(1234)
                  netnodename = names(V(g))
                  netnodechrnum = regexpr(":", netnodename)
                  netnodechr = substr(netnodename, 1, netnodechrnum - 
                    1)
                  netnodestartnum = regexpr("-", netnodename)
                  netnodestart = as.numeric(substr(netnodename, 
                    netnodechrnum + 1, netnodestartnum - 1))
                  netnodeendnum = nchar(netnodename)
                  netnodeend = as.numeric(substr(netnodename, 
                    netnodestartnum + 1, netnodeendnum))
                  netcsv = data.frame(chrom = character(0), start = numeric(0), 
                    end = numeric(0), degree = numeric(0), closeness = numeric(0), 
                    betweenness = numeric(0), Local_cluster_coefficient = numeric(0), 
                    Eigenvector_centrality = numeric(0), membership = numeric(0), 
                    stringsAsFactors = FALSE)
                  netdegree = degree(g)
                  netcloseness = closeness(g)
                  netbetweenness = betweenness(g)
                  netcoefficient = transitivity(g, type = "local")
                  netcentrality = evcent(g)$vector
                  netcsv[1:(length(netdegree)), 1] = netnodechr
                  netcsv[1:(length(netdegree)), 2] = netnodestart
                  netcsv[1:(length(netdegree)), 3] = netnodeend
                  netcsv[1:(length(netdegree)), 4] = as.data.frame(netdegree)
                  netcsv[1:(length(netdegree)), 5] = as.data.frame(netcloseness)
                  netcsv[1:(length(netdegree)), 6] = as.data.frame(netbetweenness)
                  netcsv[1:(length(netdegree)), 7] = as.data.frame(netcoefficient)
                  netcsv[1:(length(netdegree)), 8] = as.data.frame(netcentrality)
                  colors <- c("#fff7f3", "#fff7f3", "#fde0dd", 
                    "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", 
                    "#b0017e", "#7a0177", "#49006a")
                  weight_range = range(E(g)$weight)
                  E(g)$color = colors[1]
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1])/10)]$color = colors[2]
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 2/10)]$color = colors[3]
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 3/10)]$color = colors[4]
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 4/10)]$color = colors[5]
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 5/10)]$color = colors[6]
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 6/10)]$color = colors[7]
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 7/10)]$color = colors[8]
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 8/10)]$color = colors[9]
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 9/10)]$color = colors[10]
                  edge_width = 0.05
                  E(g)$width = edge_width
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1])/10)]$width = 2 * edge_width
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 2/10)]$width = 3 * edge_width
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 3/10)]$width = 4 * edge_width
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 4/10)]$width = 5 * edge_width
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 5/10)]$width = 6 * edge_width
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 6/10)]$width = 7 * edge_width
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 7/10)]$width = 8 * edge_width
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 8/10)]$width = 9 * edge_width
                  E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                    weight_range[1]) * 9/10)]$width = 10 * edge_width
                  if (NetVertexChangeSize == "degree") {
                    V(g)$deg <- netcsv[, 4]
                    deg_range = range(netcsv[, 4])
                    V(g)$size = NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1])/5)]$size = 2 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                  }
                  else if (NetVertexChangeSize == "closeness") {
                    V(g)$deg <- netcsv[, 5]
                    deg_range = range(netcsv[, 5])
                    V(g)$size = NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1])/5)]$size = 2 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                  }
                  else if (NetVertexChangeSize == "betweenness") {
                    V(g)$deg <- netcsv[, 6]
                    deg_range = range(netcsv[, 6])
                    V(g)$size = NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1])/5)]$size = 2 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                  }
                  else if (NetVertexChangeSize == "Local_cluster_coefficient") {
                    V(g)$deg <- netcsv[, 7]
                    deg_range = range(netcsv[, 7])
                    V(g)$size = NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1])/5)]$size = 2 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                  }
                  else if (NetVertexChangeSize == "Eigenvector_centrality") {
                    V(g)$deg <- netcsv[, 8]
                    deg_range = range(netcsv[, 8])
                    V(g)$size = NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1])/5)]$size = 2 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                    V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                      deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                  }
                  else {
                    V(g)$size = NetVertexSize
                  }
                  if (net_layout == "layout.fruchterman.reingold") {
                    if (NetClusterType == "NULL") {
                      plot(g, layout = layout.fruchterman.reingold, 
                        vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                      netcsv[1:(length(netdegree)), 9] = "NULL"
                    }
                    if (NetClusterType == "edgeBetweenness") {
                      system.time(ec <- edge.betweenness.community(g))
                      print(modularity(ec))
                      netcsv[1:(length(netdegree)), 9] = ec$membership
                      plot(ec, g, layout = layout.fruchterman.reingold, 
                        vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                    }
                    if (NetClusterType == "walktrap") {
                      system.time(wc <- walktrap.community(g))
                      netcsv[1:(length(netdegree)), 9] = wc$membership
                      print(modularity(wc))
                      plot(wc, g, layout = layout.fruchterman.reingold, 
                        vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                    }
                    if (NetClusterType == "multileve") {
                      system.time(mc <- multilevel.community(g, 
                        weights = NA))
                      print(modularity(mc))
                      plot(mc, g, layout = layout.fruchterman.reingold, 
                        vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                      netcsv[1:(length(netdegree)), 9] = mc$membership
                    }
                    if (NetClusterType == "labelPropagation") {
                      system.time(lc <- label.propagation.community(g))
                      print(modularity(lc))
                      plot(lc, g, layout = layout.fruchterman.reingold, 
                        vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                      netcsv[1:(length(netdegree)), 9] = lc$membership
                    }
                  }
                  else if (net_layout == "layout.circle") {
                    if (NetClusterType == "NULL") {
                      plot(g, layout = layout.circle, vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                      netcsv[1:(length(netdegree)), 9] = "NULL"
                    }
                    if (NetClusterType == "edgeBetweenness") {
                      system.time(ec <- edge.betweenness.community(g))
                      print(modularity(ec))
                      netcsv[1:(length(netdegree)), 9] = ec$membership
                      plot(ec, g, layout = layout.circle, vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                    }
                    if (NetClusterType == "walktrap") {
                      system.time(wc <- walktrap.community(g))
                      netcsv[1:(length(netdegree)), 9] = wc$membership
                      print(modularity(wc))
                      plot(wc, g, layout = layout.circle, vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                    }
                    if (NetClusterType == "multileve") {
                      system.time(mc <- multilevel.community(g, 
                        weights = NA))
                      print(modularity(mc))
                      plot(mc, g, layout = layout.circle, vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                      netcsv[1:(length(netdegree)), 9] = mc$membership
                    }
                    if (NetClusterType == "labelPropagation") {
                      system.time(lc <- label.propagation.community(g))
                      print(modularity(lc))
                      plot(lc, g, layout = layout.circle, vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                      netcsv[1:(length(netdegree)), 9] = lc$membership
                    }
                  }
                  else {
                    if (NetClusterType == "NULL") {
                      plot(g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                      netcsv[1:(length(netdegree)), 9] = "NULL"
                    }
                    if (NetClusterType == "edgeBetweenness") {
                      system.time(ec <- edge.betweenness.community(g))
                      print(modularity(ec))
                      netcsv[1:(length(netdegree)), 9] = ec$membership
                      plot(ec, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                    }
                    if (NetClusterType == "walktrap") {
                      system.time(wc <- walktrap.community(g))
                      netcsv[1:(length(netdegree)), 9] = wc$membership
                      print(modularity(wc))
                      plot(wc, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                    }
                    if (NetClusterType == "multileve") {
                      system.time(mc <- multilevel.community(g, 
                        weights = NA))
                      print(modularity(mc))
                      plot(mc, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                      netcsv[1:(length(netdegree)), 9] = mc$membership
                    }
                    if (NetClusterType == "labelPropagation") {
                      system.time(lc <- label.propagation.community(g))
                      print(modularity(lc))
                      plot(lc, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                        vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                        vertex.label.cex = NetVertexLabelCex)
                      netcsv[1:(length(netdegree)), 9] = lc$membership
                    }
                  }
                  dev.off()
                  write.csv(netcsv, paste(matrix_dir, "/", chrName, 
                    "_network.csv", sep = ""), row.names = FALSE)
                }
                bedIplot = cbind(rbind(chrBedToBedInter[, 14], 
                  chrBedToBedInter[, 15]), rbind(chrBedToBedInter[, 
                  15], chrBedToBedInter[, 14]))
                hm_dim = dim(chrCmap)[1]
                chrCmap = as.matrix(chrCmap)
                hm_mean = mean(chrCmap)
                for (j in 1:hm_dim) {
                  chrCmap[j, which(chrCmap[j, ] > 5 * hm_mean)] = 5 * 
                    hm_mean
                }
                chrhmCmap = melt(chrCmap)
                print(paste("plot ", chrName, "bed picture", 
                  sep = ""))
                if ((outputpdf == "TRUE") || (outputpdf == "true")) {
                  pdf(paste(matrix_dir, "/", chrName, "_bedplot.pdf", 
                    sep = ""), width = 8, height = 8)
                }
                else {
                  jpeg(paste(matrix_dir, "/", chrName, "_bedplot.jpeg", 
                    sep = ""), width = 1000, height = 1000, quality = 100)
                }
                grid.newpage()
                heatmapViewport <- viewport(height = 0.5, width = 0.5, 
                  x = 0.25, y = 0.5)
                scatterViewport <- viewport(height = 0.5, width = 0.5, 
                  x = 0.75, y = 0.5)
                densityViewport <- viewport(height = 0.25, width = 0.5, 
                  x = 0.75, y = 0.125)
                hmdensityViewport <- viewport(height = 0.25, 
                  width = 0.5, x = 0.25, y = 0.125)
                jit = position_jitter(width = 0.5)
                hmrange = range(chrhmCmap[, 1])
                bedIplot = bedIplot + hmrange[1]
                chrhm = ggplot(chrhmCmap, aes(x = Var1, y = Var2, 
                  fill = value)) + scale_y_discrete(breaks = seq(0, 
                  10, 5)) + xlab("chrom") + ylab("chrom") + scale_fill_gradient(low = "white", 
                  high = "red") + geom_tile() + guides(fill = FALSE)
                chrbedplot = qplot(bedIplot[1, ], bedIplot[2, 
                  ], alpha = I(1/10), size = I(1)) + xlab("chrom") + 
                  ylab("chrom") + geom_jitter(position = jit, 
                  colour = "black", alpha = 1/100)
                chrbeddensitydata = as.data.frame(c(chrBedToBedInter[, 
                  14], chrBedToBedInter[, 15]))
                colnames(chrbeddensitydata) = "bed"
                chrbeddensity = ggplot(chrbeddensitydata) + geom_density(aes(x = bed))
                chrhmdensitydata = NULL
                for (iiii in 1:chrTotSize) {
                  chrhmdensitydata = c(chrhmdensitydata, which(chrCmap[, 
                    iiii] > 0))
                }
                chrhmdensitydata = as.data.frame(chrhmdensitydata)
                colnames(chrhmdensitydata) = "chrom"
                chrcmapdensity = ggplot(chrhmdensitydata) + geom_density(aes(x = chrom))
                print(chrhm, vp = heatmapViewport)
                print(chrbedplot, vp = scatterViewport)
                print(chrbeddensity, vp = densityViewport)
                print(chrcmapdensity, vp = hmdensityViewport)
                dev.off()
                write.table(chrBedToBedInter, file = paste(matrix_dir, 
                  "/", chrName, "_BedToBedInter.txt", sep = ""), 
                  sep = "\t", row.names = FALSE, col.names = FALSE, 
                  quote = FALSE)
                rm(chrCmap)
                rm(chrBed)
                rm(chrBedBin)
                rm(chrBedMatrix)
                rm(chrBedToBedInter)
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
                if (chrend <= 0) {
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
                      resolution * 1000, chrName, chrTotSize, 
                      bedWindow)
                    chrBedToBedInter = find_bed_to_bed_interaction(chrCmap, 
                      chrBedMatrix, chrBedBin, chrBed, chrName, 
                      chrCmap)
                    if ((netplot == "true") || (netplot == "TRUE")) {
                      netgraph = data.frame(p1 = character(dim(chrBedToBedInter)[1]), 
                        p2 = character(dim(chrBedToBedInter)[1]), 
                        weight = numeric(dim(chrBedToBedInter)[1]))
                      netgraph[, 1] = as.data.frame(paste(chrBedToBedInter[, 
                        2], ":", chrBedToBedInter[, 3], "-", 
                        chrBedToBedInter[, 4], sep = ""))
                      netgraph[, 2] = as.data.frame(paste(chrBedToBedInter[, 
                        8], ":", chrBedToBedInter[, 9], "-", 
                        chrBedToBedInter[, 10], sep = ""))
                      netgraph[, 3] = as.data.frame(chrBedToBedInter[, 
                        13])
                      set.seed(1234)
                      if ((outputpdf == "TRUE") || (outputpdf == 
                        "true")) {
                        pdf(paste(matrix_dir, "/", chrName, "_netplot.pdf", 
                          sep = ""), width = 8, height = 8)
                      }
                      else {
                        jpeg(paste(matrix_dir, "/", chrName, 
                          "_netplot.jpeg", sep = ""), width = 1000, 
                          height = 1000, quality = 100)
                      }
                      g = graph.data.frame(netgraph, directed = F)
                      set.seed(1234)
                      netnodename = names(V(g))
                      netnodechrnum = regexpr(":", netnodename)
                      netnodechr = substr(netnodename, 1, netnodechrnum - 
                        1)
                      netnodestartnum = regexpr("-", netnodename)
                      netnodestart = as.numeric(substr(netnodename, 
                        netnodechrnum + 1, netnodestartnum - 
                          1))
                      netnodeendnum = nchar(netnodename)
                      netnodeend = as.numeric(substr(netnodename, 
                        netnodestartnum + 1, netnodeendnum))
                      netcsv = data.frame(chrom = character(0), 
                        start = numeric(0), end = numeric(0), 
                        degree = numeric(0), closeness = numeric(0), 
                        betweenness = numeric(0), Local_cluster_coefficient = numeric(0), 
                        Eigenvector_centrality = numeric(0), 
                        membership = numeric(0), stringsAsFactors = FALSE)
                      netdegree = degree(g)
                      netcloseness = closeness(g)
                      netbetweenness = betweenness(g)
                      netcoefficient = transitivity(g, type = "local")
                      netcentrality = evcent(g)$vector
                      netcsv[1:(length(netdegree)), 1] = netnodechr
                      netcsv[1:(length(netdegree)), 2] = netnodestart
                      netcsv[1:(length(netdegree)), 3] = netnodeend
                      netcsv[1:(length(netdegree)), 4] = as.data.frame(netdegree)
                      netcsv[1:(length(netdegree)), 5] = as.data.frame(netcloseness)
                      netcsv[1:(length(netdegree)), 6] = as.data.frame(netbetweenness)
                      netcsv[1:(length(netdegree)), 7] = as.data.frame(netcoefficient)
                      netcsv[1:(length(netdegree)), 8] = as.data.frame(netcentrality)
                      colors <- c("#fff7f3", "#fff7f3", "#fde0dd", 
                        "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", 
                        "#b0017e", "#7a0177", "#49006a")
                      weight_range = range(E(g)$weight)
                      E(g)$color = colors[1]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1])/10)]$color = colors[2]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 2/10)]$color = colors[3]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 3/10)]$color = colors[4]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 4/10)]$color = colors[5]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 5/10)]$color = colors[6]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 6/10)]$color = colors[7]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 7/10)]$color = colors[8]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 8/10)]$color = colors[9]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 9/10)]$color = colors[10]
                      edge_width = 0.05
                      E(g)$width = edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1])/10)]$width = 2 * edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 2/10)]$width = 3 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 3/10)]$width = 4 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 4/10)]$width = 5 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 5/10)]$width = 6 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 6/10)]$width = 7 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 7/10)]$width = 8 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 8/10)]$width = 9 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 9/10)]$width = 10 * 
                        edge_width
                      if (NetVertexChangeSize == "degree") {
                        V(g)$deg <- netcsv[, 4]
                        deg_range = range(netcsv[, 4])
                        V(g)$size = NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1])/5)]$size = 2 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                      }
                      else if (NetVertexChangeSize == "closeness") {
                        V(g)$deg <- netcsv[, 5]
                        deg_range = range(netcsv[, 5])
                        V(g)$size = NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1])/5)]$size = 2 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                      }
                      else if (NetVertexChangeSize == "betweenness") {
                        V(g)$deg <- netcsv[, 6]
                        deg_range = range(netcsv[, 6])
                        V(g)$size = NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1])/5)]$size = 2 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                      }
                      else if (NetVertexChangeSize == "Local_cluster_coefficient") {
                        V(g)$deg <- netcsv[, 7]
                        deg_range = range(netcsv[, 7])
                        V(g)$size = NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1])/5)]$size = 2 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                      }
                      else if (NetVertexChangeSize == "Eigenvector_centrality") {
                        V(g)$deg <- netcsv[, 8]
                        deg_range = range(netcsv[, 8])
                        V(g)$size = NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1])/5)]$size = 2 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                      }
                      else {
                        V(g)$size = NetVertexSize
                      }
                      if (net_layout == "layout.fruchterman.reingold") {
                        if (NetClusterType == "NULL") {
                          plot(g, layout = layout.fruchterman.reingold, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = "NULL"
                        }
                        if (NetClusterType == "edgeBetweenness") {
                          system.time(ec <- edge.betweenness.community(g))
                          print(modularity(ec))
                          netcsv[1:(length(netdegree)), 9] = ec$membership
                          plot(ec, g, layout = layout.fruchterman.reingold, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "walktrap") {
                          system.time(wc <- walktrap.community(g))
                          netcsv[1:(length(netdegree)), 9] = wc$membership
                          print(modularity(wc))
                          plot(wc, g, layout = layout.fruchterman.reingold, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "multileve") {
                          system.time(mc <- multilevel.community(g, 
                            weights = NA))
                          print(modularity(mc))
                          plot(mc, g, layout = layout.fruchterman.reingold, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = mc$membership
                        }
                        if (NetClusterType == "labelPropagation") {
                          system.time(lc <- label.propagation.community(g))
                          print(modularity(lc))
                          plot(lc, g, layout = layout.fruchterman.reingold, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = lc$membership
                        }
                      }
                      else if (net_layout == "layout.circle") {
                        if (NetClusterType == "NULL") {
                          plot(g, layout = layout.circle, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = "NULL"
                        }
                        if (NetClusterType == "edgeBetweenness") {
                          system.time(ec <- edge.betweenness.community(g))
                          print(modularity(ec))
                          netcsv[1:(length(netdegree)), 9] = ec$membership
                          plot(ec, g, layout = layout.circle, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "walktrap") {
                          system.time(wc <- walktrap.community(g))
                          netcsv[1:(length(netdegree)), 9] = wc$membership
                          print(modularity(wc))
                          plot(wc, g, layout = layout.circle, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "multileve") {
                          system.time(mc <- multilevel.community(g, 
                            weights = NA))
                          print(modularity(mc))
                          plot(mc, g, layout = layout.circle, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = mc$membership
                        }
                        if (NetClusterType == "labelPropagation") {
                          system.time(lc <- label.propagation.community(g))
                          print(modularity(lc))
                          plot(lc, g, layout = layout.circle, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = lc$membership
                        }
                      }
                      else {
                        if (NetClusterType == "NULL") {
                          plot(g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = "NULL"
                        }
                        if (NetClusterType == "edgeBetweenness") {
                          system.time(ec <- edge.betweenness.community(g))
                          print(modularity(ec))
                          netcsv[1:(length(netdegree)), 9] = ec$membership
                          plot(ec, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "walktrap") {
                          system.time(wc <- walktrap.community(g))
                          netcsv[1:(length(netdegree)), 9] = wc$membership
                          print(modularity(wc))
                          plot(wc, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "multileve") {
                          system.time(mc <- multilevel.community(g, 
                            weights = NA))
                          print(modularity(mc))
                          plot(mc, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = mc$membership
                        }
                        if (NetClusterType == "labelPropagation") {
                          system.time(lc <- label.propagation.community(g))
                          print(modularity(lc))
                          plot(lc, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = lc$membership
                        }
                      }
                      dev.off()
                      write.csv(netcsv, paste(matrix_dir, "/", 
                        chrName, "_network.csv", sep = ""), row.names = FALSE)
                    }
                    bedIplot = cbind(rbind(chrBedToBedInter[, 
                      14], chrBedToBedInter[, 15]), rbind(chrBedToBedInter[, 
                      15], chrBedToBedInter[, 14]))
                    hm_dim = dim(chrCmap)[1]
                    chrCmap = as.matrix(chrCmap)
                    hm_mean = mean(chrCmap)
                    for (j in 1:hm_dim) {
                      chrCmap[j, which(chrCmap[j, ] > 5 * hm_mean)] = 5 * 
                        hm_mean
                    }
                    chrhmCmap = melt(chrCmap)
                    print(paste("plot ", chrName, "bed picture", 
                      sep = ""))
                    if ((outputpdf == "TRUE") || (outputpdf == 
                      "true")) {
                      pdf(paste(matrix_dir, "/", chrName, "_bedplot.pdf", 
                        sep = ""), width = 8, height = 8)
                    }
                    else {
                      jpeg(paste(matrix_dir, "/", chrName, "_bedplot.jpeg", 
                        sep = ""), width = 1000, height = 1000, 
                        quality = 100)
                    }
                    grid.newpage()
                    heatmapViewport <- viewport(height = 0.5, 
                      width = 0.5, x = 0.25, y = 0.5)
                    scatterViewport <- viewport(height = 0.5, 
                      width = 0.5, x = 0.75, y = 0.5)
                    densityViewport <- viewport(height = 0.25, 
                      width = 0.5, x = 0.75, y = 0.125)
                    hmdensityViewport <- viewport(height = 0.25, 
                      width = 0.5, x = 0.25, y = 0.125)
                    jit = position_jitter(width = 0.5)
                    hmrange = range(chrhmCmap[, 1])
                    bedIplot = bedIplot + hmrange[1]
                    chrhm = ggplot(chrhmCmap, aes(x = Var1, y = Var2, 
                      fill = value)) + scale_y_discrete(breaks = seq(0, 
                      10, 5)) + xlab("chrom") + ylab("chrom") + 
                      scale_fill_gradient(low = "white", high = "red") + 
                      geom_tile() + guides(fill = FALSE)
                    chrbedplot = qplot(bedIplot[1, ], bedIplot[2, 
                      ], alpha = I(1/10), size = I(1)) + xlab("chrom") + 
                      ylab("chrom") + geom_jitter(position = jit, 
                      colour = "black", alpha = 1/100)
                    chrbeddensitydata = as.data.frame(c(chrBedToBedInter[, 
                      14], chrBedToBedInter[, 15]))
                    colnames(chrbeddensitydata) = "bed"
                    chrbeddensity = ggplot(chrbeddensitydata) + 
                      geom_density(aes(x = bed))
                    chrhmdensitydata = NULL
                    for (iiii in 1:chrTotSize) {
                      chrhmdensitydata = c(chrhmdensitydata, 
                        which(chrCmap[, iiii] > 0))
                    }
                    chrhmdensitydata = as.data.frame(chrhmdensitydata)
                    colnames(chrhmdensitydata) = "chrom"
                    chrcmapdensity = ggplot(chrhmdensitydata) + 
                      geom_density(aes(x = chrom))
                    print(chrhm, vp = heatmapViewport)
                    print(chrbedplot, vp = scatterViewport)
                    print(chrbeddensity, vp = densityViewport)
                    print(chrcmapdensity, vp = hmdensityViewport)
                    dev.off()
                    write.table(chrBedToBedInter, file = paste(matrix_dir, 
                      "/", chrName, "_BedToBedInter.txt", sep = ""), 
                      sep = "\t", row.names = FALSE, col.names = FALSE, 
                      quote = FALSE)
                    rm(chrCmap)
                    rm(chrBed)
                    rm(chrBedBin)
                    rm(chrBedMatrix)
                    rm(chrBedToBedInter)
                    gc()
                  }
                }
                else {
                  chrBed = choose_chr_bed(m_bed, chrName)
                  chrBed = chrBed[which(chrBed[, 3] < chrend), 
                    ]
                  chrBed[, 2] = chrBed[, 2] - chrstart
                  chrBed[, 3] = chrBed[, 3] - chrstart
                  chrBed = chrBed[which(chrBed[, 2] > 0), ]
                  chrBedNum = dim(chrBed)[1]
                  print(matrix_full_dir[i])
                  if (chrBedNum > 0) {
                    print(chrName)
                    chrCmap = read.table(file = matrix_full_dir[i], 
                      fill = TRUE, stringsAsFactors = FALSE)
                    tmpCmapStart = abs(ceiling((chrstart/(resolution * 
                      1000))))
                    tmpCmapEnd = abs(ceiling((chrend/(resolution * 
                      1000))))
                    chrTotSize = dim(chrCmap)[1]
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
                    chrBedBin = check_bed_bin(chrBed, resolution * 
                      1000)
                    chrBedMatrix = convert_bed_to_matrix(chrBed, 
                      resolution * 1000, chrName, chrTotSize, 
                      bedWindow)
                    chrBedToBedInter = find_bed_to_bed_interaction(chrCmap, 
                      chrBedMatrix, chrBedBin, chrBed, chrName, 
                      chrCmap)
                    if ((netplot == "true") || (netplot == "TRUE")) {
                      netgraph = data.frame(p1 = character(dim(chrBedToBedInter)[1]), 
                        p2 = character(dim(chrBedToBedInter)[1]), 
                        weight = numeric(dim(chrBedToBedInter)[1]))
                      netgraph[, 1] = as.data.frame(paste(chrBedToBedInter[, 
                        2], ":", (chrBedToBedInter[, 3] + chrstart), 
                        "-", (chrBedToBedInter[, 4] + chrstart), 
                        sep = ""))
                      netgraph[, 2] = as.data.frame(paste(chrBedToBedInter[, 
                        8], ":", (chrBedToBedInter[, 9] + chrstart), 
                        "-", (chrBedToBedInter[, 10] + chrstart), 
                        sep = ""))
                      netgraph[, 3] = as.data.frame(chrBedToBedInter[, 
                        13])
                      set.seed(1234)
                      if ((outputpdf == "TRUE") || (outputpdf == 
                        "true")) {
                        pdf(paste(matrix_dir, "/", chrName, "_netplot.pdf", 
                          sep = ""), width = 8, height = 8)
                      }
                      else {
                        jpeg(paste(matrix_dir, "/", chrName, 
                          "_netplot.jpeg", sep = ""), width = 1000, 
                          height = 1000, quality = 100)
                      }
                      g = graph.data.frame(netgraph, directed = F)
                      set.seed(1234)
                      netnodename = names(V(g))
                      netnodechrnum = regexpr(":", netnodename)
                      netnodechr = substr(netnodename, 1, netnodechrnum - 
                        1)
                      netnodestartnum = regexpr("-", netnodename)
                      netnodestart = as.numeric(substr(netnodename, 
                        netnodechrnum + 1, netnodestartnum - 
                          1))
                      netnodeendnum = nchar(netnodename)
                      netnodeend = as.numeric(substr(netnodename, 
                        netnodestartnum + 1, netnodeendnum))
                      netcsv = data.frame(chrom = character(0), 
                        start = numeric(0), end = numeric(0), 
                        degree = numeric(0), closeness = numeric(0), 
                        betweenness = numeric(0), Local_cluster_coefficient = numeric(0), 
                        Eigenvector_centrality = numeric(0), 
                        membership = numeric(0), stringsAsFactors = FALSE)
                      netdegree = degree(g)
                      netcloseness = closeness(g)
                      netbetweenness = betweenness(g)
                      netcoefficient = transitivity(g, type = "local")
                      netcentrality = evcent(g)$vector
                      netcsv[1:(length(netdegree)), 1] = netnodechr
                      netcsv[1:(length(netdegree)), 2] = netnodestart
                      netcsv[1:(length(netdegree)), 3] = netnodeend
                      netcsv[1:(length(netdegree)), 4] = as.data.frame(netdegree)
                      netcsv[1:(length(netdegree)), 5] = as.data.frame(netcloseness)
                      netcsv[1:(length(netdegree)), 6] = as.data.frame(netbetweenness)
                      netcsv[1:(length(netdegree)), 7] = as.data.frame(netcoefficient)
                      netcsv[1:(length(netdegree)), 8] = as.data.frame(netcentrality)
                      colors <- c("#fff7f3", "#fff7f3", "#fde0dd", 
                        "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", 
                        "#b0017e", "#7a0177", "#49006a")
                      weight_range = range(E(g)$weight)
                      E(g)$color = colors[1]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1])/10)]$color = colors[2]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 2/10)]$color = colors[3]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 3/10)]$color = colors[4]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 4/10)]$color = colors[5]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 5/10)]$color = colors[6]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 6/10)]$color = colors[7]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 7/10)]$color = colors[8]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 8/10)]$color = colors[9]
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 9/10)]$color = colors[10]
                      edge_width = 0.05
                      E(g)$width = edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1])/10)]$width = 2 * edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 2/10)]$width = 3 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 3/10)]$width = 4 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 4/10)]$width = 5 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 5/10)]$width = 6 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 6/10)]$width = 7 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 7/10)]$width = 8 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 8/10)]$width = 9 * 
                        edge_width
                      E(g)[weight >= (weight_range[1] + (weight_range[2] - 
                        weight_range[1]) * 9/10)]$width = 10 * 
                        edge_width
                      if (NetVertexChangeSize == "degree") {
                        V(g)$deg <- netcsv[, 4]
                        deg_range = range(netcsv[, 4])
                        V(g)$size = NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1])/5)]$size = 2 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                      }
                      else if (NetVertexChangeSize == "closeness") {
                        V(g)$deg <- netcsv[, 5]
                        deg_range = range(netcsv[, 5])
                        V(g)$size = NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1])/5)]$size = 2 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                      }
                      else if (NetVertexChangeSize == "betweenness") {
                        V(g)$deg <- netcsv[, 6]
                        deg_range = range(netcsv[, 6])
                        V(g)$size = NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1])/5)]$size = 2 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                      }
                      else if (NetVertexChangeSize == "Local_cluster_coefficient") {
                        V(g)$deg <- netcsv[, 7]
                        deg_range = range(netcsv[, 7])
                        V(g)$size = NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1])/5)]$size = 2 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                      }
                      else if (NetVertexChangeSize == "Eigenvector_centrality") {
                        V(g)$deg <- netcsv[, 8]
                        deg_range = range(netcsv[, 8])
                        V(g)$size = NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1])/5)]$size = 2 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 2/5)]$size = 3 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 3/5)]$size = 4 * NetVertexSize
                        V(g)[deg >= (deg_range[1] + (deg_range[2] - 
                          deg_range[1]) * 4/5)]$size = 5 * NetVertexSize
                      }
                      else {
                        V(g)$size = NetVertexSize
                      }
                      if (net_layout == "layout.fruchterman.reingold") {
                        if (NetClusterType == "NULL") {
                          plot(g, layout = layout.fruchterman.reingold, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = "NULL"
                        }
                        if (NetClusterType == "edgeBetweenness") {
                          system.time(ec <- edge.betweenness.community(g))
                          print(modularity(ec))
                          netcsv[1:(length(netdegree)), 9] = ec$membership
                          plot(ec, g, layout = layout.fruchterman.reingold, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "walktrap") {
                          system.time(wc <- walktrap.community(g))
                          netcsv[1:(length(netdegree)), 9] = wc$membership
                          print(modularity(wc))
                          plot(wc, g, layout = layout.fruchterman.reingold, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "multileve") {
                          system.time(mc <- multilevel.community(g, 
                            weights = NA))
                          print(modularity(mc))
                          plot(mc, g, layout = layout.fruchterman.reingold, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = mc$membership
                        }
                        if (NetClusterType == "labelPropagation") {
                          system.time(lc <- label.propagation.community(g))
                          print(modularity(lc))
                          plot(lc, g, layout = layout.fruchterman.reingold, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = lc$membership
                        }
                      }
                      else if (net_layout == "layout.circle") {
                        if (NetClusterType == "NULL") {
                          plot(g, layout = layout.circle, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = "NULL"
                        }
                        if (NetClusterType == "edgeBetweenness") {
                          system.time(ec <- edge.betweenness.community(g))
                          print(modularity(ec))
                          netcsv[1:(length(netdegree)), 9] = ec$membership
                          plot(ec, g, layout = layout.circle, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "walktrap") {
                          system.time(wc <- walktrap.community(g))
                          netcsv[1:(length(netdegree)), 9] = wc$membership
                          print(modularity(wc))
                          plot(wc, g, layout = layout.circle, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "multileve") {
                          system.time(mc <- multilevel.community(g, 
                            weights = NA))
                          print(modularity(mc))
                          plot(mc, g, layout = layout.circle, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = mc$membership
                        }
                        if (NetClusterType == "labelPropagation") {
                          system.time(lc <- label.propagation.community(g))
                          print(modularity(lc))
                          plot(lc, g, layout = layout.circle, 
                            vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = lc$membership
                        }
                      }
                      else {
                        if (NetClusterType == "NULL") {
                          plot(g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = "NULL"
                        }
                        if (NetClusterType == "edgeBetweenness") {
                          system.time(ec <- edge.betweenness.community(g))
                          print(modularity(ec))
                          netcsv[1:(length(netdegree)), 9] = ec$membership
                          plot(ec, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "walktrap") {
                          system.time(wc <- walktrap.community(g))
                          netcsv[1:(length(netdegree)), 9] = wc$membership
                          print(modularity(wc))
                          plot(wc, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                        }
                        if (NetClusterType == "multileve") {
                          system.time(mc <- multilevel.community(g, 
                            weights = NA))
                          print(modularity(mc))
                          plot(mc, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = mc$membership
                        }
                        if (NetClusterType == "labelPropagation") {
                          system.time(lc <- label.propagation.community(g))
                          print(modularity(lc))
                          plot(lc, g, layout = layout.auto, vertex.label.dist = NetVertexLableDist, 
                            vertex.color = NetVertexColor, edge.arrow.size = 0.05, 
                            vertex.label.cex = NetVertexLabelCex)
                          netcsv[1:(length(netdegree)), 9] = lc$membership
                        }
                      }
                      dev.off()
                      write.csv(netcsv, paste(matrix_dir, "/", 
                        chrName, "_network.csv", sep = ""), row.names = FALSE)
                    }
                    bedIplot = cbind(rbind(chrBedToBedInter[, 
                      14], chrBedToBedInter[, 15]), rbind(chrBedToBedInter[, 
                      15], chrBedToBedInter[, 14]))
                    hm_dim = dim(chrCmap)[1]
                    chrCmap = as.matrix(chrCmap)
                    hm_mean = mean(chrCmap)
                    for (j in 1:hm_dim) {
                      chrCmap[j, which(chrCmap[j, ] > 5 * hm_mean)] = 5 * 
                        hm_mean
                    }
                    chrhmCmap = melt(chrCmap)
                    print(paste("plot ", chrName, "bed picture", 
                      sep = ""))
                    if ((outputpdf == "TRUE") || (outputpdf == 
                      "true")) {
                      pdf(paste(matrix_dir, "/", chrName, "_bedplot.pdf", 
                        sep = ""), width = 8, height = 8)
                    }
                    else {
                      jpeg(paste(matrix_dir, "/", chrName, "_bedplot.jpeg", 
                        sep = ""), width = 1000, height = 1000, 
                        quality = 100)
                    }
                    grid.newpage()
                    heatmapViewport <- viewport(height = 0.5, 
                      width = 0.5, x = 0.25, y = 0.5)
                    scatterViewport <- viewport(height = 0.5, 
                      width = 0.5, x = 0.75, y = 0.5)
                    densityViewport <- viewport(height = 0.25, 
                      width = 0.5, x = 0.75, y = 0.125)
                    hmdensityViewport <- viewport(height = 0.25, 
                      width = 0.5, x = 0.25, y = 0.125)
                    jit = position_jitter(width = 0.5)
                    hmrange = range(chrhmCmap[, 1])
                    bedIplot = bedIplot + hmrange[1]
                    chrhm = ggplot(chrhmCmap, aes(x = Var1, y = Var2, 
                      fill = value)) + scale_y_discrete(breaks = seq(0, 
                      10, 5)) + xlab("chrom") + ylab("chrom") + 
                      scale_fill_gradient(low = "white", high = "red") + 
                      geom_tile() + guides(fill = FALSE)
                    chrbedplot = qplot(bedIplot[1, ], bedIplot[2, 
                      ], alpha = I(1/10), size = I(1)) + xlab("chrom") + 
                      ylab("chrom") + geom_jitter(position = jit, 
                      colour = "black", alpha = 1/100)
                    chrbeddensitydata = as.data.frame(c(chrBedToBedInter[, 
                      14], chrBedToBedInter[, 15]))
                    colnames(chrbeddensitydata) = "bed"
                    chrbeddensity = ggplot(chrbeddensitydata) + 
                      geom_density(aes(x = bed))
                    chrhmdensitydata = NULL
                    for (iiii in 1:chrTotSize) {
                      chrhmdensitydata = c(chrhmdensitydata, 
                        which(chrCmap[, iiii] > 0))
                    }
                    chrhmdensitydata = as.data.frame(chrhmdensitydata)
                    colnames(chrhmdensitydata) = "chrom"
                    chrcmapdensity = ggplot(chrhmdensitydata) + 
                      geom_density(aes(x = chrom))
                    print(chrhm, vp = heatmapViewport)
                    print(chrbedplot, vp = scatterViewport)
                    print(chrbeddensity, vp = densityViewport)
                    print(chrcmapdensity, vp = hmdensityViewport)
                    dev.off()
                    write.table(chrBedToBedInter, file = paste(matrix_dir, 
                      "/", chrName, "_BedToBedInter.txt", sep = ""), 
                      sep = "\t", row.names = FALSE, col.names = FALSE, 
                      quote = FALSE)
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
