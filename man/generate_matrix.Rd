\name{generate_matrix}
\alias{generate_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generate the matrix file
}
\description{
generate the matrix file
}
\usage{
generate_matrix(all_hic_file, all_bed_file, outputpdf = "FALSE", matrix_dir = "dm3", resolution = 5, chrom_file = "chrom_dm3.sizes")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all_hic_file}{
%%     ~~Describe \code{all_hic_file} here~~
}
  \item{all_bed_file}{
%%     ~~Describe \code{all_bed_file} here~~
}
  \item{outputpdf}{
%%     ~~Describe \code{outputpdf} here~~
}
  \item{matrix_dir}{
%%     ~~Describe \code{matrix_dir} here~~
}
  \item{resolution}{
%%     ~~Describe \code{resolution} here~~
}
  \item{chrom_file}{
%%     ~~Describe \code{chrom_file} here~~
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
function (all_hic_file, all_bed_file, outputpdf = "FALSE", matrix_dir = "dm3", 
    resolution = 5, chrom_file = "chrom_dm3.sizes") 
{
    all_bed_data = read.table(file = all_bed_file, fill = TRUE, 
        stringsAsFactors = FALSE)
    all_hic_data = read.table(file = all_hic_file, fill = TRUE, 
        stringsAsFactors = FALSE)
    chrom_info = read.table(chrom_file, fill = TRUE, stringsAsFactors = FALSE)
    chr_num = dim(chrom_info)[1]
    inter_interaction = NULL
    for (tttt in 1:chr_num) {
        chrname = chrom_info[tttt, 1]
        chr_bed = NULL
        chr_bed = rbind(chr_bed, all_bed_data[all_bed_data[, 
            1] == chrname, ])
        chr_max_length = dim(chr_bed)[1]
        chr_hic = NULL
        chr_hic = all_hic_data[all_hic_data[, 1] >= chr_bed[1, 
            4], ]
        chr_hic = chr_hic[chr_hic[, 1] <= chr_bed[chr_max_length, 
            4], ]
        inter_interaction = rbind(inter_interaction, chr_hic[chr_hic[, 
            2] < chr_bed[1, 4], ])
        inter_interaction = rbind(inter_interaction, chr_hic[chr_hic[, 
            2] > chr_bed[chr_max_length, 4], ])
        chr_hic = chr_hic[chr_hic[, 2] >= chr_bed[1, 4], ]
        chr_hic = chr_hic[chr_hic[, 2] <= chr_bed[chr_max_length, 
            4], ]
        chr_hic_data = matrix(data = 0, nrow = chr_max_length, 
            ncol = chr_max_length)
        for (i in 1:dim(chr_hic)[1]) {
            chr_hic_data[(chr_hic[i, 1] - chr_bed[1, 4]), (chr_hic[i, 
                2] - chr_bed[1, 4])] = chr_hic[i, 3]
            chr_hic_data[(chr_hic[i, 2] - chr_bed[1, 4]), (chr_hic[i, 
                1] - chr_bed[1, 4])] = chr_hic[i, 3]
        }
        tmpfilename = paste(matrix_dir, "/", chrname, ".matrix", 
            sep = "")
        write.table(chr_hic_data, tmpfilename, quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
        if ((outputpdf == "TRUE") || (outputpdf == "true")) {
            pdf(paste(matrix_dir, "/", chrname, "_heatmap.pdf", 
                sep = ""))
        }
        else {
            jpeg(paste(matrix_dir, "/", chrname, "_heatmap.jpeg", 
                sep = ""), width = 1000, height = 1000, quality = 100)
        }
        hm_dim = dim(chr_hic_data)[1]
        hm_mean = mean(chr_hic_data)
        for (n in 1:hm_dim) {
            for (nn in 1:hm_dim) {
                if (chr_hic_data[n, nn] > 10 * hm_mean) 
                  chr_hic_data[n, nn] = 10 * hm_mean
            }
        }
        m_color = heat.colors(100)
        n_color = m_color
        for (i in 1:100) {
            m_color[i] = n_color[101 - i]
        }
        tmphm = levelplot(chr_hic_data, col.regions = m_color, 
            xlim = c(0, hm_dim), ylim = c(0, hm_dim), xlab = "chrom", 
            ylab = "chrom")
        plot(tmphm)
        dev.off()
        rm(chr_hic_data)
        gc()
        print(paste(chrname, " finish!", sep = ""))
    }
    tmpfilename = paste(matrix_dir, "/inter_chrom.matrix", sep = "")
    write.table(inter_interaction, tmpfilename, quote = FALSE, 
        sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
