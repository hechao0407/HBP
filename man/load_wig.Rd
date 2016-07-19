\name{load_wig}
\alias{load_wig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
the function to load tracks
}
\description{
the function to load tracks
}
\usage{
load_wig(wigfile, resolution = 2000, chrom = "chr2L", chrTotSize, chrstart, chrend)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wigfile}{
%%     ~~Describe \code{wigfile} here~~
}
  \item{resolution}{
%%     ~~Describe \code{resolution} here~~
}
  \item{chrom}{
%%     ~~Describe \code{chrom} here~~
}
  \item{chrTotSize}{
%%     ~~Describe \code{chrTotSize} here~~
}
  \item{chrstart}{
%%     ~~Describe \code{chrstart} here~~
}
  \item{chrend}{
%%     ~~Describe \code{chrend} here~~
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
function (wigfile, resolution = 2000, chrom = "chr2L", chrTotSize, 
    chrstart, chrend) 
{
    tmp_wig <- read.table(file = wigfile, fill = TRUE, stringsAsFactors = FALSE, 
        skip = 1)
    ooo <- data.frame(seg.name = character(chrTotSize), seg.po = numeric(chrTotSize), 
        value = numeric(chrTotSize), stringsAsFactors = FALSE)
    ooo_num <- data.frame(wig_num = numeric(chrTotSize))
    sec = NULL
    sec = rbind(sec, tmp_wig[tmp_wig[, 1] == chrom, ])
    sec[, 2] = as.numeric(sec[, 2])
    sec[, 3] = as.numeric(sec[, 3])
    sec[, 4] = as.numeric(sec[, 4])
    if (chrend > 0) {
        sec = sec[which(sec[, 3] < chrend), ]
        sec[, 2] = sec[, 2] - chrstart
        sec[, 3] = sec[, 3] - chrstart
        sec = sec[which(sec[, 2] > 0), ]
    }
    secdim = dim(sec)[1]
    if (is.null(secdim)) {
        return(NULL)
    }
    else if (secdim < 2) {
        return(NULL)
    }
    else {
        ooo[, 1] = chrom
        ooo[, 2:3] = 0
        ooo_num[, 1] = 0
        wig_bin = 1
        sec_num = dim(sec)[1]
        for (i in 1:sec_num) {
            wig_bin = ceiling((sec[i, 2] + sec[i, 3])/(2 * resolution))
            if (wig_bin > chrTotSize) {
                wig_bin = chrTotSize
            }
            ooo[wig_bin, 3] = ooo[wig_bin, 3] + sec[i, 4]
            ooo_num[wig_bin, 1] = ooo_num[wig_bin, 1] + 1
            ooo[wig_bin, 2] = wig_bin
        }
        for (i in 1:chrTotSize) {
            if (ooo[i, 2] == 0) {
                ooo[i, 2] = i
            }
            if (ooo_num[i, 1] > 0) {
                ooo[i, 3] = (ooo[i, 3])/(ooo_num[i, 1])
            }
            else {
                ooo[i, 3] = 0
            }
        }
        return(ooo)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
