\name{convert_bed_to_matrix}
\alias{convert_bed_to_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convert the bed file to matrix file for the next analysis
}
\description{
convert the bed file to matrix file for the next analysis
}
\usage{
convert_bed_to_matrix(bedfile, bin = 2000, chrom_list = "chr2L", tot_size = 0, bed_window = 2000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bedfile}{
%%     ~~Describe \code{bedfile} here~~
}
  \item{bin}{
%%     ~~Describe \code{bin} here~~
}
  \item{chrom_list}{
%%     ~~Describe \code{chrom_list} here~~
}
  \item{tot_size}{
%%     ~~Describe \code{tot_size} here~~
}
  \item{bed_window}{
%%     ~~Describe \code{bed_window} here~~
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
function (bedfile, bin = 2000, chrom_list = "chr2L", tot_size = 0, 
    bed_window = 2000) 
{
    if (tot_size == 0) {
        print("please input tot_size")
    }
    else {
        cmap = matrix(data = 0, nrow = tot_size, ncol = 1)
        sec = NULL
        for (j in 1:length(chrom_list)) {
            sec = rbind(sec, bedfile[bedfile$chr == chrom_list[j], 
                ])
        }
        ordereduniques = sec
        count = 0
        for (i in 1:dim(ordereduniques)[1]) {
            st = (ceiling(ordereduniques[i, 2]/bin))
            ed = (ceiling(ordereduniques[i, 3]/bin))
            st1 = (ceiling((ordereduniques[i, 2] - bed_window)/bin))
            ed1 = (ceiling((ordereduniques[i, 3] + bed_window)/bin))
            if (st < 1) {
                st = 1
            }
            if (st1 < 1) {
                st1 = 1
            }
            if (ed < 1) {
                ed = 1
            }
            if (ed1 < 1) {
                ed1 = 1
            }
            if (st > tot_size) {
                st = tot_size
            }
            if (st1 > tot_size) {
                st1 = tot_size
            }
            if (ed > tot_size) {
                ed = tot_size
            }
            if (ed1 > tot_size) {
                ed1 = tot_size
            }
            if (st != ed) {
                cmap[st, 1] = cmap[st, 1] + 1
                count = count + 1
                cmap[ed, 1] = cmap[ed, 1] + 1
                count = count + 1
                if (st1 != st) {
                  cmap[st1, 1] = cmap[st1, 1] + 1
                  count = count + 1
                }
                if (ed1 != ed) {
                  cmap[ed1, 1] = cmap[ed1, 1] + 1
                  count = count + 1
                }
            }
            else {
                cmap[ed, 1] = cmap[ed, 1] + 1
                count = count + 1
                if (st1 != st) {
                  cmap[st1, 1] = cmap[st1, 1] + 1
                  count = count + 1
                }
                if (ed1 != ed) {
                  cmap[ed1, 1] = cmap[ed1, 1] + 1
                  count = count + 1
                }
            }
        }
    }
    cat(sprintf("Mapped Fragments: \%s\n", count))
    return(cmap)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
