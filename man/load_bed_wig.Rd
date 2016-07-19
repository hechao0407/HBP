\name{load_bed_wig}
\alias{load_bed_wig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
the function to load wig file in the specific sites
}
\description{
the function to load wig file in the specific sites
}
\usage{
load_bed_wig(wigfile, chrBed, chrom, chrstart, chrend, m_win)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wigfile}{
%%     ~~Describe \code{wigfile} here~~
}
  \item{chrBed}{
%%     ~~Describe \code{chrBed} here~~
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
  \item{m_win}{
%%     ~~Describe \code{m_win} here~~
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
function (wigfile, chrBed, chrom, chrstart, chrend, m_win) 
{
    all_wig = load_all_wig(wigfile)
    all_wig = all_wig[all_wig$chr == chrom, ]
    all_wig[, 2] = as.numeric(all_wig[, 2])
    all_wig[, 3] = as.numeric(all_wig[, 3])
    m_result = data.frame(chrom = character(0), start = numeric(0), 
        end = numeric(0), wig_value = numeric(0), stringsAsFactors = FALSE)
    dim_bed = dim(chrBed)[1]
    m_result[1:dim_bed, 1:3] = chrBed[, 1:3]
    if (chrend > 0) {
        all_wig = all_wig[which(all_wig[, 3] < chrend), ]
    }
    if (chrstart > 0) {
        all_wig[, 2] = all_wig[, 2] - chrstart
        all_wig[, 3] = all_wig[, 3] - chrstart
        all_wig = all_wig[which(all_wig[, 2] >= 0), ]
    }
    for (i in 1:dim_bed) {
        tmp_bed_wig = which(all_wig[, 2] >= chrBed[i, 2] - m_win)
        tmp_bed_wig = tmp_bed_wig[which(all_wig[tmp_bed_wig, 
            3] <= chrBed[i, 3] + m_win)]
        m_result[i, 4] = mean(as.numeric(all_wig[tmp_bed_wig, 
            4]))
    }
    return(m_result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
