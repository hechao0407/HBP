\name{calculate_omiccircos_data}
\alias{calculate_omiccircos_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
output the interaction that our bed interact with ourbed.
}
\description{
output the interaction that our bed interact with ourbed.
}
\usage{
calculate_omiccircos_data(cmap, bed_matrix, bed_bin, bedfile, chr = "chr4", st_cmap, m_threshold = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cmap}{
%%     ~~Describe \code{cmap} here~~
}
  \item{bed_matrix}{
%%     ~~Describe \code{bed_matrix} here~~
}
  \item{bed_bin}{
%%     ~~Describe \code{bed_bin} here~~
}
  \item{bedfile}{
%%     ~~Describe \code{bedfile} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{st_cmap}{
%%     ~~Describe \code{st_cmap} here~~
}
  \item{m_threshold}{
%%     ~~Describe \code{m_threshold} here~~
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
function (cmap, bed_matrix, bed_bin, bedfile, chr = "chr4", st_cmap, 
    m_threshold = 0) 
{
    n = dim(cmap)[1]
    bed_count = 0
    nm_count = 0
    all_count = 0
    cccc = 0
    ooo <- data.frame(chr_nm1 = character(0), n_start1 = numeric(0), 
        bed_seq1 = character(0), chr_nm2 = character(0), n_start2 = numeric(0), 
        bed_seq2 = character(0), read_count = numeric(0), stringsAsFactors = FALSE)
    ooo = lapply(rbind(1:(dim(bed_matrix)[1] - 1)), calculate_omiccircos_data_solo, 
        cmap = cmap, bed_matrix = bed_matrix, bed_bin = bed_bin, 
        bedfile = bedfile, chr = chr, st_cmap = st_cmap, m_threshold = m_threshold)
    final_ooo <- data.frame(chr_nm1 = character(0), n_start1 = numeric(0), 
        bed_seq1 = character(0), chr_nm2 = character(0), n_start2 = numeric(0), 
        bed_seq2 = character(0), read_count = numeric(0), stringsAsFactors = FALSE)
    for (tt in 1:(dim(bed_matrix)[1] - 1)) {
        if (!is.null(ooo[tt][[1]])) {
            final_ooo = rbind(final_ooo, as.data.frame(ooo[tt][[1]]))
        }
    }
    colnames(final_ooo) = c("chr_nm1", "n_start1", "bed_seq1", 
        "chr_nm2", "n_start2", "bed_seq2", "read_count")
    return(final_ooo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
