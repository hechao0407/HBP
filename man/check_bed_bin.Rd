\name{check_bed_bin}
\alias{check_bed_bin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
process the bed file data to bin
}
\description{
process the bed file data to bin
}
\usage{
check_bed_bin(m_bed, bin = 2000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m_bed}{
%%     ~~Describe \code{m_bed} here~~
}
  \item{bin}{
%%     ~~Describe \code{bin} here~~
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
function (m_bed, bin = 2000) 
{
    ooo = matrix(data = 0, nrow = dim(m_bed)[1], ncol = 4)
    for (i in 1:dim(m_bed)[1]) {
        st = (ceiling(m_bed[i, 2]/bin))
        ed = (ceiling(m_bed[i, 3]/bin))
        st1 = (ceiling((m_bed[i, 2] - 500)/bin))
        ed1 = (ceiling((m_bed[i, 3] + 500)/bin))
        if (st != ed) {
            ooo[i, 2] = st
            ooo[i, 3] = ed
            if (st1 != st) {
                ooo[i, 1] = st1
            }
            if (ed1 != ed) {
                ooo[i, 4] = ed1
            }
        }
        else {
            ooo[i, 2] = st
            if (st1 != st) {
                ooo[i, 1] = st1
            }
            if (ed1 != ed) {
                ooo[i, 4] = ed1
            }
        }
    }
    return(ooo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
