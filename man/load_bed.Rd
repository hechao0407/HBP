\name{load_bed}
\alias{load_bed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
the function to load bed file
}
\description{
the function to load bed file
}
\usage{
load_bed(bedfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bedfile}{
%%     ~~Describe \code{bedfile} here~~
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
function (bedfile) 
{
    uniques <- read.table(file = bedfile, fill = TRUE, stringsAsFactors = FALSE, 
        skip = 1)
    uniques <- uniques[, c(1, 2, 3, 4)]
    colnames(uniques) <- c("chr", "start", "end", "name")
    uniques <- na.omit(uniques)
    ordereduniques <- uniques[with(uniques, order(chr, as.numeric(start))), 
        ]
    if (dim(ordereduniques)[1] == 0) {
        stop("No data loaded! \n\n")
    }
    else {
        return(ordereduniques)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
