\name{run_hicpro}
\alias{run_hicpro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
help users to run the hic-pro
}
\description{
help users to run the hic-pro
}
\usage{
run_hicpro(hicpro_path = "HiC-Pro", inputfile = "rawdata", configfile = "config-hicpro", outdir = "hg19")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hicpro_path}{
%%     ~~Describe \code{hicpro_path} here~~
}
  \item{inputfile}{
%%     ~~Describe \code{inputfile} here~~
}
  \item{configfile}{
%%     ~~Describe \code{configfile} here~~
}
  \item{outdir}{
%%     ~~Describe \code{outdir} here~~
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
function (hicpro_path = "HiC-Pro", inputfile = "rawdata", configfile = "config-hicpro", 
    outdir = "hg19") 
{
    hiccmd = paste(hicpro_path, " -i ", inputfile, " -o ", outdir, 
        " -c ", configfile, sep = "")
    print(hiccmd)
    system(hiccmd)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
