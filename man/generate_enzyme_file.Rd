\name{generate_enzyme_file}
\alias{generate_enzyme_file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generate enzyme file
}
\description{
generate enzyme file
}
\usage{
generate_enzyme_file(enzyme = "HindIII", enzymesite = "AAGCTT", chrom_file = "chrom_hg19.sizes", enzymedir = "annotation", enzymeoverhangs5 = 1, genomeName = "hg19", resolution = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enzyme}{
%%     ~~Describe \code{enzyme} here~~
}
  \item{enzymesite}{
%%     ~~Describe \code{enzymesite} here~~
}
  \item{chrom_file}{
%%     ~~Describe \code{chrom_file} here~~
}
  \item{enzymedir}{
%%     ~~Describe \code{enzymedir} here~~
}
  \item{enzymeoverhangs5}{
%%     ~~Describe \code{enzymeoverhangs5} here~~
}
  \item{genomeName}{
%%     ~~Describe \code{genomeName} here~~
}
  \item{resolution}{
%%     ~~Describe \code{resolution} here~~
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
function (enzyme = "HindIII", enzymesite = "AAGCTT", chrom_file = "chrom_hg19.sizes", 
    enzymedir = "annotation", enzymeoverhangs5 = 1, genomeName = "hg19", 
    resolution = 100) 
{
    chrom_info = read.table(chrom_file, fill = TRUE, stringsAsFactors = FALSE)
    all_genomedb = read.csv("all_genome_db.csv", header = FALSE)
    enzymesitesfile = paste(enzyme, "_resfrag_", genomeName, 
        ".bed", sep = "")
    enzyme_dir = list.files(path = enzymedir, full.names = F, 
        pattern = enzymesitesfile)
    if (length(enzyme_dir) == 0) {
        for (i in 1:length(all_genomedb)) {
            dbfind = -1
            dbrequire = FALSE
            r_genomedb = NULL
            dbfind = regexpr(genomeName, all_genomedb[1, i])
            if (dbfind != -1) {
                r_genomedb = as.character(all_genomedb[1, i])
                dbrequire = require(all_genomedb[1, i], character.only = TRUE)
                break
            }
        }
        if (!is.null(r_genomedb)) {
            if (dbrequire == FALSE) {
                print(paste("try to install the R package ", 
                  r_genomedb, sep = ""))
                source("http://www.bioconductor.org/biocLite.R")
                biocLite(r_genomedb)
            }
            else {
                all_chr <- chrom_info[, 1]
                resFrag <- getRestrictionFragmentsPerChromosome(resSite = enzymesite, 
                  chromosomes = all_chr, overhangs5 = enzymeoverhangs5, 
                  genomePack = r_genomedb)
                allRF <- do.call("c", resFrag)
                names(allRF) <- unlist(sapply(resFrag, function(x) {
                  paste0("HIC_", seqlevels(x), "_", 1:length(x))
                }))
                export(allRF, format = "bed", con = paste(enzymedir, 
                  "/", enzymesitesfile, sep = ""))
            }
        }
        else {
            print(paste("can not find the genome file", sep = ""))
        }
    }
    else {
        print(paste("enzymesites file ", enzymesitesfile, " is already existed, skip to generate it"))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
