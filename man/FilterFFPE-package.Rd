\name{FilterFFPE-package}
\alias{FilterFFPE-package}
\alias{FilterFFPE}
\docType{package}
\title{
\packageTitle{FilterFFPE}
}
\description{
\packageDescription{FilterFFPE}
}
\details{

\packageDESCRIPTION{FilterFFPE}
The next-generation sequencing (NGS) reads from formalin-fixed paraffin-embedded
(FFPE) samples contain numerous artifact chimeric reads, which can lead to a
large number of false positive structural variation (SV) calls. This package
finds and filters these artifact chimeric reads from BAM files of FFPE samples.

\packageIndices{FilterFFPE}
There are three available functions to find and/or filter artifact chimeric
reads of FFPE samples:

1. \code{\link{findArtifactChimericReads}}: Find artifact chimeric reads in BAM
file of FFPE sample.

2. \code{\link{filterBamByReadNames}}: Filter reads from BAM file by read names.

3. \code{\link{FFPEReadFilter}}: Find and filter artifact chimeric reads in BAM
file of FFPE sample.
}
\author{
Lanying Wei <lanying.wei@uni-muenster.de>
}

\keyword{package}
\seealso{
\code{\link{FilterFFPE}}, \code{\link{filterBamByReadNames}},
\code{\link{FFPEReadFilter}}
}
\examples{
file <- system.file("extdata", "example.bam", package = "FilterFFPE")
outFolder <- tempdir()
FFPEReadsFile <- paste0(outFolder, "/example.FFPEReads.txt")
dupChimFile <- paste0(outFolder, "/example.dupChim.txt")
destination <- paste0(outFolder, "/example.FilterFFPE.bam")
FFPEReadFilter(file = file, threads = 2, destination = destination,
               overwrite = TRUE, FFPEReadsFile = FFPEReadsFile,
               dupChimFile = dupChimFile)
}
