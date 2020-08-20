\name{filterBamByReadNames}
\alias{filterBamByReadNames}

\title{
Filter reads from BAM file by read names
}
\description{
Generate filtered BAM file that does not contain reads with the input read
names.
}
\usage{
filterBamByReadNames(file, readsToFilter, index=file,
destination=sub("\\\.bam(\\\.gz)?", ".FilterFFPE.bam", file), overwrite=FALSE)
}
\arguments{
  \item{file}{
  Path to the input BAM file.
}
  \item{readsToFilter}{
  A character vector of read names to filter.
}
  \item{index}{
  Path of the index file of the input BAM file.
}
  \item{destination}{
  Path of the output filtered BAM file.
}
  \item{overwrite}{
  Boolean value indicating whether the destination can be over-written if it
  already exists.
}
}
\details{
Generate filtered BAM file that does not contain reads with the input read
names, index file is also created.
}
\value{
The file name of the created destination file.
}
\author{
Lanying Wei <lanying.wei@uni-muenster.de>
}
\seealso{
\code{\link{FilterFFPE}}, \code{\link{findArtifactChimericReads}},
\code{\link{FFPEReadFilter}}
}
\examples{
file <- system.file("extdata", "example.bam", package = "FilterFFPE")
outFolder <- tempdir()
FFPEReadsFile <- paste0(outFolder, "/example.FFPEReads.txt")
dupChimFile <- paste0(outFolder, "/example.dupChim.txt")
destination <- paste0(outFolder, "/example.FilterFFPE.bam")
artifactReads <- findArtifactChimericReads(file = file, threads = 2,
                                           FFPEReadsFile = FFPEReadsFile,
                                           dupChimFile = dupChimFile)
dupChim <- readLines(dupChimFile)
readsToFilter <- c(artifactReads, dupChim)
filterBamByReadNames(file = file, readsToFilter = readsToFilter,
                     destination = destination, overwrite=TRUE)
}
