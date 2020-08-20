\name{findArtifactChimericReads}
\alias{findArtifactChimericReads}

\title{
Find artifact chimeric reads in BAM file of FFPE sample
}
\description{
Artifact chimeric reads are enriched in NGS data of FFPE samples, these reads
can lead to a large number of false positive SV calls. This function finds
these artifact chimeric reads.
}
\usage{
findArtifactChimericReads(file, maxReadsOfSameBreak=2, minMapBase=1,
threads=1, FFPEReadsFile=sub("\\\.bam(\\\.gz)?", ".FFPEReads.txt", file),
dupChimFile=sub("\\\.bam(\\\.gz)?", ".dupChim.txt", file))
}
\arguments{
  \item{file}{
  Path to the BAM file.
}
  \item{maxReadsOfSameBreak}{
  The maximum allowed number of artifact chimeric reads sharing a false positive
  breakpoint. If the number of reads sharing the same breakpoint exceeds this
  number, these reads are not recognized as artifact chimeric reads. Reads
  marked as PCR or optical duplicates are excluded from the calculation.
}
  \item{minMapBase}{
  The minimum required length (bp) of a short complementary mapping for an
  artifact chimeric read. Artifact chimeric reads are derived from the
  combination of two single-stranded DNA fragments linked by short reverse
  complementary mapping. Reads with such complementary sequences shorter than
  this length are not recognized as artifact chimeric reads.
}
  \item{threads}{
  Number of threads to use. Multi-threading can speed up the process.
}
  \item{FFPEReadsFile}{
  Path of the output txt file with artifact chimeric read names.
}
  \item{dupChimFile}{
  Path of the output txt file with read names of PCR or optical duplicates of
  all chimeric reads.
}
}
\details{
The next-generation sequencing (NGS) reads from formalin-fixed paraffin-embedded
(FFPE) samples contain numerous artifact chimeric reads, which can lead to a
large number of false positive structural variation (SV) calls. This function
finds the read names of these artifact chimeric reads. To further filter these
reads, \code{\link{filterBamByReadNames}} can be applied.
}
\value{
A character vector of artifact chimeric read names.
}
\author{
Lanying Wei <lanying.wei@uni-muenster.de>
}
\seealso{
\code{\link{FilterFFPE}}, \code{\link{filterBamByReadNames}},
\code{\link{FFPEReadFilter}}
}
\examples{
file <- system.file("extdata", "example.bam", package = "FilterFFPE")
outFolder <- tempdir()
FFPEReadsFile <- paste0(outFolder, "/example.FFPEReads.txt")
dupChimFile <- paste0(outFolder, "/example.dupChim.txt")
artifactReads <- findArtifactChimericReads(file = file, threads = 2,
                                           FFPEReadsFile = FFPEReadsFile,
                                           dupChimFile = dupChimFile)
head(artifactReads)
}

