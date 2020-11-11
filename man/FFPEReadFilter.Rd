\name{FFPEReadFilter}
\alias{FFPEReadFilter}
\title{
Find and filter artifact chimeric reads in BAM file of FFPE sample
}
\description{
Artifact chimeric reads are enriched in NGS data of FFPE samples, these reads
can lead to a large number of false positive SV calls. This function finds
and filters these artifact chimeric reads.
}
\usage{
FFPEReadFilter(file, maxReadsOfSameBreak=2, minMapBase=1, threads=1,
index=file, destination=sub("\\\.bam(\\\.gz)?", ".FilterFFPE.bam", file),
overwrite=FALSE, FFPEReadsFile=sub("\\\.bam(\\\.gz)?", ".FFPEReads.txt", file),
dupChimFile=sub("\\\.bam(\\\.gz)?", ".dupChim.txt", file), filterdupChim=TRUE)
}
\arguments{
  \item{file}{
  Path to the BAM file.
}
  \item{maxReadsOfSameBreak}{
  The maximum allowed number of artifact chimeric reads sharing a false positive
  breakpoint. If the number of reads sharing the same breakpoint exceeds this
  number, these reads are not recognized as artifact chimeric reads. Reads
  marked as PCR or optical duplicates are excluded from the calculation. For 
  paired-end sequencing, a read pair of artifact chimeric fragments may both 
  contain the artifact breakpoints; thereby, the defalut is set to 2.
}
  \item{minMapBase}{
  The minimum required length (bp) of a short complementary mapping for an
  artifact chimeric read. Artifact chimeric reads are derived from the
  combination of two single-stranded DNA fragments linked by short reverse
  complementary regions (SRCR). Reads with SRCR shorter than this length are 
  not recognized as artifact chimeric reads. Note: sequence errors and mutations
  might influence the detection of the existence and length of SRCR. Suggested
  range: 0-3. When it is set to 0 or any value below 1, this step will be 
  skipped.
}
  \item{threads}{
  Number of threads to use. Multi-threading can speed up the process.
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
  \item{FFPEReadsFile}{
  Path of the output txt file with artifact chimeric read names.
}
  \item{dupChimFile}{
  Path of the output txt file with supplementary reads that are marked as PCR
  or optical duplicates.
}
  \item{filterdupChim}{
  Filter PCR or optical duplicates of all chimeric reads when set to true.
  These reads may contain duplicates of artifact chimeric reads; therefore, it
  is recommended to also remove these reads.
}
}
\details{
The next-generation sequencing (NGS) reads from formalin-fixed paraffin-embedded
(FFPE) samples contain numerous artifact chimeric reads, which can lead to a
large number of false positive structural variation (SV) calls. This function
finds and filters these artifact chimeric reads. An index file is also generated
for the created filtered BAM file.
}
\value{
The file name of the created destination file.
}
\author{
Lanying Wei <lanying.wei@uni-muenster.de>
}
\seealso{
\code{\link{FilterFFPE}}, \code{\link{findArtifactChimericReads}},
\code{\link{filterBamByReadNames}}
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
