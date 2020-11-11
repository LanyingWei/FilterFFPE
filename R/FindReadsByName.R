findReadsWithID <- function(bamFile, readNames, what=scanBamWhat(),
                            which=IRangesList(), isSupplementaryAlignment=NA,
                            reverseComplement=NA) {

    message('Finding reads with read name...')
    if (!'qname' %in% what) {
        what = c('qname', what)
    }
    f <- scanBamFlag(isNotPassingQualityControls = FALSE,
                     isDuplicate = FALSE,
                     isSupplementaryAlignment = isSupplementaryAlignment)

    p <- ScanBamParam(flag = f,
                      simpleCigar = FALSE,
                      reverseComplement = reverseComplement,
                      what = what,
                      which = which)

    bam <- scanBam(bamFile, param = p)[[1]]
    bam <- lapply(bam, function(x) x[bam$qname %in% readNames])
    message(paste("Found", length(bam$qname), "reads."))
    return(bam)
}

findAllReadsWithID <- function(bamFilePath, readNames,
                               what=c("qname", "flag", "cigar"),
                               isSupplementaryAlignment=NA,
                               reverseComplement=FALSE, threads=1) {

    chunkSize <- 1e+07
    if (!'qname' %in% what) {
        what = c('qname', what)
    }

    bamFile <- BamFile(bamFilePath)
    seqInfo <- scanBamHeader(bamFile)$targets
    chrRegions <- GRanges(names(seqInfo), IRanges(1, seqInfo))
    chrRegions <- tile(chrRegions, width = chunkSize)
    chrRegions <- unlist(chrRegions)
    cl <- makeCluster(threads)
    registerDoParallel(cl)
    i <- 0
    allReads <- foreach(i=seq_along(chrRegions),
                        .combine = function(x, y) Map("c", x, y),
                        .inorder = TRUE,
                        .verbose = FALSE,
                        .packages = c("Biostrings", "Rsamtools"),
                        .export = c("findReadsWithID")) %dopar% {
                            findReadsWithID(bamFile = bamFile,
                                            readNames = readNames,
                                            what = what,
                                            which = chrRegions[i],
                                            isSupplementaryAlignment = NA,
                                            reverseComplement = FALSE)
                        }
    stopCluster(cl)
    return(allReads)
}
