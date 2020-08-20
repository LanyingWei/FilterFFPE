findArtifactChimericReads <- function(file, maxReadsOfSameBreak=2,
                                      minMapBase=1, threads=1,
                                      FFPEReadsFile=sub("\\.bam(\\.gz)?",
                                                        ".FFPEReads.txt",
                                                          file),
                                      dupChimFile=sub("\\.bam(\\.gz)?",
                                                      ".dupChim.txt", file)) {
    if(missing(file)) stop("file is required")
    if(maxReadsOfSameBreak < 1) {
        stop("maxReadsOfSameBreak should be larger than 1")
    }
    if(minMapBase < 1) {
        stop("minMapBase should be larger than 0")
    }
    message('Finding reads with supplementary alignments...')
    supAlign <- findAllSupplymentaryReads(bamFilePath = file,
                                          mapqFilter = 1,
                                          writeDupToFile = TRUE,
                                          dupFile = dupChimFile,
                                          threads = threads)
    message(paste0('Found ', length(supAlign$qname), " reads."))
    message("Excluding chimeric reads from real SVs...")
    supAlignNotSV <- filterReadsOfSV(supAlign,
                                     maxReadsOfSameBreak = maxReadsOfSameBreak)
    message(paste0(length(supAlignNotSV$qname), " reads remained."))
    readNames <-  unique(supAlignNotSV$qname)
    message("Filtering reads based on short complementary sequences...")
    allChimericReadsNotSV <- findAllReadsWithID(bamFilePath = file,
                                                readNames = readNames,
                                                threads = threads)
    artifactChimericReads <-
        filterAllReadsWithoutShortMapping(chimericReads = allChimericReadsNotSV,
                                          threads = threads)
    FFPEReads <- unique(artifactChimericReads$qname)
    write(FFPEReads, FFPEReadsFile)
    message(paste0(length(FFPEReads), " reads remained and reported."))
    return(FFPEReads)
}


filterBamByReadNames <- function(file, readsToFilter, index=file,
                                 destination=sub("\\.bam(\\.gz)?",
                                                 ".FilterFFPE.bam", file),
                                 overwrite=FALSE) {
    if(missing(file)) stop("file is required")
    if(missing(readsToFilter)) stop("readsToFilter is required")
    if (overwrite == FALSE & file.exists(destination)) {
        stop("Destination file already exists!")
    }
    filterReadsByNames <- function(readsToFilter) {
        list(keepQname = function(x) ! x$qname %in% readsToFilter)
    }
    filter <- FilterRules(filterReadsByNames(readsToFilter = readsToFilter))
    outFile <- filterBam(file = file, destination = destination, index = index,
                         filter = filter,
                         indexDestination = TRUE,
                         overwrite = overwrite,
                         param = ScanBamParam(what="qname"))
    return (outFile)
}

FFPEReadFilter <- function(file, maxReadsOfSameBreak=2, minMapBase=1,
                           threads=1, index=file,
                           destination=sub("\\.bam(\\.gz)?",
                                           ".FilterFFPE.bam", file),
                           overwrite=FALSE,
                           FFPEReadsFile=sub("\\.bam(\\.gz)?",
                                           ".FFPEReads.txt", file),
                           dupChimFile=sub("\\.bam(\\.gz)?",
                                           ".dupChim.txt", file),
                           filterdupChim=TRUE) {
    if(missing(file)) stop("file is required")
    if(maxReadsOfSameBreak < 1) {
        stop("maxReadsOfSameBreak should be larger than 1")
    }
    if(minMapBase < 1) {
        stop("minMapBase should be larger than 0")
    }
    if (overwrite == FALSE & file.exists(destination)) {
        stop("Destination file already exists!")
    }
    artifactReads <-
        findArtifactChimericReads(file = file,
                                  maxReadsOfSameBreak = maxReadsOfSameBreak,
                                  minMapBase = minMapBase,
                                  threads = threads,
                                  FFPEReadsFile = FFPEReadsFile,
                                  dupChimFile = dupChimFile)
    if (filterdupChim == TRUE) {
        dupChim <- readLines(dupChimFile)
        message(paste0("Filtering ", length(unique(artifactReads)),
                       " artifact chimeric reads and ", length(unique(dupChim)),
                       " PCR or optical duplicates of all chimeric reads from ",
                       "BAM file..."))
        readsToFilter <- unique(c(artifactReads, dupChim))
    } else {
        readsToFilter <- artifactReads
        message(paste0("Filtering ", length(readsToFilter),
                       " artifact chimeric reads on BAM file..."))
    }

    outFile <- filterBamByReadNames(file = file,
                                    readsToFilter = readsToFilter,
                                    destination = destination,
                                    overwrite = overwrite,
                                    index = index)
    message("Done.")
    return(outFile)
}


