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
    message('Finding reads with supplementary alignments...')
    supAlign <- findAllSupplymentaryReads(bamFilePath = file,
                                          writeDupToFile = TRUE,
                                          dupFile = dupChimFile,
                                          threads = threads)
    message(paste0('Found ', length(unique(supAlign$qname)), " read names"))
    message("Excluding chimeric reads from real SVs...")
    message("Filtering reads based on the number of shared breakpoints...")
    supAlignNotSV <- filterReadsOfSV(supAlign,
                                     maxReadsOfSameBreak = maxReadsOfSameBreak)
    readNames <-  unique(supAlignNotSV$qname)
    message(paste0(length(readNames), " read names remained."))
    if (minMapBase < 1) {
        message("minMapBase < 1, skip filtering with minMapBase.")
        message(paste0(length(readNames), " read names remained and reported."))
        write(readNames, FFPEReadsFile)
        return(readNames)
    } else {
        message(paste0("Filtering reads based on the existence and the length ",
                       "of short reverse complementary regions..."))
        allChimericReads <- findAllReadsWithID(bamFilePath = file,
                                               readNames = readNames,
                                               threads = threads)
        artifactChimericReads <-
            filterAllReadsWithoutShortMapping(chimericReads = allChimericReads,
                                              threads = threads, 
                                              minMapBase = minMapBase)
        FFPEReads <- unique(artifactChimericReads)
        write(FFPEReads, FFPEReadsFile)
        message(paste0(length(FFPEReads), " read names remained and reported."))
        return(FFPEReads)
    }
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
        dupChim <- dupChim[dupChim != ""]
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


