findSupplymentaryReads <- function(bamFile, which=IRangesList(),
                                   what=scanBamWhat(), mapqFilter=1,
                                   writeDupToFile=FALSE, dupFile=NA,
                                   append=TRUE) {
    message('Finding reads with supplymentary alignments...')
    f <- scanBamFlag(isNotPassingQualityControls = FALSE,
                     isDuplicate = FALSE,
                     isSupplementaryAlignment = TRUE)

    p <- ScanBamParam(flag = f,
                      simpleCigar = FALSE,
                      what = what,
                      which = which,
                      mapqFilter = mapqFilter)

    supAlign <- scanBam(bamFile, param = p)[[1]]
    message(paste("Found ", length(supAlign$qname),
                  " reads with supplymentary alignments."))

    f <- scanBamFlag(isNotPassingQualityControls = FALSE,
                     isDuplicate = TRUE,
                     isSupplementaryAlignment = NA)

    p <- ScanBamParam(flag = f,
                      simpleCigar = FALSE,
                      what = "qname",
                      which = which,
                      mapqFilter = mapqFilter)

    dupReadName <- unique(scanBam(bamFile, param = p)[[1]]$qname)
    dupIndex <- which(supAlign$qname %in% dupReadName)
    if (length(dupIndex) > 0) {
        if (writeDupToFile) {
            write(unique(supAlign$qname[dupIndex]), dupFile, append = append)
        }
        message(paste("Found", length(dupIndex), "duplicates."))
        supAlign <- lapply(supAlign, function(x) x[-dupIndex])
        message(paste("Duplicates removed,", length(supAlign$qname),
                      "reads with supplymentary alignment remained."))
    }
    return(supAlign)
}

findAllSupplymentaryReads <- function(bamFilePath, writeDupToFile=FALSE,
                                      dupFile=NA, threads=1, mapqFilter=1) {

    bamFile <- BamFile(bamFilePath)
    seqInfo <- scanBamHeader(bamFile)$targets
    chrRegions <- GRanges(names(seqInfo), IRanges(1, seqInfo))

    cl <- makeCluster(threads)
    registerDoParallel(cl)
    what=c("qname", "rname", "pos", "qwidth")
    i <- 0
    if (writeDupToFile == TRUE) {
        write("", dupFile)
    }
    supAlign <- foreach(i=seq_along(chrRegions),
                        .combine = function(x, y) Map("c", x, y),
                        .inorder = TRUE,
                        .verbose = FALSE,
                        .packages = c("Biostrings", "Rsamtools"),
                        .export = c("findSupplymentaryReads")
    ) %dopar% {
        findSupplymentaryReads(bamFile = bamFile,
                               which = chrRegions[i],
                               what = what,
                               mapqFilter = mapqFilter,
                               writeDupToFile = writeDupToFile,
                               dupFile = dupFile,
                               append = TRUE)
                        }
    stopCluster(cl)
    return(supAlign)
}
