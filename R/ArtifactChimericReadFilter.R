filterReadsOfSV <- function(supAlign, maxReadsOfSameBreak) {
    chrs <- unique(supAlign$rname)
    for (i in seq_along(chrs)) {
        chr <- chrs[i]
        temp <- lapply(supAlign, function(x) x[which(supAlign$rname == chr)])
        startPosFreq <- table(temp$pos)
        endPosFreq <- table(temp$pos + temp$qwidth)
        startPos <- names(startPosFreq[startPosFreq <= maxReadsOfSameBreak])
        endPos <- names(endPosFreq[endPosFreq <= maxReadsOfSameBreak])
        keepIndex <-
            which(temp$pos %in% startPos & (temp$pos + temp$qwidth) %in% endPos)
        temp <- lapply(temp, function(x) x[keepIndex])
        if (i == 1) {
            res = temp
        } else {
            res <- Map("c", res, temp)
        }
    }
    return(res)
}

filterReadsWithoutShortMapping <- function(chimericReads, minMapBase) {
    
    firstRead <-
        vapply(chimericReads$flag, function(x) as.integer(intToBits(x))[7],
               FUN.VALUE = numeric(1))
    isSupplementary <-
        vapply(chimericReads$flag, function(x) as.integer(intToBits(x))[12],
               FUN.VALUE = numeric(1))
    cigarChar <-
        lapply(chimericReads$cigar, function(x)
            unlist(strsplit(as.character(x), split = "[0-9]+"))[-1])
    cigarNum <-
        lapply(chimericReads $cigar, function(x)
            unlist(strsplit(as.character(x), split = "[A-Z]+")))
    
    chimericReads[[length(chimericReads ) + 1]] <- firstRead
    chimericReads[[length(chimericReads ) + 1]] <- isSupplementary
    chimericReads[[length(chimericReads ) + 1]] <- cigarChar
    chimericReads[[length(chimericReads ) + 1]] <- cigarNum
    
    names(chimericReads )[(length(chimericReads ) - 3):
                              length(chimericReads ) ] <-
        c("firstRead", "isSupplementary", "cigarChar", "cigarNum")
    
    readNames <- unique(chimericReads$qname)
    keepedReadNames <- NULL
    
    for (readName in readNames) {
        shortMapping <- FALSE
        missInfo <- c(FALSE, FALSE)
        readIndex <- which(chimericReads$qname == readName)
        if (length(readIndex) >= 2) {
            allReads <- lapply(chimericReads, function(x) x[readIndex])
            for (i in c(0, 1)) {
                read <- (allReads$firstRead == i) & !allReads$isSupplementary
                readSup <- (allReads$firstRead == i) & allReads$isSupplementary
                softClipLen <- NA
                supMatchLen <- NA
                if (any(read)) {
                    softClipLen <-
                        suppressWarnings(max(as.integer(
                            allReads$cigarNum[read][[1]][
                                allReads$cigarChar[read][[1]] == "S"])))
                }
                if (any(readSup)) {
                    supMatchLen <-
                        suppressWarnings(sum(as.integer(
                            allReads$cigarNum[readSup][[1]][
                                allReads$cigarChar[readSup][[1]] == "M"])))
                }
                if (is.na(softClipLen) | is.na(supMatchLen)) {
                    missInfo[i] <- TRUE
                }
                if(!is.na(softClipLen) & !is.na(supMatchLen)) {
                    if (supMatchLen - softClipLen >= minMapBase) {
                        shortMapping <- TRUE
                    }
                }
            }
        } else {
            shortMapping <- TRUE
        }
        if(all(missInfo)) {
            shortMapping <- TRUE
        }
        if (shortMapping) {
            keepedReadNames <- c(keepedReadNames, readName)
        }
    }
    
    return (unique(keepedReadNames))
}

filterAllReadsWithoutShortMapping <- function(chimericReads, minMapBase, 
                                              threads) {
    nReads <- 5e+4
    starts <- seq(1, length(chimericReads$qname), by = nReads)
    ends <- starts + nReads - 1
    ends[length(ends)] <- length(chimericReads$qname)
    nameOrder <- order(chimericReads$qname)
    chimericReads <- lapply(chimericReads, function(x) x[nameOrder])
    starts <- seq(1, length(chimericReads$qname), by = nReads)
    ends <- starts + nReads - 1
    ends[length(ends)] <- length(chimericReads$qname)
    
    cl <- makeCluster(threads)
    registerDoParallel(cl)
    res <- foreach(start=starts,
                   end=ends,
                   .combine = "c",
                   .verbose = FALSE,
                   .packages = c("Biostrings", "Rsamtools","BiocGenerics"),
                   .export = c("filterReadsWithoutShortMapping")
    ) %dopar% {
        filterReadsWithoutShortMapping(chimericReads = 
                                           lapply(chimericReads, 
                                                  function(x) x[start:end]), 
                                       minMapBase = minMapBase)
    }
    stopCluster(cl)
    return(res)
}

