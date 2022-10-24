#' @importFrom methods is
#' @importFrom dplyr %>%
#' @importFrom stringr str_c
checkParameters <- function(tomoObj, query) {
    if (is(tomoObj, "tomoSeq") == FALSE) {
        stopMessage <- str_c(
            "invalid class: ",
            class(tomoObj),
            "\nFirst argument must be tomoSeq class object."
        )
        stop(stopMessage)
    }
    if (is.element(query, names(tomoObj[["results"]])) %>% min() == 0) {
        queryNotIn <- query[is.element(query, tomoObj[["results"]]) == FALSE]
        queryNotInStr <- str_c(queryNotIn, sep=", ")
        stopMessage <- str_c(queryNotInStr, ' is not in data.')
        stop(stopMessage)
    }
}

## This function is used in estimate3dExpression().
#' @importFrom dplyr %>%
getGeneExpression <- function (tomoSeqData, geneID) {
    retvalMatrix <- tomoSeqData[as.vector(tomoSeqData[, 1] == geneID), ]
    retvalMatrix <- retvalMatrix[, -1] %>% as.matrix()
    return(retvalMatrix)
}

repMat <- function (targetVector, nTimesRepeat) {
    lenX <- length(targetVector) * nTimesRepeat[1]
    lenY <- nTimesRepeat[2]
    lenZ <- nTimesRepeat[3]
    rep2d <- matrix(targetVector, nrow=lenX, ncol=lenY, byrow = FALSE)
    rep3d <- array(rep2d, dim=c(lenX, lenY, lenZ))
    return(rep3d)
}

#' @importFrom dplyr %>%
singleEstimate <- function (
    geneID,
    dataX,
    dataY,
    dataZ,
    mask,
    normalize,
    numIter
) {
    sumPerSectionX <- dataX[, -1] %>% colSums()
    sumPerSectionY <- dataY[, -1] %>% colSums()
    sumPerSectionZ <- dataZ[, -1] %>% colSums()

    xk <- dataX %>% getGeneExpression(geneID)
    yk <- dataY %>% getGeneExpression(geneID)
    zk <- dataZ %>% getGeneExpression(geneID)
    xLen <- length(xk)
    yLen <- length(yk)
    zLen <- length(zk)

    maskSumX <- mask %>% apply(1, sum)
    maskSumY <- mask %>% apply(2, sum)
    maskSumZ <- mask %>% apply(3, sum)

    if (normalize) {
        xk<- xk / sumPerSectionX * maskSumX
        yk<- yk / sumPerSectionY * maskSumY
        zk<- zk / sumPerSectionZ * maskSumZ
    }

    xk[1, is.nan(xk)] <- 0
    yk[1, is.nan(yk)] <- 0
    zk[1, is.nan(zk)] <- 0
    xk[is.infinite(xk)] <- max(xk[xk < Inf])
    yk[is.infinite(yk)] <- max(zk[yk < Inf])
    zk[is.infinite(zk)] <- max(yk[zk < Inf])

    m <- mean(c(sum(sumPerSectionX), sum(sumPerSectionY), sum(sumPerSectionZ)))

    xk <- xk / sum(xk) * m
    yk <- yk / sum(yk) * m
    zk <- zk / sum(zk) * m
    reconstArray <- mask

    errFuncVal <- rep_len(0, numIter)
    dimMask <- dim(mask)
    for (i in seq_len(numIter)) {
        recArrX <- reconstArray %>% apply(1, sum)
        reconstArray <- reconstArray *
            repMat(xk / recArrX, c(1, dimMask[2], dimMask[3]))
        reconstArray[is.nan(reconstArray)] <- 0

        recArrY <- reconstArray %>% apply(2, sum)
        reconstArray <- reconstArray * aperm(
            repMat(yk / recArrY, c(1, dimMask[1], dimMask[3])),
            perm = c(2, 1, 3)
        )
        reconstArray[is.nan(reconstArray)] <- 0

        recArrZ <- reconstArray %>% apply(3, sum)
        reconstArray <- reconstArray * aperm(
            repMat(zk / recArrZ, c(1, dimMask[1], dimMask[2])),
            perm = c(2, 3, 1)
        )
        reconstArray[is.nan(reconstArray)] <- 0

        errFuncVal[i] <- sum((recArrX - xk)^2) +
            sum((recArrY - yk)^2) +
            sum((recArrZ - zk)^2)
    }

    retList <- list(
        "geneID" = geneID,
        "reconst" = reconstArray,
        "errFunc" = errFuncVal,
        "x" = xk[1, ],
        "y" = yk[1, ],
        "z" = zk[1, ]
    )
    return(retList)
}

#' @importFrom grDevices hcl.colors
colFunc <- function (n) {
    return(c("#000000", hcl.colors(n - 1, "Blues", rev = TRUE)))
}

#' @importFrom graphics filled.contour
#' @importFrom stringr str_c
basePlot <- function (
    sectionNumber,
    sourceArray,
    main,
    xlab,
    ylab,
    aspectRatio
) {
    filled.contour(
        sourceArray[, , sectionNumber],
        main=str_c(main, "_", sectionNumber),
        xlab=xlab,
        ylab=ylab,
        asp=aspectRatio,
        frame.plot=FALSE,
        zlim = c(-1, max(sourceArray)),
        nlevels = 50,
        color.palette = colFunc
    )
}

makePlotArray <- function (
    reconstArray,
    maskArray,
    zlim,
    type
) {
    if (type == "expression") {
        plotArray <- reconstArray
        plotArray[plotArray < zlim[1]] <- 0
        plotArray[zlim[2] < plotArray] <- 0
    } else if (type == "mask") {
        plotArray <- maskArray * 2 - 1
    } else if (type == "unite") {
        plotArray <- reconstArray
        plotArray[plotArray < zlim[1]] <- 0
        plotArray[zlim[2] < plotArray] <- 0
        plotArray <- plotArray + maskArray - 1
    } else {
        stop("'type' must be 'expression', 'mask' or 'unite'.")
    }

    return(plotArray)
}

#' @importFrom shiny
#' incProgress
animateForGIF <- function (
    reconstArray,
    maskArray,
    main,
    xlab,
    ylab,
    zlim,
    aspectRatio,
    type,
    forShiny
) {
    if (length(aspectRatio) < 2) {
        ## Dim of reconstructed matrix should be equal to
        ## that of mask.
        plotDim <- dim(reconstArray)
        asp <- plotDim[2] / plotDim[1]
    } else {
        asp <- aspectRatio[2] / aspectRatio[1]
    }

    plotArray <- makePlotArray(
        reconstArray=reconstArray,
        maskArray=maskArray,
        zlim=zlim,
        type=type
    )
    if (forShiny == FALSE) {
        message("Plotting", appendLF=FALSE)
    }
    ind <- seq_along(plotArray[1, 1, ])
    lapply(
        ind,
        basePlot,
        sourceArray=plotArray,
        main=main,
        xlab=xlab,
        ylab=ylab,
        aspectRatio=asp
    )
    if (forShiny == TRUE) {
        incProgress(1, detail = "Converting to GIF...")
    } else {
        message("")
        message("Converting to GIF...")
    }
}

#' @export
print.tomoSeq <- function (x, ...) {
    cat("Gene list:\n")
    cat(names(x[["results"]]))
}

#' @importFrom dplyr %>%
matrixToDataFrame <- function (reconst) {
    vecReconst <- as.vector(reconst)
    dim <- dim(reconst)
    xlen <- dim[1]
    ylen <- dim[2]
    zlen <- dim[3]
    xIndex <- rep(seq_len(xlen), ylen * zlen)
    yIndex <- seq_len(ylen) %>%
        vapply(function (p) {rep(p, xlen)}, FUN.VALUE = seq_len(xlen)) %>%
        rep(zlen)
    zIndex <- seq_len(zlen) %>%
        vapply(
            function (p) {rep(p, xlen * ylen)},
            FUN.VALUE = seq_len(xlen * ylen)
        ) %>%
        as.vector()
    data.frame(x=xIndex, y=yIndex, z=zIndex, value=vecReconst) %>%
        return()
}

isPeakGene <- function (v, threshold = 10) {
    maxV <- max(v)
    indexOfMaxV <- which.max(v)
    meanExeptMaxV <- mean(v[-1 * indexOfMaxV])
    if((maxV > threshold) && (maxV > 10 * meanExeptMaxV)) {
        return(c(maxV, meanExeptMaxV, indexOfMaxV))
    } else {
        return(c(maxV, meanExeptMaxV, 0))
    }
}

#' @importFrom tibble tibble
findAxialGenesInner <- function (targetData) {
    retMat <- targetData %>%
        select(-1) %>%
        apply(MARGIN = 1, FUN = isPeakGene) %>%
        t()
    tibble(
        geneID = as.vector(t(targetData[, 1])),
        max = retMat[, 1],
        meanExeptMax = retMat[, 2],
        isPeakGene = retMat[, 3]
    ) %>%
        return()
}

#' @importFrom tools R_user_dir
#' @importFrom BiocFileCache BiocFileCache
getTomoseqrCache <- function() {
    cache <- R_user_dir("tomoseqrCache", which="cache")
    BiocFileCache(cache)
}

#' @importFrom BiocFileCache bfcquery
#' @importFrom BiocFileCache bfcadd
#' @importFrom BiocFileCache bfcneedsupdate
#' @importFrom BiocFileCache bfcdownload
downloadData <- function (bfc, rname, URL, verbose=FALSE) {
    rid <- bfcquery(bfc, rname, "rname")$rid
    if (!length(rid)) {
        if (verbose) {
            message("Downloading file")
        }
        rid <- names(bfcadd(bfc, rname, URL))
    }
    if (isTRUE(bfcneedsupdate(bfc, rid))) {
        bfcdownload(bfc, rid)
    }
}