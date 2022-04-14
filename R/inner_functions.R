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
    normCount,
    normMask,
    numIter
) {
    sumX <- dataX[, -1] %>% colSums()
    sumY <- dataY[, -1] %>% colSums()
    sumZ <- dataZ[, -1] %>% colSums()

    x0 <- dataX %>% getGeneExpression(geneID)
    y0 <- dataY %>% getGeneExpression(geneID)
    z0 <- dataZ %>% getGeneExpression(geneID)
    xLen <- length(x0)
    yLen <- length(y0)
    zLen <- length(z0)

    maskX <- mask %>% apply(1, sum)
    maskY <- mask %>% apply(2, sum)
    maskZ <- mask %>% apply(3, sum)

    if (normCount == "countSum") {
        x <- x0 / sumX
        y <- y0 / sumY
        z <- z0 / sumZ
    }
    
    if (is.list(normCount)) {
        x <- x0 / normCount[["x"]]
        y <- y0 / normCount[["y"]]
        z <- z0 / normCount[["z"]]
    }

    if (normCount == "none") {
        x <- x0
        y <- y0
        z <- z0
    }

    if (normMask) {
        x <- x * maskX
        y <- y * maskY
        z <- z * maskZ
    }

    x[1, is.nan(x)] <- 0
    y[1, is.nan(y)] <- 0
    z[1, is.nan(z)] <- 0
    x[is.infinite(x)] <- max(x[x < Inf])
    y[is.infinite(y)] <- max(z[y < Inf])
    z[is.infinite(z)] <- max(y[z < Inf])

    xRaw <- sumX
    yRaw <- sumY
    zRaw <- sumZ

    m <- mean(c(sum(xRaw), sum(yRaw), sum(zRaw)))

    x <- x / sum(x) * m
    y <- y / sum(y) * m
    z <- z / sum(z) * m
    a <- mask

    er <- c()

    for (i in seq_len(numIter)) {
        xa <- a %>% apply(1, sum)
        a <- a * repMat(x / xa, c(1, dim(mask)[2], dim(mask)[3]))
        a[is.nan(a)] <- 0
        ya <- a %>% apply(2, sum)
        a <- a * aperm(
            repMat(y / ya, c(1, dim(mask)[1], dim(mask)[3])),
            perm = c(2, 1, 3)
        )
        a[is.nan(a)] <- 0
        za <- a %>% apply(3, sum)
        a <- a * aperm(
            repMat(z / za, c(1, dim(mask)[1], dim(mask)[2])),
            perm = c(2, 3, 1)
        )
        a[is.nan(a)] <- 0
        er <- append(er, sum((xa - x)^2) + sum((ya - y)^2) + sum((za - z)^2))
    }

    retList <- list(
        "geneID" = geneID,
        "reconst" = a,
        "errFunc" = er,
        "x" = x[1, ],
        "y" = y[1, ],
        "z" = z[1, ]
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
#' withProgress
#' incProgress
animateForGIF <- function (
    reconstArray,
    maskArray,
    main,
    xlab,
    ylab,
    zlim,
    aspectRatio,
    type
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
    message("Plotting", appendLF=FALSE)
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
    message("")
    message("Converting to GIF...")
}

plotForImageViewer <- function (
    reconstArray,
    maskArray,
    main,
    xlab,
    ylab,
    zlim,
    aspectRatio,
    type,
    sectionNumber
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
    basePlot(
        plotArray,
        sectionNumber=sectionNumber,
        main,
        xlab,
        ylab,
        asp
    )
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