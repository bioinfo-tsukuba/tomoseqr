#' @importFrom methods is
CheckParameters <- function(tomoObj, query) {
    if (is(tomoObj, "tomoSeq") == FALSE) {
        stop(
            paste(
                "invalid class:",
                class(tomoObj),
                "\nFirst argument must be tomoSeq class object."
            )
        )
    }
    if (is.element(query, tomoObj$geneList) %>% min() == 0) {
        queryNotIn <- query[is.element(query, tomoObj$geneList) == FALSE]
        queryNotInStr <- paste(queryNotIn, collapse=", ")
        stop(paste(queryNotInStr, ' is not in data.', sep=''))
    }
}

## This function is used in Estimate3dExpression().
GetGeneExpression <- function (tomoSeqData, geneID) {
    retvalMatrix <- tomoSeqData[tomoSeqData[, 1] == geneID, ]
    retvalMatrix <- retvalMatrix[, -1] %>% as.matrix()
    return(retvalMatrix)
}

RepMat <- function (targetVector, nTimesRepeat) {
    lenX <- length(targetVector) * nTimesRepeat[1]
    lenY <- nTimesRepeat[2]
    lenZ <- nTimesRepeat[3]
    rep2d <- matrix(targetVector, nrow=lenX, ncol=lenY, byrow = F)
    rep3d <- array(NA, dim=c(lenX, lenY, lenZ))
    for (i in 1:lenZ) {
        rep3d[, , i] <- rep2d
    }
    return(rep3d)
}

SingleEstimate <- function (
    X,
    Y,
    Z,
    geneID,
    mask,
    normCount,
    normMask,
    numIter
) {
    sumX <- X[, -1] %>% colSums()
    sumY <- Y[, -1] %>% colSums()
    sumZ <- Z[, -1] %>% colSums()

    x0 <- X %>% GetGeneExpression(geneID)
    y0 <- Y %>% GetGeneExpression(geneID)
    z0 <- Z %>% GetGeneExpression(geneID)
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

    for (i in 1:numIter) {
        xa <- a %>% apply(1, sum)
        a <- a * RepMat(x / xa, c(1, dim(mask)[2], dim(mask)[3]))
        a[is.nan(a)] <- 0
        ya <- a %>% apply(2, sum)
        a <- a * aperm(
            RepMat(y / ya, c(1, dim(mask)[1], dim(mask)[3])),
            perm = c(2, 1, 3)
        )
        a[is.nan(a)] <- 0
        za <- a %>% apply(3, sum)
        a <- a * aperm(
            RepMat(z / za, c(1, dim(mask)[1], dim(mask)[2])),
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

ColFunc <- function (n) {
    return(c("#000000", hcl.colors(n - 1, "Blues", rev = TRUE)))
}

BasePlot <- function (
    sourceArray,
    sectionNumber,
    main,
    xlab,
    ylab,
    aspectRatio
) {
    filled.contour(
        sourceArray[, , sectionNumber],
        main=paste(main, "_", sectionNumber, sep=""),
        xlab=xlab,
        ylab=ylab,
        asp=aspectRatio,
        frame.plot=F,
        zlim = c(-1, max(sourceArray)),
        nlevel = 50,
        color.palette = ColFunc
    )
}

MakePlotArray <- function (
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
AnimateForGIF <- function (
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

    plotArray <- MakePlotArray(
        reconstArray=reconstArray,
        maskArray=maskArray,
        zlim=zlim,
        type=type
    )
    message("generating", appendLF=FALSE)
    for (i in seq_along(plotArray[1, 1, ])) {
        message(".", appendLF=FALSE)
        BasePlot(
            plotArray,
            sectionNumber=i,
            main,
            xlab,
            ylab,
            asp
        )
    }
    message("")
}

PlotForImageViewer <- function (
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

    plotArray <- MakePlotArray(
        reconstArray=reconstArray,
        maskArray=maskArray,
        zlim=zlim,
        type=type
    )
    BasePlot(
        plotArray,
        sectionNumber=sectionNumber,
        main,
        xlab,
        ylab,
        asp
    )
}

#' @export
print.tomoSeq <- function (tomoObj) {
    cat("Gene list:\n")
    cat(names(tomoObj[["results"]]))
}


MatrixToDataFrame <- function (reconst) {
    vecReconst <- as.vector(reconst)
    dim <- dim(reconst)
    xlen <- dim[1]
    ylen <- dim[2]
    zlen <- dim[3]
    xIndex <- rep(1:xlen, ylen * zlen)
    yIndex <- 1:ylen %>%
        sapply(function (p) {rep(p, xlen)}) %>%
        rep(zlen)
    zIndex <- 1:zlen %>%
        sapply(function (p) {rep(p, xlen * ylen)}) %>%
        as.vector()
    data.frame(x=xIndex, y=yIndex, z=zIndex, value=vecReconst) %>%
        return()
}