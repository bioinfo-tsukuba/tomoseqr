
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