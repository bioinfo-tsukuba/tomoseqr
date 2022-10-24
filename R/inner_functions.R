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
    return(c("#000000", "#FFFFFF", hcl.colors(n - 2, "Blues", rev = TRUE)))
}

#' @importFrom stringr str_c
make2DPlot <- function (
    sectionNumber,
    basePlot,
    expDf,
    along,
    xAxis,
    yAxis,
    main
) {
    xAxisParameter <- sym(xAxis)
    yAxisParameter <- sym(yAxis)
    value <- sym("value")
    print(
        basePlot +
            geom_tile(
                data=expDf[expDf[[along]] == sectionNumber,],
                aes(
                    x = !!xAxisParameter,
                    y = !!yAxisParameter,
                    fill=!!value
                )
            ) +
            ggtitle(str_c(main, " : ", sectionNumber))
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
    basePlot,
    expDf,
    dimOrder,
    along,
    xAxis,
    yAxis,
    main,
    forShiny
) {
    if (forShiny == FALSE) {
        message("Plotting", appendLF=FALSE)
    }
    ind <- seq_len(dimOrder[[along]][3])
    lapply(
        ind,
        make2DPlot,
        basePlot = basePlot,
        expDf = expDf,
        along = along,
        xAxis = xAxis,
        yAxis = yAxis,
        main = main
    )
    if (forShiny == TRUE) {
        incProgress(1, detail = "Converting to GIF. It takes long time...")
    } else {
        message("")
        message("Converting to GIF. It takes long time...")
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


show3d <- function(
    df,
    mask,
    threshold,
    clim,
    xlab,
    ylab,
    zlab,
    addMask,
    exp_size,
    exp_opacity,
    mask_size,
    mask_opacity,
    mask_color,
    aspX,
    aspY,
    aspZ
) {
    xlim <- max(df$x)
    ylim <- max(df$y)
    zlim <- max(df$z)
    maxLen <- max(xlim, ylim, zlim)
    plotResult <- plot_ly(
        df[df$value > threshold, ],
        x=~x,
        y=~y,
        z=~z,
        alpha_stroke = 0,
        size = I(exp_size),
        marker = list(
        color =~value,
        cmin = min(df$value),
        cmax = clim[2],
        showscale = TRUE,
        opacity = exp_opacity * 0.01,
        colorscale = "Viridis"
        )

    ) %>%
        add_markers() %>%
        layout(
            autosize = TRUE,
            showlegend = FALSE,
            scene = list(
            xaxis=list(
                title = xlab,
                range = c(0, xlim),
                gridcolor="#ffffff",
                zerolinecolor ="#ffffff"
            ),
            yaxis = list(
                title = ylab,
                range=c(0,ylim),
                gridcolor="#ffffff",
                zerolinecolor ="#ffffff"
            ),
            zaxis = list(
                title = zlab,
                range=c(0,zlim),
                gridcolor="#ffffff",
                zerolinecolor ="#ffffff"
            ),
            aspectratio = list(
                x=xlim / maxLen * aspX,
                y=ylim/maxLen * aspY,
                z=zlim/maxLen * aspZ
                )
            ),
            paper_bgcolor = "#000000",
            font = list(color="#ffffff")
        ) %>%
        config(
            toImageButtonOptions = list(
                format = "svg",
                filename = "myplot",
                width = 1000,
                height = 1000,
                scale = 3
            )
        )
    if (addMask == TRUE) {
        maskWithoutExp <- mask[(mask$value == 1 & df$value <= threshold),]
        plotResult %>%
            add_markers(
                data = maskWithoutExp,
                x = ~x,
                y=~y,
                z=~z,
                inherit = FALSE,
                alpha = mask_opacity * 0.01,
                # alpha = 0.2,
                alpha_stroke = 0,
                size = I(mask_size),
                color = I(mask_color)
            ) %>%
            return()
    } else {
        return(plotResult)
    }
}

showExp <- function(
    tomoObj,
    geneID,
    threshold,
    xlab,
    ylab,
    zlab,
    addMask,
    exp_size,
    exp_opacity,
    mask_size,
    mask_opacity,
    mask_color,
    aspX,
    aspY,
    aspZ
    ) {
    df <- toDataFrame(tomoObj, geneID)
    mask <- matrixToDataFrame(tomoObj[["mask"]])
    show3d(
        df,
        mask,
        threshold,
        c(threshold, max(df$value)),
        xlab = xlab,
        ylab = ylab,
        zlab = zlab,
        addMask = addMask,
        exp_size = exp_size,
        exp_opacity = exp_opacity,
        mask_size = mask_size,
        mask_opacity = mask_opacity,
        mask_color = mask_color,
        aspX,
        aspY,
        aspZ
    )
}
plotsurface <- function(arr) {
    dimarr <- dim(arr)
    arrSurZero <- array(
        0,
        dim = c(
            dimarr[1] + 2,
            dimarr[2] + 2,
            dimarr[3] + 2
        )
    )
    dimArrS <- dim(arrSurZero)
    xCenter <- 2:(dimArrS[1] - 1)
    xLeftShift <- seq_len(dimArrS[1] - 2)
    xRightShift <- 3:dimArrS[1]
    yCenter <- 2:(dimArrS[2] - 1)
    yLeftShift <- seq_len(dimArrS[2] - 2)
    yRightShift <- 3:dimArrS[2]
    zCenter <- 2:(dimArrS[3] - 1)
    zLeftShift <- seq_len(dimArrS[3] - 2)
    zRightShift <- 3:dimArrS[3]

    arrSurZero[xCenter, yCenter, zCenter] <- arr
    return(
        arr * (
            arrSurZero[xLeftShift, yCenter, zCenter] *
            arrSurZero[xRightShift, yCenter, zCenter] *
            arrSurZero[xCenter, yLeftShift, zCenter] *
            arrSurZero[xCenter, yRightShift, zCenter] *
            arrSurZero[xCenter, yCenter, zLeftShift] *
            arrSurZero[xCenter, yCenter, zRightShift]
            == 0
        )
    )
}

#' @importFrom ggplot2
#'      xlab
#'      ylab
#'      xlim
#'      ylim
makeBasePlot <- function(
    expDf,
    xAxis,
    yAxis, 
    xMax,
    yMax,
    xAsp,
    yAsp,
    xlabel,
    ylabel,
    zlim
) {
    xAxisParameter <- sym(xAxis)
    yAxisParameter <- sym(yAxis)
    value <- sym("value")
    return(
        ggplot(
            expDf,
            aes(
                x = !!xAxisParameter,
                y = !!yAxisParameter,
                fill = !!value
            )
        ) +
            scale_fill_gradientn(
                colors = hcl.colors(50, "Blues", rev = TRUE),
                limits = zlim
            ) +
            xlim(0.0, xMax) +
            ylim(0.0, yMax) +
            xlab(xlabel) +
            ylab(ylabel) +
            theme_minimal() +
            theme(
                plot.background= element_rect(fill="black"),
                text = element_text(color = "white", size=20),
                axis.line.x.bottom = element_line(color = "white"),
                axis.line.y.left = element_line(color = "white"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks=element_line(colour = "white"),
                axis.text=element_text(colour = "white", size=20)
                ) +
            theme(aspect.ratio = yAsp / xAsp)
    )
}