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
    return(c("#000000", "#FFFFFF", hcl.colors(n - 2, "Blues", rev = TRUE)))
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
    # hoge <- viridis::plasma(100, alpha = 0.3)
    # maskWithoutExp <- mask[(mask$value == 1 & df$value <= threshold),]
    xlim <- max(df$x)
    ylim <- max(df$y)
    zlim <- max(df$z)
    maxLen <- max(xlim, ylim, zlim)
    plotResult <- plot_ly(
        df[df$value > threshold, ],
        x=~x,
        y=~y,
        z=~z,
        #color = ~value,
        #colors = viridis::plasma(1000000, alpha = 0.3),
        # alpha = exp_opacity * 0.01,
        alpha_stroke = 0,
        size = I(exp_size),
        marker = list(
        #color = rgb(0.267, 0.071, 0.341, 0.3),
        color =~value,
        # autocolorscale = FALSE,
        # cmin = clim[1],
        cmin = min(df$value),
        cmax = clim[2],
        # showscale = FALSE,
        showscale = TRUE,
        # opacity = exp_opacity * 0.01,
        opacity = exp_opacity * 0.01,
        colorscale = "Viridis"
        # colorscale = list(list(0, "#F4FAFE"), list(0.5, "#7FABD3"), list(1, "#273871"))
        # colorscale = list(list(0, "#4B005580"), list(0.5, "#009B0580"), list(1, "#FDE33380"))
        #colorscale = list(list(0, "#4401544D"), list(0.5, "#21908C4D"), list(1, "#FDE7254D"))
        # colorscale = list(list(0, "#44015466"), list(0.5, "#21908C66"), list(1, "#FDE72566"))
        )
        #colorscale  = c(rgb(0.267, 0.071, 0.341), rgb(0.984,0.894,0.161)))

    ) %>%
        add_markers()%>%
        # add_markers(
        #     data = maskWithoutExp,
        #     x = ~x,
        #     y=~y,
        #     z=~z,
        #     inherit = F,
        #     alpha = 0.2,
        #     alpha_stroke = 0,
        #     size = I(20),
        #     color = I("#12121233")
        # ) %>%
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
            aspectratio = list(x=xlim / maxLen * aspX, y=ylim/maxLen * aspY, z=zlim/maxLen * aspZ)
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
                inherit = F,
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
    mask <- tomoseqr:::matrixToDataFrame(tomoObj$mask)
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
    xLeftShift <- 1:(dimArrS[1] - 2)
    xRightShift <- 3:dimArrS[1]
    yCenter <- 2:(dimArrS[2] - 1)
    yLeftShift <- 1:(dimArrS[2] - 2)
    yRightShift <- 3:dimArrS[2]
    zCenter <- 2:(dimArrS[3] - 1)
    zLeftShift <- 1:(dimArrS[3] - 2)
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