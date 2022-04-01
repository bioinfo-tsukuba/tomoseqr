#' Extract geneIDs to which hoge can be applied.
#' @param x A data.frame object containing a Tomo-seq data
#' for x-axis sections. The rows represent genes. The first column
#' contains gene IDs and the second and subsequent columns contain
#' gene expression levels in sections.
#' @param y A data.frame object containing a Tomo-seq data for y-axis
#' sections. The rows represent genes. The first column contains gene IDs and
#' the second and subsequent columns contain gene expression levels in sections.
#' @param z A data.frame object containing a Tomo-seq data for
#' z-axis sections. The rows represent genes. The first column
#' contains gene IDs and the second and subsequent columns contain
#' gene expression levels in sections.
#' @return A vector that contains genes which can be used for
#' `Estimate3dExpressions`.
#' @examples
#' data("testx", "testy", "testz")
#' ExtractGeneList(testx, testy, testz)
#' @importFrom dplyr %>%
#' @export
ExtractGeneList <- function (x, y, z) {
    xGene <- x[, 1] %>% t()
    yGene <- y[, 1] %>% t()
    zGene <- z[, 1] %>% t()
    xAndy <- intersect(xGene, yGene)
    return(intersect(xAndy, zGene))
}

#' Estimate 3d expressions
#' @param x A data.frame object containing a simulated Tomo-seq data
#' for x-axis sections. The rows represent genes. The first column
#' contains gene IDs and the second and subsequent columns contain
#' gene expression levels in sections.
#' @param y A data.frame object containing a simulated Tomo-seq data for y-axis
#' sections. The rows represent genes. The first column contains gene IDs and
#' the second and subsequent columns contain gene expression levels in sections.
#' @param z A data.frame object containing a simulated Tomo-seq data for
#' z-axis sections. The rows represent genes. The first column
#' contains gene IDs and the second and subsequent columns contain
#' gene expression levels in sections.
#' @param mask A 3D array that represents if each boxel is included to sample.
#' You can make a mask using `masker`.
#' @param query Vector of gene IDs
#' @param numIter How many times iterate
#' @param normCount Specifies the method to normalize
#' the expression amount data.
#' @param normMask Whether to normalize by mask or not
#' @return tomoSeq object
#' @importFrom purrr list_along
#' @examples
#' data("testx", "testy", "testz", "mask")
#' Estimate3dExpressions(
#'     testx,
#'     testy,
#'     testz,
#'     mask = mask,
#'     query = c("gene1"), 
#'     normCount = "countSum",
#'     normMask = TRUE
#' )
#' @export
Estimate3dExpressions <- function (
    x,
    y,
    z,
    mask,
    query,
    numIter = 100,
    normCount="countSum",
    normMask=TRUE
) {
    recList <- list_along(query)
    names(recList) <- query
    for (geneID in query) {
        recList[[geneID]] <- SingleEstimate(
            x,
            y,
            z,
            mask=mask,
            geneID=geneID,
            normCount=normCount,
            normMask=normMask,
            numIter=numIter
        )
    }
    retList <- list(
        "mask" = mask,
        "results" = recList
    )
    class(retList) <- "tomoSeq"
    return(retList)
}


#' Plot the trend of the value of the loss function.
#' @param tomoObj tomoSeq object
#' @param geneID single gene ID (string)
#' @importFrom dplyr %>%
#' @return NA
#' @examples
#' data(tomoObj)
#' PlotLossFunction(tomoObj, "gene2")
#' @export
PlotLossFunction <- function (tomoObj, geneID) {
    CheckParameters(tomoObj, geneID)
    tomoObj[["results"]][[geneID]][["errFunc"]] %>%
        plot(
            type = "l",
            main = geneID,
            xlab = "Iteration number",
            ylab = "Loss"
        )
}

#' Animate 2D expressions along one axis and generate GIF file.
#' @param tomoObj tomoSeq object
#' @param geneID single gene ID (string)
#' @param target "expression", "mask" or "unite" (combination of expression and
#' mask). Default is `expression`.
#' @param xaxis Number to specify as x-axis (1, 2 or 3). Default is `1`.
#' @param yaxis Number to specify as y-axis (1, 2 or 3). Default is `2`.
#' @param main A string used for the title of the plot. Default is `geneID`.
#' @param xlab Label of x axis. Default is `xaxis`.
#' @param ylab Label of y axis. Default is `yaxis`.
#' @param file Path of GIF file.
#' @param zlim Limit of value of heatmap. If target="mask", it is ignored.
#' @param interval interval of GIF animation.
#' @param aspectRatio A 2D vector that represents the ratio of figure. You can
#' specify the ratio as `c(width, height)`. If you don't specify the value of
#' this parameter, the ratio is calculated based on the number of sections
#' along each axis.
#' @importFrom stringr str_c
#' @importFrom dplyr %>%
#' @importFrom animation saveGIF
#' @return It generate GIF file.
#' @examples
#' if(interactive()) {
#'     data(tomoObj)
#'     Animate2d(tomoObj, "gene2", target = "expression", file = "example.gif")
#' }
#' @export
Animate2d <- function (
    tomoObj,
    geneID,
    target="expression",
    xaxis=1,
    yaxis=2,
    main=geneID,
    xlab=xaxis,
    ylab=yaxis,
    file=str_c(geneID, "_", target, "_", xaxis, "_", yaxis, ".gif"),
    zlim=NA,
    interval=0.1,
    aspectRatio=c()
) {
    CheckParameters(tomoObj, geneID)
    if (length(aspectRatio) != 0 & length(aspectRatio) != 2) {
        stop("`aspectRatio` should be a 2D vector.")
    }
    if (target == "mask" & is.na(zlim[1]) == FALSE) {
        warning('If target = "mask", parameter "zlim" is ignored.')
    }

    reconstArray <- tomoObj[["results"]][[geneID]][["reconst"]] %>%
        aperm(perm=c(xaxis, yaxis, 6 - (xaxis + yaxis)))
    maskArray <- tomoObj[["mask"]] %>%
        aperm(perm=c(xaxis, yaxis, 6 - (xaxis + yaxis)))

    saveGIF(
        AnimateForGIF(
            reconstArray = reconstArray,
            maskArray = maskArray,
            main = main,
            xlab = xlab,
            ylab = ylab,
            zlim = zlim,
            aspectRatio = aspectRatio,
            type = target
        ),
        movie.name=file,
        interval=interval,
        autobrowse=FALSE
    )
}

#' Plot expression of single gene along an axis
#' @param tomoObj tomoSeq object
#' @param geneID single gene ID (string)
#' @param axes axis ("x", "y" or "z")
#' @importFrom graphics par
#' @return NA
#' @examples
#' data(tomoObj)
#' Plot1dExpression(tomoObj, "gene2", "x")
#' @export
Plot1dExpression <- function (tomoObj, geneID, axes) {
    CheckParameters(tomoObj, geneID)
    convertList <- list("x" = 1, "y" = 2, "z" = 3)
    oldpar <- par(no.readonly=TRUE)
    on.exit(par(oldpar))
    marginal <- tomoObj[["results"]][[geneID]][[axes]]
    plot(marginal, type="l", lty=3, axes=FALSE, ann=FALSE)
    par(new=TRUE)
    plot(
        apply(
            tomoObj[["results"]][[geneID]][["reconst"]],
            convertList[[axes]],
            sum
        ),
        type="l",
        lty=2,
        ylim=range(marginal),
        col="red",
        xlab = "Section",
        ylab = "Expression level"
    )
}

#' Plot expressions of all genes along an axis
#' @param tomoSeqData A data.frame object containing a Tomo-seq data
#' for any sections. The rows represent genes. The first column
#' contains gene IDs and the second and subsequent columns contain
#' gene expression levels in sections.
#' @param ... Arguments which are related to plot parameters.
#' Prease refer to \code{\link[graphics]{plot}}.
#' @importFrom dplyr %>%
#' @return NA
#' @examples
#' data("testx")
#' Plot1dAllExpression(testx)
#' @export
Plot1dAllExpression <- function (tomoSeqData, ...) {
    tomoSeqData[, -1] %>% colSums() %>% plot(type="l", ...)
}

#' Convert reconstructed matrix to data.frame.
#' @param tomoObj tomoSeq object
#' @param geneID single gene ID
#' @importFrom dplyr %>%
#' @return Reconstruction result converted to dataframe.
#' @examples
#' data(tomoObj)
#' ToDataFrame(tomoObj, "gene2")
#' @export
ToDataFrame <- function (tomoObj, geneID) {
    CheckParameters(tomoObj, geneID)
    tomoObj[["results"]][[geneID]][["reconst"]] %>%
        MatrixToDataFrame() %>%
        return()
}

#' Get reconstructed matrix
#' @param tomoObj tomoSeq object
#' @param geneID single gene ID
#' @return Reconstruction result as matrix
#' @examples
#' data(tomoObj)
#' GetReconstructedResult(tomoObj, "gene2")
#' @export
GetReconstructedResult <- function (tomoObj, geneID) {
    CheckParameters(tomoObj, geneID)
    return(tomoObj[["results"]][[geneID]][["reconst"]])
}

#' Find peak genes on axial
#' @param tomoSeqData tomo-seq data of any axis
#' @param genes If run for all genes
#' @return A data frame consisting of gene ID, max of expression levels of the
#' gene, mean of expression levels calculated by excluding the maximum value and
#' section number showing the maximum expression level (0 means that there is
#' no such section).
#' @importFrom dplyr select
#' @examples
#' data(testx)
#' FindAxialGenes(testx)
#' @export
FindAxialGenes <- function (tomoSeqData, genes = "all") {
    if (length(genes) == 1 && genes == "all") {
        tomoSeqData %>%
            FindAxialGenesInner() %>%
            return()
    } else {
        tomoSeqData[as.vector(t(tomoSeqData[, 1]) %in% genes), ] %>%
            FindAxialGenesInner() %>%
            return()
    }
}

#' Download  part of the Tomo-seq data published by Junker et al.
#' @param verbose If you want to force downloads with or without cache,
#' set this to TRUE.
#' @return BiocFileCache object.
#' @export 
DownloadJunker2014 <- function ( verbose = FALSE ) {
    sheldAVURL <- "https://figshare.com/ndownloader/files/34384922"
    sheldVDURL <- "https://figshare.com/ndownloader/files/34384928"
    sheldLRURL <- "https://figshare.com/ndownloader/files/34384925"
    maskURL <- "https://figshare.com/ndownloader/files/34523069"

    bfc <- GetTomoseqrCache()
    DownloadData(bfc=bfc, rname="sheld_AV", URL=sheldAVURL)
    DownloadData(bfc=bfc, rname="sheld_VD", URL=sheldVDURL)
    DownloadData(bfc=bfc, rname="sheld_LR", URL=sheldLRURL)
    DownloadData(bfc=bfc, rname="mask97x97x97", URL=maskURL)
    return(bfc)
}

#' Load data of Junker2014 from cache.
#' @param tomoseqrCache Cache of tomoseqr. You can get it using
#' `DownloadJunker2014`.
#' @return List of tomo-seq data in cache.
#' @importFrom BiocFileCache bfcquery
#' @importFrom readr read_tsv
#' @export
LoadJunker2014 <- function (tomoseqrCache) {
    ridAV <- bfcquery(tomoseqrCache, "sheld_AV", "rname")$rid
    ridVD <- bfcquery(tomoseqrCache, "sheld_VD", "rname")$rid
    ridLR <- bfcquery(tomoseqrCache, "sheld_LR", "rname")$rid
    ridMask <- bfcquery(tomoseqrCache, "mask97x97x97", "rname")$rid
    sheldAV <- read_tsv(tomoseqrCache[[ridAV]])
    sheldVD <- read_tsv(tomoseqrCache[[ridVD]])
    sheldLR <- read_tsv(tomoseqrCache[[ridLR]])
    mask <- get(load(tomoseqrCache[[ridMask]]))
    return(
        list(
            sheldAV=sheldAV,
            sheldVD=sheldVD,
            sheldLR=sheldLR,
            mask=mask
        )
    )
}