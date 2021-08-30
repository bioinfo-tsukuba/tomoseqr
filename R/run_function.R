#' Make set of tomoseq objects
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
#' @param maskShape shape of mask. You can choose from "rectangle", "round" or
#' "halfround". The default is "rectangle"
#' @export
MakeTomoObjSet <- function (x, y, z, maskShape="rectangle") {
    if (is.element(maskShape, c("rectangle", "round", "halfround")) == FALSE) {
        stop('maskShape should be "rectangle", "round" or "halfround".')
    }
        return(tomoSeq$new(x=x, y=y,z=z, maskShape=maskShape))
    }

#' @importFrom methods is
CheckParameters <- function(tomoObj, query) {
    if (is(tomoObj, "tomoSeq") == FALSE) {
        stop(paste("invalid class:", class(tomoObj),
                   "\nFirst argument must be tomoSeq class object."
             )
        )
    }
    if (is.element(query, tomoObj$geneList) == FALSE) {
        stop(paste('gene "', query, '" is not in data.', sep=''))
    }
}
#' Estimate 3d expressions
#' @param tomoObj tomoSeq object
#' @param query Vector of gene IDs
#' @return tomoSeq object
#' @export
#' @note  You can do the same things with
#' `tomoObj$Estimate3dExpressions(query)`.
Estimate3dExpressions <- function (tomoObj, query) {
    CheckParameters(tomoObj, query)
    tomoObj$Estimate3dExpressions(queries=query)
    return(invisible(tomoObj))
    }


#' Plot the trend of the value of the loss function.
#' @param tomoObj tomoSeq object
#' @param geneID single gene ID (string)
#' @export
#' @note  You can do the same things with
#' `tomoObj$PlotLossFunction(geneID)`.
PlotLossFunction <- function (tomoObj, geneID) {
    CheckParameters(tomoObj, geneID)
    tomoObj$PlotLossFunction(geneID=geneID)
    return(invisible(tomoObj))
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
#' @export
#' @note  You can do the same thing with `tomoObj$Animate2d(geneID, ...)`.
Animate2d <- function (tomoObj,
                       geneID,
                       target="expression",
                       xaxis=1, yaxis=2,
                       main=geneID,
                       xlab=xaxis, ylab=yaxis,
                       file=paste(geneID, "_", target, "_",
                                  xaxis, "_", yaxis, ".gif",
                                  sep=""
                       ),
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
    tomoObj$Animate2d(geneID=geneID,
                       target=target,
                       xaxis=xaxis, yaxis=yaxis,
                       main=main,
                       xlab=xlab, ylab=ylab,
                       file=file,
                       zlim=zlim,
                       interval=interval,
                       aspectRatio=aspectRatio
    )
    return(invisible(tomoObj))
    }

#' Plot expression of single gene along an axis
#' @param tomoObj tomoSeq object
#' @param geneID single gene ID (string)
#' @param axes axis (1, 2 or 3)
#' @export
#' @note  You can do the same thing with
#' `tomoObj$Plot1dExpression(geneID, axes)`.
Plot1dExpression <- function (tomoObj, geneID, axes) {
    CheckParameters(tomoObj, geneID)
    tomoObj$Plot1dExpression(geneID, axes)
    return(invisible(tomoObj))
    }

#' Plot expressions of all genes along an axis
#' @param tomoObj tomoSeq object
#' @param axes axis (1, 2 or 3)
#' @export
#' @note  You can do the same thing with
Plot1dAllExpression <- function (tomoObj, axes) {
    tomoObj$Plot1dAllExpression(axes)
    return(invisible(tomoObj))
    }

#' Convert reconstructed matrix to data.frame.
#' @param tomoObj tomoSeq object
#' @param geneID single gene ID
#' @export
#' @note  You can do the same thing with
#' `tomoObj$ToDataFrame(geneID)`.
ToDataFrame <- function (tomoObj, geneID) {
    CheckParameters(tomoObj, geneID)
    tomoObj$ToDataFrame(geneID)
    }

#' Get reconstructed matrix
#' @param tomoObj tomoSeq object
#' @param geneID single gene ID
#' @export
#' @note  You can do the same thing with
#' `tomoObj$GetReconstructedResult(geneID)`.
GetReconstructedResult <- function (tomoObj, geneID) {
    CheckParameters(tomoObj, geneID)
    tomoObj$GetReconstructedResult(geneID)
    }
