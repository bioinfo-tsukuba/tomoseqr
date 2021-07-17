#' Make set of tomoseq objects
#' @param x A data.frame object containing a simulated Tomo-seq data for x-axis sections.
#'          The rows represent genes. The first column contains gene IDs and the second and
#'          subsequent columns contain gene expression levels in sections.
#' @param y A data.frame object containing a simulated Tomo-seq data for y-axis sections.
#'          The rows represent genes. The first column contains gene IDs and the second and
#'          subsequent columns contain gene expression levels in sections.
#' @param z A data.frame object containing a simulated Tomo-seq data for z-axis sections.
#'          The rows represent genes. The first column contains gene IDs and the second and
#'          subsequent columns contain gene expression levels in sections.
#' @param mask_shape shape of mask.
#' @export
makeTomoObjSet <- function (x, y, z, mask_shape="rectangle") {
  return(tomo_seq$new(x=x, y=y,z=z, mask_shape=mask_shape))
}

#' Estimate 3d expression
#' @param tomo_obj tomo_seq object
#' @param query Vector of gene ID
#' @return tomo_seq object
#' @export
#' @note  You can do the same things with
#' `tomo_obj$estimate3dExpressions(query)`.
estimate3dExpressions <- function(tomo_obj, query) {
  tomo_obj$estimate3dExpressions(queries=query)
  return(invisible(tomo_obj))
}


#' Plot the trend of the value of the loss function.
#' @param tomo_obj tomo_seq object
#' @param gene_ID single gene ID (string)
#' @export
#' @note  You can do the same things with
#' `tomo_obj$plotLossFunction(gene_ID)`.
plotLossFunction <- function(tomo_obj, gene_ID) {
  tomo_obj$plotLossFunction(gene_ID = gene_ID)
  return(invisible(tomo_obj))
}

#' Animate 2D expressions along one axis and generate GIF file.
#' @param tomo_obj tomo_seq object
#' @param gene_ID single gene ID (string)
#' @param target "expression", "mask" or "unite" (combination of expression and mask). Default is `expression`.
#' @param axes1 Number to specify as x-axis (1, 2 or 3). Default is `1`.
#' @param axes2 Number to specify as y-axis (1, 2 or 3). Default is `2`.
#' @param main A string used for the title of the plot. Default is `gene_ID`.
#' @param xlab Label of x axis. Default is `axes1`.
#' @param ylab Label of y axis. Default is `axes2`.
#' @param file Path of GIF file.
#' @param zlim Limit of value of heatmap. If target="mask", it is ignored.
#' @param interval interval of GIF animation.
#' @export 
#' @note  You can do the same things with
#' `tomo_obj$animate2d(gene_ID, ...)`.
animate2d <- function (tomo_obj, gene_ID, target="expression", axes1=1, axes2=2, main=gene_ID, xlab=axes1, ylab=axes2,
                       file=paste(gene_ID, "_", target, "_", axes1, "_", axes2, ".gif", sep=""), zlim=NA, interval=0.1)
                       {
                         if (target == "mask" & is.na(zlim[1]) == FALSE) {
                           warning('If target = "mask", parameter "zlim" is ignored.')
                         }
                         tomo_obj$animate2d(gene_ID=gene_ID, target=target, axes1=axes1, axes2=axes2, main=main, xlab=xlab, ylab=ylab,
                                            file=file, zlim=zlim, interval=interval)
                                            return(invisible(tomo_obj))
}

#' Plot expression of single gene along axes
#' @param tomo_obj tomo_seq object
#' @param gene_ID single gene ID (string)
#' @param axes axes (1, 2 or 3)
#' @export
#' @note  You can do the same things with
#' `tomo_obj$plot1dExpression(gene_ID, axes)`.
plot1dExpression <- function (tomo_obj, gene_ID, axes) {
  tomo_obj$plot1dExpression(gene_ID, axes)
  return(invisible(tomo_obj))
}

#' Plot expressions of all genes along axes
#' @param tomo_obj tomo_seq object
#' @param axes axes (1, 2 or 3)
#' @export
#' @note  You can do the same things with
plot1dAllExpression <- function(tomo_obj, axes) {
  tomo_obj$plot1dAllExpression(axes)
  return(invisible(tomo_obj))
}

#' Convert reconstructed matrix to data.frame
#' @param tomo_obj tomo_seq object
#' @param gene_ID single gene ID
#' @export
#' @note  You can do the same things with
#' `tomo_obj$toDataFrame(gene_ID)`.
toDataFrame <- function(tomo_obj, gene_ID) {
  tomo_obj$toDataFrame(gene_ID)
}

#' Get reconstructed matrix
#' @param tomo_obj tomo_seq object
#' @param gene_ID single gene ID
#' @export
#' @note  You can do the same things with
#' `tomo_obj$getReconstructedResult(gene_ID)`.
getReconstructedResult <- function(tomo_obj, gene_ID) {
  tomo_obj$getReconstructedResult(gene_ID)
}