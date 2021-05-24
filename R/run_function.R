#' Make set of tomoseq objects
#' @param x marginal data (x axes)
#' @param y marginal data (y axes)
#' @param z marginal data (z axes)
#' @param mask_shape shape of mask.
#' @param species scientific name
#' @export
makeSet <- function (x, y, z, mask_shape="rectangle", species="") {
  return(tomo_seq$new(x=x, y=y,z=z, mask_shape=mask_shape, species=species))
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

#' Animate 2D expressions along axes
#' @param tomo_obj tomo_seq object
#' @param gene_ID single gene ID (string)
#' @param  target "expression" or "mask"
#' @export 
#' @note  You can do the same things with
#' `tomo_obj$animate2d(gene_ID, target)`.
animate2d <- function (tomo_obj, gene_ID, target="expression") {
  tomo_obj$animate2d(gene_ID, target)
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