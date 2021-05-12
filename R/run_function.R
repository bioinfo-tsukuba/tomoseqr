#' make set of tomoseq objects
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


#' plots the trend of the value of the loss function.
#' @param tomo_obj tomo_seq object
#' @param gene_ID single gene ID (string)
#' @export
plotLossFunction <- function(tomo_obj, gene_ID) {
  tomo_obj$plotLossFunction(gene_ID = gene_ID)
  return(invisible(tomo_obj))
}

#' animate 2D expressions along axes
#' @param tomo_obj tomo_seq object
#' @param gene_ID single gene ID (string)
#' @param  target "expression" or "mask"
#' @export 
animate2d <- function (tomo_obj, gene_ID, target="expression") {
  tomo_obj$animate2d(gene_ID, target)
  return(invisible(tomo_obj))
}

#' plot expression of single gene along axes
#' @param tomo_obj tomo_seq object
#' @param gene_ID single gene ID (string)
#' @param axes axes (1, 2 or 3)
#' @export
plot1dExpression <- function (tomo_obj, gene_ID, axes) {
  tomo_obj$plot1dExpression(gene_ID, axes)
  return(invisible(tomo_obj))
}

#' plot expressions of all genes along axes
#' @param tomo_obj tomo_seq object
#' @param axes axes (1, 2 or 3)
#' @export
plot1dAllExpression <- function(tomo_obj, axes) {
  tomo_obj$plot1dAllExpression(axes)
  return(invisible(tomo_obj))
}

#' convert reconstructed matrix to data.frame
#' @param tomo_obj tomo_seq object
#' @param gene_ID single gene ID
#' @export
toDataFrame <- function(tomo_obj, gene_ID) {
  tomo_obj$toDataFrame(gene_ID)
}