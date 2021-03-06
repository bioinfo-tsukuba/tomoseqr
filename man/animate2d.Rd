% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_function.R
\name{animate2d}
\alias{animate2d}
\title{Animate 2D expressions along one axis and generate GIF file.}
\usage{
animate2d(
  tomoObj,
  geneID,
  target = "expression",
  xaxis = 1,
  yaxis = 2,
  main = geneID,
  xlab = xaxis,
  ylab = yaxis,
  file = str_c(geneID, "_", target, "_", xaxis, "_", yaxis, ".gif"),
  zlim = NA,
  interval = 0.1,
  aspectRatio = c()
)
}
\arguments{
\item{tomoObj}{tomoSeq object}

\item{geneID}{single gene ID (string)}

\item{target}{"expression", "mask" or "unite" (combination of expression and
mask). Default is \code{expression}.}

\item{xaxis}{Number to specify as x-axis (1, 2 or 3). Default is \code{1}.}

\item{yaxis}{Number to specify as y-axis (1, 2 or 3). Default is \code{2}.}

\item{main}{A string used for the title of the plot. Default is \code{geneID}.}

\item{xlab}{Label of x axis. Default is \code{xaxis}.}

\item{ylab}{Label of y axis. Default is \code{yaxis}.}

\item{file}{Path of GIF file.}

\item{zlim}{Limit of value of heatmap. If target="mask", it is ignored.}

\item{interval}{interval of GIF animation.}

\item{aspectRatio}{A 2D vector that represents the ratio of figure. You can
specify the ratio as \code{c(width, height)}. If you don't specify the value of
this parameter, the ratio is calculated based on the number of sections
along each axis.}
}
\value{
It generate GIF file.
}
\description{
Animate 2D expressions along one axis and generate GIF file.
}
\examples{
if(interactive()) {
    data(tomoObj)
    animate2d(tomoObj, "gene2", target = "expression", file = "example.gif")
}
}
