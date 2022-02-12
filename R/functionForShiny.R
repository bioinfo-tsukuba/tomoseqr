#' @importFrom stringr str_c
#' @importFrom dplyr %>%
#' @importFrom animation saveGIF
#' @importFrom shiny withProgress
ShinyAnimate2d <- function (
    tomoObj,
    geneID,
    target,
    xaxis,
    yaxis,
    main,
    xlab,
    ylab,
    file,
    zlim,
    interval,
    aspectRatio
) {
    reconstArray <- tomoObj[["results"]][[geneID]][["reconst"]] %>%
        aperm(perm=c(xaxis, yaxis, 6 - (xaxis + yaxis)))
    maskArray <- tomoObj[["mask"]] %>%
        aperm(perm=c(xaxis, yaxis, 6 - (xaxis + yaxis)))

    withProgress(message='generating GIF', value=0, {
        saveGIF(
            ShinyAnimateForGIF(
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
    })
}

#' @importFrom shiny
#' incProgress
ShinyAnimateForGIF <- function (
    reconstArray,
    maskArray,
    main,
    xlab,
    ylab,
    zlim,
    aspectRatio,
    type
) {
    asp <- aspectRatio[2] / aspectRatio[1]
    plotArray <- MakePlotArray(
        reconstArray=reconstArray,
        maskArray=maskArray,
        zlim=zlim,
        type=type
    )
    nTimesRepeat <- length(plotArray[1, 1, ])
    for (i in seq_along(plotArray[1, 1, ])) {
        BasePlot(
            plotArray,
            sectionNumber=i,
            main,
            xlab,
            ylab,
            asp
        )
        incProgress(
            1 / nTimesRepeat,
            detail = str_c("Plot:", i, " / ", nTimesRepeat)
        )
    }
    incProgress(1, detail = "Converting to GIF...")
}
