#' tomoSeq object
#' @description The class that contains a list of information about each genes.
#' @importFrom R6 R6Class
#' @importFrom dplyr %>%
#' @importFrom animation saveGIF
#' @importFrom stringr str_c
tomoSeq <- R6Class(
    classname = "tomoSeq",
    public = list(

#' @description Make tomoSeq object. This method is not available to users.
#' Instead of it, [MakeTomoObjSet()] is available.
#' @param x A data.frame object containing a simulated Tomo-seq data for
#' x-axis sections. The rows represent genes. The first column contains gene IDs
#' and the second and subsequent columns contain gene expression levels
#' in sections.
#' @param y A data.frame object containing a simulated Tomo-seq data for y-axis
#' sections. The rows represent genes. The first column contains gene IDs and
#' the second and subsequent columns contain gene expression levels in sections.
#' @param z A data.frame object containing a simulated Tomo-seq data for z-axis
#' sections. The rows represent genes. The first column contains gene IDs and
#' the second and subsequent columns contain gene expression levels in sections.
#' @param mask A 3D array that represents if each boxel is included to sample.
#' You can make a mask using `masker`.
        initialize = function (x, y, z, mask) {
            ## x, y, z: Tomo-seq data (about all genes) of each axis.
            private$x <- x
            private$y <- y
            private$z <- z
            private$valGeneList <- private$ExtractGeneList()
            private$valMask <- mask

            ## make each singleGene objects and compile as dictionary.
            for (gene in private$valGeneList) {
                singleGeneObject <- private$singleGene$new(x, y, z, gene)
                private$objectsEachGene <- private$objectsEachGene %>%
                    append(singleGeneObject)
            }
            names(private$objectsEachGene) <- private$valGeneList
        },

#' @description Reconstructs 3D expression patterns of genes which are
#' specified. See also [Estimate3dExpressions()].
#' @param queries A vector consists of gene IDs.
#' @param normCount Specifies the method to normalize
#' the expression amount data.
#' @param normMask Whether to normalize by mask or not
        Estimate3dExpressions = function (
            queries=c(),
            normCount="countSum",
            normMask=TRUE
        ) {
            ## Check normCount and normMask
            stopMsg <- str_c(
                'normCount must be either "countSum",',
                '"none", or a list of length 3.'
            )
            if (is.character(normCount)) {
                if (normCount != "countSum" & normCount != "none") {
                    stop(stopMsg)
                }
            } else if (is.list(normCount)) {
                if (length(normCount) != 3) {
                    stop(stopMsg)
                }
            } else {
                stop(stopMsg)
            }

            if (is.logical(normMask) == FALSE) {
                stop("normMask must be a boolean.")
            }

            for (geneID in queries) {
                estimatedResult <- Estimate3dExpression(
                    private$x,
                    private$y,
                    private$z,
                    geneID,
                    private$valMask,
                    normCount=normCount,
                    normMask=normMask,
                    numIter=100
                )
                private$objectsEachGene[[geneID]]$reconst <- estimatedResult[["reconst"]]
                private$objectsEachGene[[geneID]]$alreadyReconstructed <- TRUE
                private$objectsEachGene[[geneID]]$loss <- estimatedResult[["errFunc"]]
                private$objectsEachGene[[geneID]]$marginalDist <- estimatedResult[["marginalDist"]]
            }
        },

#' @description Plot the transition of the value of the loss function.
#' @param geneID A gene ID as string.
        PlotLossFunction = function (geneID) {
            private$objectsEachGene[[geneID]]$PlotLossFunction()
        },

#' @description Export a result of reconstruction as data.frame.
#' @param geneID A gene ID as string.
        ToDataFrame = function (geneID) {
            private$objectsEachGene[[geneID]]$ToDataFrame()
        },

#' @description Export a result of reconstruction as 3D matrix.
#' @param geneID A gene ID as string.
        GetReconstructedResult = function (geneID) {
            private$objectsEachGene[[geneID]]$GetReconstructedResult()
        },

#' @description Reconstruct a 3D expression pattern of a gene.
#' @param geneID A gene ID as string.
#' @param target A target of exportation. You can choose from "expression",
#' "mask" or "unite". The default value is "expression".
#' @param xaxis One of the axes of reconstructed matrix,
#' which is used the x-axis of the animation. You can choose from 1, 2 or 3.
#' The default value is 1.
#' @param yaxis One of the axes of reconstructed matrix,
#' which is used the y-axis of the animation. You can choose from 1, 2 or 3.
#' The default value is 2.
#' @param main The title of animation. It is NOT the file name.
#' The default value is same as geneID.
#' @param xlab A string, that is a label of x-axis. The default value is same
#' as `xaxis`.
#' @param ylab A string, that is a label of y-axis. The default value is same
#' as `yaxis`.
#' @param file A name of GIF file that is exported. The default value is
#' generated using geneID, xaxis and yaxis.
#' @param zlim Limits of expression levels that is displayed. You can specify
#' it as `c(min, max)`. The default value is automatically calculated using
#' the result of reconstruction.
#' @param interval An interval of GIF animation. The default value is 0.1.
#' @param aspectRatio A 2D vector of aspect ratio of animation.
#' You can specify the ratio as `c(width, height)`. If you don't specify
#' the value of this parameter, the ratio is calculated based on
#' the number of sections along each axis.
        Animate2d = function (
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
            aspectRatio=c()
        ) {
            if (length(aspectRatio) != 0 & length(aspectRatio) != 2) {
                stop("`aspectRatio` should be a 2D vector.")
            }
            reconstArray <- private$objectsEachGene[[geneID]]$reconst %>%
                aperm(perm=c(xaxis, yaxis, 6 - (xaxis + yaxis)))

            maskArray <- private$valMask %>%
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
        },

#' @description Plot expression of single gene along axis
#' @param geneID A gene ID as string
#' @param axes An axis of which you want to plot expression (1, 2 or 3).
        Plot1dExpression = function (geneID, axes) {
            private$objectsEachGene[[geneID]]$Plot1dExpression(axes)
        },

#' @description Plot expressions of all genes along an axis
#' @param axes An axis of which you want to plot expression (1, 2 or 3).
        Plot1dAllExpression = function (axes) {
            if (axes==1) {
                private$x[, -1] %>% colSums() %>% plot(type="l")
            } else if (axes==2) {
                private$y[, -1] %>% colSums() %>% plot(type="l")
            } else if (axes==3) {
                private$z[, -1] %>% colSums() %>% plot(type="l")
            } else {
                stop("axes must be 1, 2 or 3.\n")
            }
        }
    ),
    active = list(

        ## getter ----------------------

#' @field geneList A list of gene ID. It's read-only.
        geneList = function (value) {
            if (missing(value)) {
                return(private$valGeneList)
            } else {
                stop("geneList is read-only.")
            }
        },

#' @field expMatOriginal Original expression matrices along each axis.
#' It's read-only.
        expMatOriginal = function (value) {
            if (missing(value)) {
                return(list(private$x, private$y, private$z))
            } else {
                stop("Original matrices are read-only.")
            }
        },

#' @field mask 3D matrix, whose elements are 0 (not included in sample)
#' or 1 (included in sample). It's read-only.
        mask = function (value) {
            if (missing(value)) {
                return(private$valMask)
            } else {
                stop("mask is read-only.")
            }
        }
    ),

    private = list(
        x = data.frame(),
        y = data.frame(),
        z = data.frame(),
        valMask = array(0, dim=c(1, 1, 1)),
        valGeneList = c(),
        objectsEachGene = c(),

        ExtractGeneList = function () {
            xGene <- private$x[, 1] %>% t()
            yGene <- private$y[, 1] %>% t()
            zGene <- private$z[, 1] %>% t()
            xAndy <- intersect(xGene, yGene)
            return(intersect(xAndy, zGene))
        },

        ## This object has reconstruction result about a gene.
        singleGene = R6Class(
            classname = "singleGene",
            public = list(
                geneID = "gene ID",
                reconst = array(0, dim = c(1, 1, 1)),
                loss = c(),
                marginalDist = array(0, dim=c(1, 1, 1)),
                alreadyReconstructed = FALSE,

                initialize = function (x, y, z, geneID) {
                    self$geneID <- geneID
                },

                CheckReconstructed = function () {
                    if (self$alreadyReconstructed == FALSE) {
                        stop(
                            paste(
                                "No reconstructed results.",
                                "You need to run Estimate3dExpressions()"
                            )
                        )
                    }
                },

                PlotLossFunction = function () {
                    self$CheckReconstructed()
                    plot(
                        self$loss,
                        type="l",
                        main=self$geneID,
                        xlab="Iteration number",
                        ylab="Loss"
                    )
                },

                Plot1dExpression = function (axes) {
                    self$CheckReconstructed()
                    oldpar <- par(no.readonly=T)
                    marginal <- self$marginalDist[[axes]]
                    plot(marginal, type="l", lty=3, axes=F, ann=F)
                    par(new=T)
                    plot(
                        apply(self$reconst, axes, sum),
                        type="l",
                        lty=2,
                        ylim=range(marginal), col="red"
                    )
                    par(oldpar)
                },

                ToDataFrame = function () {
                    self$CheckReconstructed()
                    if (self$alreadyReconstructed == TRUE) {
                        vecReconst <- as.vector(self$reconst)
                        dim <- dim(self$reconst)
                        xlen <- dim[1]
                        ylen <- dim[2]
                        zlen <- dim[3]
                        xIndex <- rep(1:xlen, ylen * zlen)
                        yIndex <- 1:ylen %>%
                            sapply(function (p) {rep(p, xlen)}) %>%
                            rep(zlen)
                        zIndex <- 1:zlen %>%
                            sapply(function (p) {rep(p, xlen * ylen)}) %>%
                            as.vector()
                        data.frame(
                            x=xIndex,
                            y=yIndex,
                            z=zIndex,
                            value=vecReconst
                        ) %>%
                            return()
                    } else {
                        stop("No result of reconstruction.")
                    }
                },

                GetReconstructedResult = function () {
                    self$CheckReconstructed()
                    if (self$alreadyReconstructed == TRUE) {
                        return(self$reconst)
                    } else {
                        stop("No result of reconstruction.")
                    }
                }
            )
        )
    )
)
