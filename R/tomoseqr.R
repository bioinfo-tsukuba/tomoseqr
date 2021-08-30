#' tomoSeq object
#' @description The class that contains a list of information about each genes.
#' @importFrom R6 R6Class
#' @importFrom dplyr %>%
#' @importFrom animation saveGIF
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
#' @param maskShape Shape of mask. You can choose the value from
#' `"rectangle"`, `"round"` or `"halfround"`. The default is `"rectangle"`.
        initialize = function (x, y, z, maskShape="rectangle") {
            ## x, y, z: Tomo-seq data (about all genes) of each axis.
            ## maskShape: The shape of mask
            ## ("rectangle", "round" or "halfround").
            maskShapeList <- list(private$MakeRectangle,
                               private$MakeRound,
                               private$MakeHalfRound
            )
            names(maskShapeList) <- c("rectangle", "round", "halfround")

            private$x <- x
            private$y <- y
            private$z <- z
            private$valGeneList <- private$ExtractGeneList()

            ## make each singleGene objects and compile as dictionary.
            for (gene in private$valGeneList) {
                singleGeneObject <- private$singleGene$new(x, y, z, gene)
                private$objectsEachGene <- private$objectsEachGene %>%
                    append(singleGeneObject)
            }
            names(private$objectsEachGene) <- private$valGeneList

            ## Create mask.
            ## Each length must be 1 shorter because first column (gene ID) is
            ## excluded from reconstruction.
            private$valMask <- maskShapeList[[maskShape]](
                                                         length(x) - 1,
                                                         length(y) - 1,
                                                         length(z) - 1
                                )
        },

#' @description Reconstructs 3D expression patterns of genes which are
#' specified. See also [Estimate3dExpressions()].
#' @param queries A vector consists of gene IDs.
        Estimate3dExpressions = function (queries=c()) {
            for (geneID in queries) {
                private$objectsEachGene[[geneID]]$Estimate3dExpression(
                    private$x, private$y, private$z, private$valMask
                )
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
        Animate2d = function (geneID,
                              target,
                              xaxis, yaxis,
                              main,
                              xlab, ylab,
                              file,
                              zlim,
                              interval,
                              aspectRatio=c()
                    ) {
            if (length(aspectRatio) != 0 & length(aspectRatio) != 2) {
                stop("`aspectRatio` should be a 2D vector.")
            }
            if (target == "expression") {
                private$objectsEachGene[[geneID]]$Animate2dExpression(
                    xaxis=xaxis, yaxis=yaxis,
                    main=main,
                    xlab=xlab, ylab=ylab,
                    zlim=zlim,
                    file=file,
                    interval=interval,
                    aspectRatio=aspectRatio
                )
            } else if (target == "mask") {
                private$objectsEachGene[[geneID]]$Animate2dMask(
                    xaxis=xaxis, yaxis=yaxis,
                    main=main,
                    xlab=xlab, ylab=ylab,
                    file=file,
                    interval=interval,
                    aspectRatio=aspectRatio
                )
            } else if (target == "unite") {
                private$objectsEachGene[[geneID]]$Animate2dUnite(
                    xaxis=xaxis, yaxis=yaxis,
                    main=main,
                    xlab=xlab, ylab=ylab,
                    zlim=zlim,
                    file=file,
                    interval=interval,
                    aspectRatio=aspectRatio
                )
            } else {
                stop(paste("Invalid option: target =", target, "\n"))
            }
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
        valMask = array(0, dim=c(1,1,1)),
        valGeneList = c(),
        objectsEachGene = c(),

        ExtractGeneList = function () {
            xGene <- private$x[, 1] %>% t()
            yGene <- private$y[, 1] %>% t()
            zGene <- private$z[, 1] %>% t()
            xAndy <- intersect(xGene, yGene)
            return(intersect(xAndy, zGene))
        },

        ## Mask-create functions -----------------------------------------------
        MakeRectangle = function (xLen, yLen, zLen) {
            return(array(1, dim=c(xLen, yLen, zLen)))
        },

        MakeRound = function (xLen, yLen, zLen) {
            mask <- array(0, dim=c(xLen, yLen, zLen))
            r <- mean(xLen, yLen, zLen) / 2
            for (x in 1:xLen) {
                for (y in 1:yLen) {
                    for (z in 1:zLen) {
                        d <- (x - r)^2 + (y - r)^2 + (z - r)^2
                        if (d >= (r * 0.75)^2 && d <= (r * 0.9)^2) {
                            mask[x, y, z] <- 1
                        }
                    }
                }
            }
            return(mask)
        },

        MakeHalfRound = function (xLen, yLen, zLen) {
            mask <- array(0, dim=c(xLen, yLen, zLen))
            r <- mean(xLen, yLen, zLen) / 2
            for (x in 1:xLen) {
                for (y in 1:yLen) {
                    for (z in 1:zLen) {
                        d <- (x - r)^2 + (y - r)^2 + (z - r)^2
                        if (d >= (r * 0.75)^2 &&
                            d <= (r * 0.9)^2 &&
                            x <= xLen / 2
                            ) {
                            mask[x, y, z] <- 1
                        }
                    }
                }
            }
            return(mask)
        },
        ## ---------------------------------------------------------------------

        ## This object has reconstruction result about a gene.
        singleGene = R6Class(
            classname = "singleGene",
            public = list(
                geneID = "gene ID",
                X = matrix(0, nrow = 1, ncol = 1),
                Y = matrix(0, nrow = 1, ncol = 1),
                Z = matrix(0, nrow = 1, ncol = 1),
                reconst = array(0, dim = c(1, 1, 1)),
                mask = array(0, dim = c(1, 1, 1)),
                loss = c(),
                marginalDist = array(0, dim=c(1,1,1)),
                alreadyReconstructed = FALSE,

                initialize = function (x, y, z, geneID) {
                    self$geneID <- geneID
                },

                ## This function is used in Estimate3dExpression().
                GetGeneExpression = function (tomoSeqData, geneID) {
                    retvalMatrix <- tomoSeqData[
                        tomoSeqData[, 1] == geneID,
                    ]
                    retvalMatrix <- retvalMatrix[, -1] %>% as.matrix()
                    return(retvalMatrix)
                },

                RepMat = function (targetVector, nTimesRepeat) {
                    lenX <- length(targetVector) * nTimesRepeat[1]
                    lenY <- nTimesRepeat[2]
                    lenZ <- nTimesRepeat[3]
                    rep2d <- matrix(targetVector, nrow=lenX,
                                     ncol=lenY, byrow = F
                              )
                    rep3d <- array(NA, dim=c(lenX, lenY, lenZ))
                    for (i in 1:lenZ) {
                        rep3d[, , i] <- rep2d
                    }
                    return(rep3d)
                },

                ## Reconstruct 3D expression pattern.
                Estimate3dExpression = function (X, Y, Z, mask, numIter = 100) {
                    self$mask <- mask
                    sumX <- X[, -1] %>% colSums()
                    sumY <- Y[, -1] %>% colSums()
                    sumZ <- Z[, -1] %>% colSums()

                    x0 <- X %>% self$GetGeneExpression(self$geneID)
                    y0 <- Y %>% self$GetGeneExpression(self$geneID)
                    z0 <- Z %>% self$GetGeneExpression(self$geneID)
                    xLen <- length(x0)
                    yLen <- length(y0)
                    zLen <- length(z0)

                    maskX <- self$mask %>% apply(1, sum)
                    maskY <- self$mask %>% apply(2, sum)
                    maskZ <- self$mask %>% apply(3, sum)

                    x <- x0 / sumX * maskX
                    y <- y0 / sumY * maskY
                    z <- z0 / sumZ * maskZ
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
                    a <- self$mask

                    er <- c()

                    for (i in 1:numIter) {
                        xa <- a %>% apply(1, sum)
                        a <- a * self$RepMat(x / xa, c(1,
                                                       dim(self$mask)[2],
                                                       dim(self$mask)[3]
                                                     )
                                 )
                        a[is.nan(a)] <- 0
                        ya <- a %>% apply(2, sum)
                        a <- a * aperm(self$RepMat(y / ya, c(1,
                                                             dim(self$mask)[1],
                                                             dim(self$mask)[3]
                                                           )
                                       ),
                                       perm = c(2, 1, 3)
                                 )
                        a[is.nan(a)] <- 0
                        za <- a %>% apply(3, sum)
                        a <- a * aperm(self$RepMat(z / za, c(1,
                                                             dim(self$mask)[1],
                                                             dim(self$mask)[2]
                                                           )
                                        ),
                                        perm = c(2, 3, 1))
                        a[is.nan(a)] <- 0
                        er <- append(er, sum((xa - x)^2) +
                                         sum((ya - y)^2) +
                                         sum((za - z)^2)
                              )
                    }
                    self$reconst <- a
                    self$loss <- er
                    self$marginalDist <- list(x[1, ], y[1, ], z[1, ])
                    self$alreadyReconstructed <- TRUE
                },

                CheckReconstructed = function () {
                    if (self$alreadyReconstructed == FALSE) {
                        stop(paste("No reconstructed results.",
                                   "You need to run Estimate3dExpressions()"
                             )
                        )
                    }
                },

                PlotLossFunction = function () {
                    self$CheckReconstructed()
                    plot(self$loss, type="l", main=self$geneID,
                         xlab="Iteration number", ylab="Loss"
                    )
                },

                Animate2dExpression = function (xaxis, yaxis,
                                                main,
                                                xlab, ylab,
                                                file,
                                                zlim,
                                                interval,
                                                aspectRatio
                                      ) {
                    self$CheckReconstructed()
                    if (is.na(zlim[1]) == TRUE) {
                        realZlim <- range(self$reconst)
                    } else {
                        realZlim <- zlim
                    }
                    self$Animate2d(self$reconst,
                                   xaxis=xaxis, yaxis=yaxis,
                                   main=main,
                                   xlab=xlab, ylab=ylab,
                                   file=file,
                                   zlim=realZlim,
                                   interval=interval,
                                   aspectRatio=aspectRatio
                    )
                },

                Animate2dMask = function (xaxis, yaxis,
                                          main,
                                          xlab, ylab,
                                          file,
                                          interval,
                                          aspectRatio
                                ) {
                    self$CheckReconstructed()
                    self$Animate2d(self$mask,
                                   xaxis=xaxis, yaxis=yaxis,
                                   main=main,
                                   xlab=xlab, ylab=ylab,
                                   file=file,
                                   zlim=c(0, 1),
                                   interval=interval,
                                   aspectRatio=aspectRatio
                    )
                },

                Animate2dUnite = function (xaxis, yaxis,
                                           main,
                                           xlab, ylab,
                                           file,
                                           zlim,
                                           interval,
                                           aspectRatio
                                 ) {
                    self$CheckReconstructed()
                    if (is.na(zlim[1]) == TRUE) {
                        realZlim <- range(self$reconst)
                    } else {
                        realZlim <- zlim
                    }
                    self$AnimateMaskAndExpression(xaxis=xaxis, yaxis=yaxis,
                                                  main=main,
                                                  xlab=xlab, ylab=ylab,
                                                  file=file,
                                                  zlim=realZlim,
                                                  interval=interval,
                                                  aspectRatio=aspectRatio
                    )
                },

                ContourForAnimate = function (array3d,
                                              main,
                                              xlab, ylab,
                                              zlim,
                                              aspectRatio
                                    ) {
                    if (length(aspectRatio) < 2) {
                        arrayDim <- dim(array3d)
                        asp <- arrayDim[2] / arrayDim[1]
                    } else {
                        asp <- aspectRatio[2] / aspectRatio[1]
                    }
                    message("generating", appendLF=FALSE)
                    for (i in seq_along(array3d[1, 1, ])) {
                        message(".", appendLF=FALSE)
                        filled.contour(array3d[, , i],
                                       main=paste(main, "_", i, sep=""),
                                       xlab=xlab,
                                       ylab=ylab, zlim=zlim, asp=asp,
                                       frame.plot=F
                        )
                    }
                    message("")
                },

                ContourMaskAndExpression = function (maskApermed,
                                                     reconstApermed,
                                                     main,
                                                     xlab, ylab,
                                                     zlim,
                                                     aspectRatio
                                           ) {
                    if (length(aspectRatio) < 2) {
                        ## Dim of reconstructed matrix should be equal to
                        ## that of mask.
                        plotDim <- dim(maskApermed)
                        asp <- plotDim[2] / plotDim[1]
                    } else {
                        asp <- aspectRatio[2] / aspectRatio[1]
                    }
                    labelList <- seq(zlim[1], floor(zlim[2]), length=6) %>%
                                      round()
                    position_list <- label_list / zlim[2]
                    message("generating", appendLF=FALSE)
                    collist <- hcl.colors(floor(zlim[2])-1, palette="Oslo")
                    ColorRamp<-colorRampPalette(collist)(100)
                    ColorLevels<-seq(from=zlim[1], to=zlim[2], length=100)
                    for (i in seq_along(mask_apermed[1, 1, ])) {
                        message(".", appendLF=FALSE)
                        par(mar=c(2,3,2,2), oma=c(0,0,0,0))
                        layout(matrix(seq(2), nrow=2, ncol=1), widths=c(1),
                               heights=c(3, 0.5)
                        )
                        image(reconstApermed[, , i], zlim=zlim, xlab=xlab,
                              ylab=ylab, breaks=seq(zlim[1], zlim[2],
                              length=floor(zlim[2])),
                              col=hcl.colors(floor(zlim[2])-1, palette="Oslo"),
                              asp=asp, axes=F
                        )
                        axis(1, seq(0, 1.0, by=0.2), seq(0, 1, by=0.2))
                        axis(2, seq(0, 1.0, by=0.2), seq(0, 1, by=0.2), pos=0)
                        mtext(xlab, side = 1, line = 2)
                        mtext(ylab, side = 2, line = 1)
                        par(new=T)
                        image(maskApermed[,,i], col=c("#000000", "#FFFFFF00"),
                              main=paste(main, "_", i, seq=""), xlab=xlab,
                              ylab=ylab, asp=asp, axes=F
                        )
                        image(as.matrix(ColorLevels),col=ColorRamp, xlab="",
                              ylab="", cex.axis=1, xaxt="n", yaxt="n"
                        )
                        axis(1, positionList, labelList)
                    }
                    message("")
                },

                Animate2d = function (array3d,
                                      xaxis, yaxis,
                                      main,
                                      xlab, ylab,
                                      file,
                                      zlim,
                                      interval,
                                      aspectRatio
                            ) {
                    array3dApermed <- aperm(array3d,
                                             perm=c(xaxis, yaxis,
                                                    6 - (xaxis + yaxis)
                                             )
                                       )
                    saveGIF(self$ContourForAnimate(array3d=array3dApermed,
                                                   main=main,
                                                   xlab=xlab,
                                                   ylab=ylab,
                                                   zlim=zlim,
                                                   aspectRatio=aspectRatio
                            ),
                            movie.name=file,
                            interval=interval,
                            autobrowse=FALSE
                    )
                },

                AnimateMaskAndExpression = function (xaxis, yaxis,
                                                     main,
                                                     xlab, ylab,
                                                     file,
                                                     zlim,
                                                     interval,
                                                     aspectRatio
                                           ) {
                    maskApermed <- aperm(self$mask,
                                          perm=c(xaxis, yaxis,
                                                 6 - (xaxis + yaxis)
                                          )
                                    )
                    reconstApermed <- aperm(self$reconst,
                                             perm=c(xaxis, yaxis,
                                                    6 - (xaxis + yaxis)
                                             )
                                       )
                    saveGIF(self$ContourMaskAndExpression(
                                maskApermed = maskApermed,
                                reconstApermed = reconstApermed,
                                main=main,
                                xlab=xlab, ylab=ylab,
                                zlim=zlim,
                                aspectRatio=aspectRatio
                            ),
                            movie.name=file,
                            interval=interval,
                            autobrowse=FALSE
                    )
                },

                Plot1dExpression = function (axes) {
                    self$CheckReconstructed()
                    oldpar <- par(no.readonly=T)
                    marginal <- self$marginalDist[[axes]]
                    plot(marginal, type="l", lty=3, axes=F, ann=F)
                    par(new=T)
                    plot(apply(self$reconst, axes, sum), type="l", lty=2,
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
                        data.frame(x=xIndex, y=yIndex, z=zIndex,
                                   value=vecReconst
                        ) %>% return()
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
