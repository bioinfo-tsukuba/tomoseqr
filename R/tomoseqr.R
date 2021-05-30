#' @importFrom R6 R6Class
#' @importFrom dplyr %>%
#' @importFrom animation saveGIF
tomo_seq <- R6Class(
  classname = "tomoSeq",
  public = list(
    initialize = function (x, y, z, mask_shape="rectangle", species="") {
      # x, y, z: 各軸のtomo-seqデータ（全遺伝子分）
      # mask_shape: maskの形状の指定（"rectangle", "round", "halfround")
      # species: 種名（属名頭文字 + 種小名、いずれも小文字）
      MASK_SHAPE <- list(private$makeRectangle, private$makeRound, private$makeHalfRound)
      names(MASK_SHAPE) <- c("rectangle", "round", "halfround")

      private$x <- x
      private$y <- y
      private$z <- z
      private$gene_list <- private$extractGeneList()

      # 各single_geneオブジェクトを作成・辞書化
      for (gene in private$gene_list) {
        single_gene_object <- private$single_gene$new(x, y, z, gene)
        private$objects_each_gene <- private$objects_each_gene %>% append(single_gene_object)
      }
      names(private$objects_each_gene) <- private$gene_list

      # maskの作成
      # 1列目の遺伝子IDを除外するから、長さはそれぞれ -1
      private$mask <- MASK_SHAPE[[mask_shape]](length(x) - 1, length(y) - 1, length(z) - 1)

      if (species != "") {
        private$species <- species
        cat("Fetching from ensembl.....\n")
        data_name <- paste(species, "_gene_ensembl", sep="")
        data_base <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL")
        dataset <- biomaRt::useDataset(data_name, mart=data_base)
        private$gene_ID_name <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id",
                                               value = private$gene_list, mart = dataset)
        private$has_species_name <- TRUE
      }
    },

    # 各種getter ----------------------
    getGeneList = function () {
      return(private$gene_list)
    },

    getGene = function () {
      return(private$objects_each_gene)
    },

    getRaw = function () {
      return(list(private$x, private$y, private$z))
    },

    getMask = function () {
      return(private$mask)
    },

    getSpecies = function () {
      return(private$species)
    },

    getIDandName = function() {
      return(private$gene_ID_name)
    },
    # ---------------------------------

    estimate3dExpressions = function (queries=c()) {
      for (gene_ID in queries) {
        private$objects_each_gene[[gene_ID]]$estimate3dExpression(private$x, private$y, private$z, private$mask)
      }
    },

    plotLossFunction = function(gene_ID) {
      private$objects_each_gene[[gene_ID]]$plotLossFunction()
    },

    toDataFrame = function(gene_ID) {
      private$objects_each_gene[[gene_ID]]$toDataFrame()
    },

    animate2d = function (gene_ID, target="expression") {
      if (target == "expression") {
        private$objects_each_gene[[gene_ID]]$animate2dExpression()
      } else if (target == "mask") {
        private$objects_each_gene[[gene_ID]]$animate2dMask()
      } else if (target == "unite") {
        private$objects_each_gene[[gene_ID]]$animate2dUnite()
      } else {
        cat("ERROR: animate2d\n")
        cat(paste("Invalid option: target =", target, "\n"))
      }
    },

    plot1dExpression = function (gene_ID, axes) {
      private$objects_each_gene[[gene_ID]]$plot1dExpression(axes)
    },

    plot1dAllExpression = function (axes) {
      if (axes==1) {
        private$x[, -1] %>% colSums() %>% plot(type="l")
      } else if (axes==2) {
        private$y[, -1] %>% colSums() %>% plot(type="l")
      } else if (axes==3) {
        private$z[, -1] %>% colSums() %>% plot(type="l")
      } else {
        cat("axes must be 1, 2 or 3.\n")
      }
    }
  ),
  private = list(
    x = data.frame(),
    y = data.frame(),
    z = data.frame(),
    mask = array(0, dim=c(1,1,1)),
    gene_list = c(),
    objects_each_gene = c(),
    gene_ID_name = c(),
    has_species_name = FALSE,
    species = "",


    extractGeneList = function () {
      x_gene <- private$x[, 1] %>% t()
      y_gene <- private$y[, 1] %>% t()
      z_gene <- private$z[, 1] %>% t()
      x_and_y <- intersect(x_gene, y_gene)
      return(intersect(x_and_y, z_gene))
    },

    # 各種mask作成-------------------------------------------------------
    makeRectangle = function (x_len, y_len, z_len) {
      return(array(1, dim=c(x_len, y_len, z_len)))
    },

    makeRound = function (x_len, y_len, z_len) {
      mask <- array(0, dim=c(x_len, y_len, z_len))
      r <- mean(x_len, y_len, z_len) / 2
      for (x in 1:x_len){
        for (y in 1:y_len) {
          for (z in 1:z_len) {
            d <- (x - r)^2 + (y - r)^2 + (z - r)^2
            if (d >= (r * 0.75)^2 && d <= (r * 0.9)^2) {
              mask[x, y, z] <- 1
            }
          }
        }
      }
      return(mask)
    },

    makeHalfRound = function (x_len, y_len, z_len) {
      mask <- array(0, dim=c(x_len, y_len, z_len))
      r <- mean(x_len, y_len, z_len) / 2
      for (x in 1:x_len){
        for (y in 1:y_len) {
          for (z in 1:z_len) {
            d <- (x - r)^2 + (y - r)^2 + (z - r)^2
            if (d >= (r * 0.75)^2 && d <= (r * 0.9)^2 && x <= x_len / 2) {
              mask[x, y, z] <- 1
            }
          }
        }
      }
      return(mask)
    },
    # ---------------------------------------------------------------------

    # 各遺伝子についての推定結果等を保持するオブジェクト
    single_gene = R6Class(
      classname = "singleGene",
      public = list(
        gene_ID = "gene ID",
        X = matrix(0, nrow = 1, ncol = 1),
        Y = matrix(0, nrow = 1, ncol = 1),
        Z = matrix(0, nrow = 1, ncol = 1),
        reconst = array(0, dim = c(1, 1, 1)),
        mask = array(0, dim = c(1, 1, 1)),
        loss = c(),
        marginal_dist = array(0, dim=c(1,1,1)),
        already_reconstructed = FALSE,


        initialize = function (x, y, z, gene_ID) {
          self$gene_ID <- gene_ID
        },

        # estimate3dExpression()内で用いる関数
        getGeneExpression = function (tomoseq_data, gene_ID) {
          retval_matrix <- tomoseq_data[tomoseq_data[, 1] == gene_ID,]
          retval_matrix <- retval_matrix[, -1] %>% as.matrix()
          return(retval_matrix)
        },

        repMat = function (target_vector, n_times_repeat) {
          len_x <- length(target_vector) * n_times_repeat[1]
          len_y <- n_times_repeat[2]
          len_z <- n_times_repeat[3]
          rep_2d <- matrix(target_vector, nrow=len_x, ncol=len_y, byrow = F)
          rep_3d <- array(NA, dim=c(len_x, len_y, len_z))
          for (i in 1:len_z) {
            rep_3d[, , i] <- rep_2d
          }
          return(rep_3d)
        },

        # 3次元発現パターンを推定
        estimate3dExpression = function(X, Y, Z, mask, num_iter = 100) {
          self$mask <- mask
          sum_x <- X[, -1] %>% colSums()
          sum_y <- Y[, -1] %>% colSums()
          sum_z <- Z[, -1] %>% colSums()

          x0 <- X %>% self$getGeneExpression(self$gene_ID)
          y0 <- Y %>% self$getGeneExpression(self$gene_ID)
          z0 <- Z %>% self$getGeneExpression(self$gene_ID)
          x_len <- length(x0)
          y_len <- length(y0)
          z_len <- length(z0)

          mask_xproj <- self$mask %>% apply(1, sum)
          mask_yproj <- self$mask %>% apply(2, sum)
          mask_zproj <- self$mask %>% apply(3, sum)

          x <- x0 / sum_x * mask_xproj
          y <- y0 / sum_y * mask_yproj
          z <- z0 / sum_z * mask_zproj
          x[1, is.nan(x)] <- 0
          y[1, is.nan(y)] <- 0
          z[1, is.nan(z)] <- 0
          x[is.infinite(x)] <- max(x[x < Inf])
          y[is.infinite(y)] <- max(z[y < Inf])
          z[is.infinite(z)] <- max(y[z < Inf])

          x_raw <- sum_x
          y_raw <- sum_y
          z_raw <- sum_z

          m <- mean(c(sum(x_raw), sum(y_raw), sum(z_raw)))

          x <- x / sum(x) * m
          y <- y / sum(y) * m
          z <- z / sum(z) * m
          a <- self$mask

          er <- c()

          for (i in 1:num_iter) {
            xa <- a %>% apply(1, sum)
            a <- a * self$repMat(x / xa, c(1, dim(self$mask)[2], dim(self$mask)[3]))
            a[is.nan(a)] <- 0
            ya <- a %>% apply(2, sum)
            a <- a * aperm(self$repMat(y / ya, c(1, dim(self$mask)[1], dim(self$mask)[3])), perm = c(2, 1, 3))
            a[is.nan(a)] <- 0
            za <- a %>% apply(3, sum)
            a <- a * aperm(self$repMat(z / za, c(1, dim(self$mask)[1], dim(self$mask)[2])), perm = c(2, 3, 1))
            a[is.nan(a)] <- 0
            er <- append(er, sum((xa - x)^2) + sum((ya - y)^2) + sum((za - z)^2))
          }
          self$reconst <- a
          self$loss <- er
          self$marginal_dist <- list(x[1, ], y[1, ], z[1, ])
          self$already_reconstructed <- TRUE
        },

        plotLossFunction = function () {
          plot(self$loss, type="l", main=self$gene_ID, xlab="Iteration number", ylab="Loss")
        },

        animate2dExpression = function(axes1=1, axes2=2, main=self$gene_ID, xlab=axes1, ylab=axes2,
                                       file=paste(self$gene_ID, "_", axes1, "_", axes2, ".gif", sep=""),
                                       zlim=range(self$reconst), interval=0.1)
        {
          # if (self$already_reconstructed == TRUE) {
          #   if (private$has_species_name == TRUE) {
          #     gene_name <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id",
          #                       values = private$query, mart = private$gene_list)
          #     gene_name <- gene_name$external_gene_name
          #   } else {
              gene_name <- ""
          #   }
            self$animate2d(self$reconst, axes1=axes1, axes2=axes2,
                           main=paste(main, "(", gene_name, ")", sep=""), xlab=xlab, ylab=ylab, file=file,
                           zlim=zlim, interval=interval)
          # } else {
          #   cat("Before animate, please run estimateExpression().")
          #   return(1)
          # }
        },

        animate2dMask = function(axes1=1, axes2=2, main=self$gene_ID, xlab=axes1, ylab=axes2,
                                 file=paste(self$gene_ID, "_MASK_", axes1, "_", axes2, ".gif", sep=""),
                                 zlim=c(0, 1), interval=0.1)
        {
          self$animate2d(self$mask, axes1=axes1, axes2=axes2,
                         main=main, xlab=xlab, ylab=ylab, file=file,
                         zlim=zlim, interval=interval)
        },

        animate2dUnite = function(axes1=1, axes2=2, main=self$gene_ID, xlab=axes1, ylab=axes2,
                                 file=paste(self$gene_ID, "_unite_", axes1, "_", axes2, ".gif", sep=""),
                                 zlim=range(self$reconst), interval=0.1)
        {
          self$animateMaskAndExpression(axes1=axes1, axes2=axes2,
                         main=main, xlab=xlab, ylab=ylab, file=file,
                         zlim=zlim, interval=interval)
        },

        contourForAnimate = function (array_3d, main, xlab, ylab, zlim) {
          cat("generating")
          for (i in seq_along(array_3d[1, 1, ])) {
            cat("...")
            filled.contour(array_3d[, , i], main=paste(main, "_", i, sep=""), xlab=xlab,
                           ylab=ylab, zlim=zlim, asp=1, frame.plot=F)
          }
          cat("\n")
        },

        contourMaskAndExpression = function (mask_apermed, reconst_apermed, main, xlab, ylab, zlim) {
          label_list <- seq(zlim[1], floor(zlim[2]), length=6) %>% round()
          position_list <- label_list / zlim[2]
          cat("generating")
          for (i in seq_along(mask_apermed[1, 1, ])) {
            cat("...")
            # image.plot(reconst_apermed[, , i], zlim=zlim, breaks=seq(zlim[1], zlim[2], length=floor(zlim[2])), col=hcl.colors(floor(zlim[2])-1, palette="Oslo"))
            collist <- hcl.colors(floor(zlim[2])-1, palette="Oslo")
            ColorRamp<-colorRampPalette(collist)(100)
            ColorLevels<-seq(from=zlim[1], to=zlim[2], length=10000)
            par(mar=c(2,3,2,2), oma=c(0,0,0,0))
            layout(matrix(seq(2),nrow=2,ncol=1),widths=c(1),heights=c(3,0.5))
            image(reconst_apermed[, , i], zlim=zlim, xlab=xlab, ylab=ylab, breaks=seq(zlim[1], zlim[2], length=floor(zlim[2])), col=hcl.colors(floor(zlim[2])-1, palette="Oslo"), asp=1, axes=F)
            axis(1, seq(0, 1.0, by=0.2), seq(0, 1, by=0.2))
            axis(2, seq(0, 1.0, by=0.2), seq(0, 1, by=0.2), pos=0)
            mtext(xlab, side = 1, line = 2)
            mtext(ylab, side = 2, line = 1)
            par(new=T)
            image(mask_apermed[,,i], col=c("#000000", "#FFFFFF00"), main=paste(main, "_", i, seq=""), xlab=xlab, ylab=ylab, asp=1, axes=F)
            image(as.matrix(ColorLevels),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
            axis(1, position_list, label_list)
          }
          cat("\n")
        },

        animate2d = function (array3d, axes1, axes2, main, xlab, ylab, file, zlim, interval) {
          array3d_apermed <- aperm(array3d, perm=c(axes1, axes2, 6 - (axes1 + axes2)))
          saveGIF(self$contourForAnimate(array_3d=array3d_apermed, main=main, xlab=xlab, ylab=ylab, zlim=zlim),
                  movie.name=file, interval=interval)
        },

        animateMaskAndExpression = function (axes1, axes2, main, xlab, ylab, file, zlim, interval) {
          mask_apermed <- aperm(self$mask, perm=c(axes1, axes2, 6 - (axes1 + axes2)))
          reconst_apermed <- aperm(self$reconst, perm=c(axes1, axes2, 6 - (axes1 + axes2)))
          saveGIF(self$contourMaskAndExpression(mask_apermed = mask_apermed, reconst_apermed = reconst_apermed, main=main, xlab=xlab, ylab=ylab, zlim=zlim),
                  movie.name=file, interval=interval)
        },

        plot1dExpression = function (axes) {
          oldpar <- par(no.readonly=T)
          marginal <- self$marginal_dist[[axes]]
          plot(marginal, type="l", lty=3, axes=F, ann=F)
          par(new=T)
          plot(apply(self$reconst, axes, sum), type="l", lty=2, ylim=range(marginal), col="red")
          par(oldpar)
        },

        toDataFrame = function () {
          if (self$already_reconstructed == TRUE) {
            vec_reconst <- as.vector(self$reconst)
            dim <- dim(self$reconst)
            xlen <- dim[1]
            ylen <- dim[2]
            zlen <- dim[3]
            x_index <- rep(1:xlen, ylen * zlen)
            y_index <- 1:ylen %>%
                       sapply(function(p){rep(p, xlen)}) %>%
                       rep(zlen)
            z_index <- 1:zlen %>%
              sapply(function(p){rep(p, xlen * ylen)}) %>%
              as.vector()
            data.frame(x=x_index, y=y_index, z=z_index, value=vec_reconst) %>% return()
          } else {
             cat("have not reconstructed yet.")
          }
        }
      )
    )
  )
)
