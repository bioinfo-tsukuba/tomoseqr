rm(list = ls())
library(tidyverse)
generateTomoseqData <- function(
    xLen,
    yLen,
    zLen,
    func,
    tValue = 1,
    fValue = 0
) {
    mask <- array(0, dim = c(xLen, yLen, zLen))
    indexList <- expand.grid(1:xLen, 1:yLen, 1:zLen)
    for (i in seq_len(length(indexList[,1]))) {
        x <- indexList[i, 1]
        y <- indexList[i, 2]
        z <- indexList[i, 3]
        if (func(x, y, z) == TRUE) {
            mask[x, y, z] <- tValue
        } else {
            mask[x, y, z] <- fValue
        }
    }
    return(mask)
}

maskFunc <- function(x, y, z) {
    return(
        (x - 25)^2 + 5 * (y - 25)^2 + 10 * (z - 25)^2 <= 25^2
    )
}
func1 <- function(x, y, z) {
    return(
        10 * (x - 25)^2 + 50 * (y - 25)^2 + 100 * (z - 25)^2 <= 25^2
    )
}
func3 <- function(x, y, z) {
    return(
        (x - 12.5)^2 + (y - 25)^2 + (z - 25)^2 <= (25 / 16)^2
    )
}
exp1 <- generateTomoseqData(50, 50, 50, func1, tValue = 100, fValue = 0)
expTmp <- generateTomoseqData(50, 50, 50, maskFunc, tValue = 10)
exp2 <- exp1 + expTmp
exp2[exp2 == 110] <- 100
exp3 <- generateTomoseqData(50, 50, 50, func3, tValue = 100, fValue = 0)
mask <- generateTomoseqData(50, 50, 50, maskFunc)

x1 <- apply(exp1, 1, sum)
x2 <- apply(exp2, 1, sum)
x3 <- apply(exp3, 1, sum)
x4 <- 30000 - (x1 + x2 + x3)
y1 <- apply(exp1, 2, sum)
y2 <- apply(exp2, 2, sum)
y3 <- apply(exp3, 2, sum)
y4 <- 30000 - (y1 + y2 + y3)
z1 <- apply(exp1, 3, sum)
z2 <- apply(exp2, 3, sum)
z3 <- apply(exp3, 3, sum)
z4 <- 30000 - (z1 + z2 + z3)

testx <- rbind(x1, x2, x3, x4)
testy <- rbind(y1, y2, y3, y4)
testz <- rbind(z1, z2, z3, z4)
rownames(testx) <- c("gene1", "gene2", "gene3", "gene4")
rownames(testy) <- c("gene1", "gene2", "gene3", "gene4")
rownames(testz) <- c("gene1", "gene2", "gene3", "gene4")
testx <- testx %>% as_tibble(rownames = "geneID")
testy <- testy %>% as_tibble(rownames = "geneID")
testz <- testz %>% as_tibble(rownames = "geneID")

usethis::use_data(testx, testy, testz, mask, overwrite = TRUE)
