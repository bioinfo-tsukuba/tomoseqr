# tomoseqr

`tomoseqr` is an R package for analyzing Tomo-seq (a method to obtain genome-wide expression data with spatial resolution) data.

## Installation

```{r}
library(devtools)
install_github("bioinfo-tsukuba/tomoseqr")
```

## Example usage

```{r}
library(tomoseqr)
tomo_x <- tomoseqr::tomoseq_test_x
tomo_y <- tomoseqr::tomoseq_test_y
tomo_z <- tomoseqr::tomoseq_test_z
tomo_obj <- makeTomoObjSet(tomo_x, tomo_y, tomo_z)
estimate3dExpressions(tomo_obj, c("gene1", "gene2"))
animate2d(tomo_obj, "gene1")
```

![example](./inst/gene1_expression_1_2.gif)

```{r}
tomo_obj2 <- makeTomoObjSet(tomo_x, tomo_y, tomo_z, mask_shape = "round")
estimate3dExpressions(tomo_obj2, c("gene1", "gene2"))
animate2d(tomo_obj2, "gene1", target="unite")
```

![example](./inst/gene1_unite_1_2.gif)

## Contact

Ryosuke Matsuzawa / [shingenmochi](https://github.com/shingenmochi)
