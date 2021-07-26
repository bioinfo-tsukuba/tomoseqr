# tomoseqr

`tomoseqr` is an R package for analyzing Tomo-seq (a method to obtain genome-wide expression data with spatial resolution) data.

## Usage

### Installation

```{r}
library(devtools)
install_github("bioinfo-tsukuba/tomoseqr")
```

### Data preparation

Please prepare Tomo-seq data that meets the following requirements.

1. It is a `data.frame` object for each axis.
1. Its **first cloumn** has gene ID. It's not enough that only row
names indicate gene ID.
1. The order of the second and subsequent columns should be the same as
the order of the sections.
1. It has a header.

#### Data example

```{r}
  gene_ID   section1 section2   section3   section4
1   gene1   0.000000   0.0000 2867.75420 9086.81135
2   gene2 440.599448 531.7915   36.91591  484.06813
3   gene3  75.446821 833.9432  736.82367  559.89157
4   gene4 506.865166 930.0414  880.26654   52.85974
5   gene5   2.159842 271.6788  210.06446  445.08979
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
