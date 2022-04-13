test_that("Extract peak genes at axial", {
    load("results_FindAxialGenes.rda")
    data("testx")
    testAll <- findAxialGenes(testx)
    testPartial <- findAxialGenes(testx, genes = c("gene1", "gene3"))
    expect_identical(testAll, resultAll)
    expect_identical(testPartial, resultPartial)
})

