test_that("Extract peak genes at axial", {
    load("results_FindAxialGenes.rda")
    data("testx")
    testAll <- findAxialGenes(testx)
    testPartial <- findAxialGenes(testx, genes = c("gene1", "gene3"))
    expect_equal(testAll, resultAll, tolerance=1e-6)
    expect_equal(testPartial, resultPartial, tolerance=1e-6)
})

