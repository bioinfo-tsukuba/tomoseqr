test_that("Reconstruction of 3D expression pattern", {
    load("trueResults.rda")
    data("testx", "testy", "testz", "mask")
    tomoSeqObject <- Estimate3dExpressions(
        testx,
        testy,
        testz,
        mask = mask,
        query = c("gene1", "gene2", "gene3", "gene4")
    )
    result1 <- GetReconstructedResult(tomoSeqObject, "gene1")
    result2 <- GetReconstructedResult(tomoSeqObject, "gene2")
    result3 <- GetReconstructedResult(tomoSeqObject, "gene3")
    result4 <- GetReconstructedResult(tomoSeqObject, "gene4")

    true1 <- GetReconstructedResult(trueResults, "gene1")
    true2 <- GetReconstructedResult(trueResults, "gene2")
    true3 <- GetReconstructedResult(trueResults, "gene3")
    true4 <- GetReconstructedResult(trueResults, "gene4")
    expect_identical(sum((result1 - true1)^2), 0)
    expect_identical(sum((result2 - true2)^2), 0)
    expect_identical(sum((result3 - true3)^2), 0)
    expect_identical(sum((result4 - true4)^2), 0)
})
