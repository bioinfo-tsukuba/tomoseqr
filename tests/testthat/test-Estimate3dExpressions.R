test_that("Reconstruction of 3D expression pattern", {
    load("trueResults.rda")
    data("testx", "testy", "testz", "mask")
    tomoSeqObject <- estimate3dExpressions(
        testx,
        testy,
        testz,
        mask = mask,
        query = c("gene1", "gene2", "gene3", "gene4")
    )
    result1 <- getReconstructedResult(tomoSeqObject, "gene1")
    result2 <- getReconstructedResult(tomoSeqObject, "gene2")
    result3 <- getReconstructedResult(tomoSeqObject, "gene3")
    result4 <- getReconstructedResult(tomoSeqObject, "gene4")

    true1 <- getReconstructedResult(trueResults, "gene1")
    true2 <- getReconstructedResult(trueResults, "gene2")
    true3 <- getReconstructedResult(trueResults, "gene3")
    true4 <- getReconstructedResult(trueResults, "gene4")
    expect_equal(sum((result1 - true1)^2), 0, tolerance=1e-6)
    expect_equal(sum((result2 - true2)^2), 0, tolerance=1e-6)
    expect_equal(sum((result3 - true3)^2), 0, tolerance=1e-6)
    expect_equal(sum((result4 - true4)^2), 0, tolerance=1e-6)
})
