test_that("Reconstruction of 3D expression pattern", {
    load("simulateData.rda")
    load("rec_results.rda")
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
    expect_identical(sum((result1 - gene1_rec)^2), 0)
    expect_identical(sum((result2 - gene2_rec)^2), 0)
    expect_identical(sum((result3 - gene3_rec)^2), 0)
    expect_identical(sum((result4 - gene4_rec)^2), 0)
})
