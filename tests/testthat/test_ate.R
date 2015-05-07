library("temporalEF")

context("Testing features of ate()")

tp <- seq_len(20)

test_that("ate() returns correct object", {
    tefs <- ate(tp)
    expect_is(tefs, "ate")
})

test_that("ate handles inputs correctly; length", {
    expect_is(ate(1:2), "ate")
    expect_error(ate(1), "At least two time points are required for ATE.")
})
