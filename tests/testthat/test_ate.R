library("temporalEF")

context("Testing features of ate()")

tp <- seq_len(20)

test_that("ate() returns correct object", {
    tefs <- ate(tp)
    expect_is(tefs, "ate")
})
