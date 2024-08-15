test_that("basic normalization works", {
    expect_equal(
        normalize(matrix(c(1, 3, 3, 9, 12, 0, 2, 2), nrow = 2, ncol = 4)),
        matrix(c(2, 6, 2, 6, 8, 0, 4, 4), nrow = 2, ncol = 4)
    )
})

test_that("sums of normalized columns are equal", {
    d <- matrix(c(1, 3, 3, 9, 12, 0, 2, 2), nrow = 2, ncol = 4)
    col_sums <- colSums(normalize(d))
    expect_true(var(col_sums) < 1e-7) # check that variance is nearly zero
})

test_that("normalizing zeroed column fails", {
    expect_error(
        normalize(matrix(c(1, 3, 3, 9, 12, 0, 0, 0), nrow = 2, ncol = 4)),
        "totals must be positive.\n"
    )
})

test_that("normalizing zeroed matrix fails", {
    expect_error(normalize(matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)), "totals must be positive.\n")
})

test_that("normalizing negative values fails", {
    expect_error(
        normalize(matrix(c(-1, 2, 2, 2), nrow = 2, ncol = 2)),
        "values must be non-negative.\n"
    )
})
