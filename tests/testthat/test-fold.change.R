test_that("basic fold.change works", {
    a <- c(1, 2, 3)
    b <- c(2, 1, 3)
    c <- c(2, -2, 1)
    expect_equal(fold.change(a, b), matrix(c, length(c), 1))
})

test_that("fold.change avoids dividing by zero", {
    a <- c(0, 2, 3)
    b <- c(2, 0, 3)
    c <- c(10000, -10000, 1)
    expect_equal(fold.change(a, b), matrix(c, length(c), 1))
})

test_that("fold.change allows overriding BIG", {
    a <- c(0, 2, 3)
    b <- c(2, 0, 3)
    c <- c(100, -100, 1)
    expect_equal(fold.change(a, b, 100), matrix(c, length(c), 1))
})

test_that("fold.change works with vectors", {
    a <- c(1, 2, 3)
    b <- matrix(c(1, 2, 3, 3, 6, 21), nrow = length(a))
    c <- matrix(c(2, 2, 4), nrow = length(a))
    expect_equal(fold.change(a, b), c)
})

test_that("fold.change works with scalars", {
    a <- 2
    b <- matrix(c(1), nrow = length(a))
    c <- matrix(c(-2), nrow = length(a))
    expect_equal(fold.change(a, b), c)
})

test_that("fold.change fails when the number of rows differ", {
    a <- c(1, 2, 3)
    b <- c(1, 2, 3, 4)
    expect_error(fold.change(a, b))
})
