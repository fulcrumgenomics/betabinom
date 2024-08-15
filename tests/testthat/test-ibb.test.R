test_that("basic ibb.test works", {
    x <- c(1, 3)
    tx <- c(4, 6)
    group <- c("group1", "group2")
    result <- ibb.test(x, tx, group)
    expect_equal(result[["p.value"]], 0.5276871, tolerance = 1e-3)
    expect_equal(result[["fc"]], 1.999985, tolerance = 1e-3)
})

test_that("ibb.test alternative hypotheses work", {
    x <- c(1, 3)
    tx <- c(4, 6)
    group <- c("group1", "group2")
    # default should be two.sided
    expect_equal(ibb.test(x, tx, group), ibb.test(x, tx, group, "two.sided"), tolerance = 1e-3)
    # two.sided p-value should be twice 'greater' (since group2 > group1 in this example)
    expect_equal(ibb.test(x, tx, group, "greater")[["p.value"]] * 2,
        ibb.test(x, tx, group, "two.sided")[["p.value"]],
        tolerance = 1e-3
    )
    # less and greater should sum to one
    expect_equal(ibb.test(x, tx, group, "greater")[["p.value"]],
        1 - ibb.test(x, tx, group, "less")[["p.value"]],
        tolerance = 1e-3
    )
})

test_that("ibb.test param 'tx' must be positive", {
    x <- c(1, 3)
    tx <- c(0, 6)
    group <- c("group1", "group2")
    expect_error(ibb.test(x, tx, group), "tx must be positive.\n")
})

test_that("ibb.test param 'x' must be non-negative", {
    x <- c(-1, 3)
    tx <- c(4, 6)
    group <- c("group1", "group2")
    expect_error(ibb.test(x, tx, group), "x must be non-negative.\n")
})

test_that("ibb.test param 'group' must be length ncol(x)", {
    x <- c(1, 3)
    tx <- c(4, 6)
    group <- c("group1")
    expect_error(
        ibb.test(x, tx, group),
        "length of 'group' must be equal to the number of columns of 'x'"
    )
})

test_that("ibb.test param 'tx' must be length ncol(x)", {
    x <- c(1, 3)
    tx <- c(4, 5, 6)
    group <- c("group1", "group2")
    expect_error(
        ibb.test(x, tx, group),
        "The number of columns of 'tx' must be equal to the number of columns of 'x'"
    )
})

test_that("ibb.test must test exactly two groups", {
    x <- c(1, 2, 3, 4)
    tx <- c(4, 5, 6, 7)
    group1 <- c("group1", "group1", "group1", "group1")
    group3 <- c("group1", "group2", "group2", "group3")
    expect_error(
        ibb.test(x, tx, group1),
        "Paired sample test of only one group is not supported."
    )
    expect_error(
        ibb.test(x, tx, group3),
        "Paired sample test of more than two groups is not supported."
    )
})

test_that("ibb.test must test groups with equal numbers of samples", {
    x <- c(1, 2, 3)
    tx <- c(4, 5, 6)
    group <- c("group1", "group2", "group2")
    expect_error(ibb.test(x, tx, group), "The number of samples in each group is not equal.")
})
