test_that("basic bb.test works", {
    x <- c(1, 2, 3)
    tx <- c(4, 5, 6)
    group <- c("group1", "group1", "group2")
    expect_equal(bb.test(x, tx, group)[["p.value"]], 0.51927813, tolerance = 1e-3)
})

test_that("bb.test alternative hypotheses work", {
    x <- c(1, 2, 3)
    tx <- c(4, 5, 6)
    group <- c("group1", "group1", "group2")
    # default should be two.sided
    expect_equal(bb.test(x, tx, group), bb.test(x, tx, group, "two.sided"), tolerance = 1e-3)
    # two.sided p-value should be twice 'greater' (since group2 > group1 in this example)
    expect_equal(bb.test(x, tx, group, "greater")[["p.value"]] * 2,
        bb.test(x, tx, group, "two.sided")[["p.value"]],
        tolerance = 1e-3
    )
    # less and greater should sum to one
    expect_equal(bb.test(x, tx, group, "greater")[["p.value"]],
        1 - bb.test(x, tx, group, "less")[["p.value"]],
        tolerance = 1e-3
    )
})

test_that("bb.test param 'tx' must be positive", {
    x <- c(1, 2, 3)
    tx <- c(0, 5, 6)
    group <- c("group1", "group1", "group2")
    expect_error(bb.test(x, tx, group), "tx must be positive.\n")
})

test_that("bb.test param 'x' must be non-negative", {
    x <- c(-1, 2, 3)
    tx <- c(4, 5, 6)
    group <- c("group1", "group1", "group2")
    expect_error(bb.test(x, tx, group), "x must be non-negative.\n")
})

test_that("bb.test param 'group' must be length ncol(x)", {
    x <- c(1, 2, 3)
    tx <- c(4, 5, 6)
    group <- c("group1", "group2")
    expect_error(
        bb.test(x, tx, group),
        "length of 'group' must be equal to the number of columns of 'x'"
    )
})

test_that("bb.test param 'tx' must be length ncol(x)", {
    x <- c(1, 2, 3)
    tx <- c(4, 6)
    group <- c("group1", "group1", "group2")
    expect_error(
        bb.test(x, tx, group),
        "The number of columns of 'tx' must be equal to the number of columns of 'x'"
    )
})

test_that("bb.test is testing at least two groups", {
    x <- c(1, 2, 3)
    tx <- c(4, 5, 6)
    group <- c("group1", "group1", "group1")
    expect_error(bb.test(x, tx, group), "Please compare at least two groups !!!")
})

test_that("bb.test must test exactly two groups for one-sided tests", {
    x <- c(1, 2, 3)
    tx <- c(4, 5, 6)
    group <- c("group1", "group2", "group3")
    expect_error(bb.test(x, tx, group, "less"), "One-sided test for two groups only !!!")
    expect_error(bb.test(x, tx, group, "greater"), "One-sided test for two groups only !!!")
})

test_that("bb.test succeeds with empty and full counts", {
    x <- c(25, 34, 26, 3, 14, 15)
    tx <- c(25, 34, 26, 18, 32, 38)
    group <- c(rep("na", 3), rep("group1", 3))
    expect_equal(bb.test(x, tx, group)[["p.value"]], 5.62228e-05, tolerance = 1e-3)
    expect_equal(bb.test(tx - x, tx, group)[["p.value"]], 5.62228e-05, tolerance = 1e-3)
})

test_that("bb.test succeeds when p = 0.5 exactly", {
    x <- rep(1, 10)
    tx <- c(7265, 6871, 8452, 8715, 8370, 7955, 7870, 7925, 7885, 8040)
    group <- c(rep("PRE", 5), rep("POST", 5))
    expect_equal(bb.test(x, tx, group, alternative = "less")[["p.value"]], 0.5, tolerance = 1e-3)
})
