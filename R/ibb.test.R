#######################################################################
#
# Author: Thang V. Pham, t.pham@vumc.nl
# 2013-06-25
#
# Citation for the inverted beta-binomial test:
#
# Pham TV, Jimenez CR (2012) An accurate paired sample test for count data.
# Bioinformatics, 28(18):i596-i602.
#
#######################################################################

ibb.test <- function(x,
                     tx,
                     group,
                     alternative = c("two.sided", "less", "greater"),
                     n.threads = -1,
                     BIG = 1e4,
                     verbose = TRUE) {
    # check x & tx
    if (any(tx <= 0)) {
        stop("tx must be positive.\n")
    }

    if (any(x < 0)) {
        stop("x must be non-negative.\n")
    }

    # making matrices
    if (is.vector(x)) {
        x <- matrix(x, nrow = 1)
    }

    K <- nrow(x)

    if (is.vector(tx)) {
        tx <- matrix(rep(tx, K), byrow = TRUE, nrow = K)
    }

    if (nrow(tx) == 1 && K > 1) {
        tx <- matrix(rep(tx, K), byrow = TRUE, nrow = K)
    }

    # check group
    if (length(group) != ncol(x)) {
        stop("length of 'group' must be equal to the number of columns of 'x'")
    }

    if (ncol(tx) != ncol(x)) {
        stop("The number of columns of 'tx' must be equal to the number of columns of 'x'")
    }


    gf <- factor(group, levels = unique(group))

    group.values <- levels(gf)
    if (length(group.values) == 1) {
        stop("Paired sample test of only one group is not supported.")
    } else if (length(group.values) != 2) {
        stop("Paired sample test of more than two groups is not supported.")
    }

    ia <- which(gf == group.values[1])
    ib <- which(gf == group.values[2])

    if (length(ia) != length(ib)) {
        stop("The number of samples in each group is not equal.")
    }


    ## prepare C call
    lower.bound <- 5.0
    N <- length(ia)
    q.size <- 2^15

    init_block <- 3 * q.size # qz + qlw + sqrt2phi
    block <- 4 * N + N * q.size + 5 * q.size # magic
    factor <- 1 # sizeof(TYPE) / sizeof(double)

    max.threads <- .set.no.thread(n.threads, nrow(x))

    mem <- double((max.threads * block + init_block) * factor)

    mem[1] <- lower.bound

    alternative <- match.arg(alternative)

    if (alternative == "two.sided") {
        mem[2] <- 0
    } else {
        if (alternative == "less") {
            mem[2] <- -1
        } else {
            mem[2] <- 1
        }
    }

    a <- double(K * N)
    b <- double(K * N)
    ta <- double(K * N)
    tb <- double(K * N)

    ind_beg <- 1
    ind_end <- N

    for (i in 1:K) {
        a[ind_beg:ind_end] <- x[i, ia]
        b[ind_beg:ind_end] <- x[i, ib]
        ta[ind_beg:ind_end] <- tx[i, ia]
        tb[ind_beg:ind_end] <- tx[i, ib]

        ind_beg <- ind_beg + N
        ind_end <- ind_end + N
    }

    if (!verbose) {
        max.threads <- -max.threads
    }

    out <- .C("ibb",
        as.integer(K),
        as.double(a),
        as.double(b),
        as.double(ta),
        as.double(tb),
        as.integer(N),
        as.double(mem),
        as.integer(max.threads),
        p.value = double(K),
        fc = double(K)
    )

    p.value <- out$p.value
    fc <- out$fc

    fc[fc < 1.0] <- -1.0 / fc[fc < 1.0]

    if (nrow(x) > 1) {
        total_g1 <- if (length(ia) > 1) rowSums(x[, ia]) else x[, ia]
        total_g2 <- if (length(ib) > 1) rowSums(x[, ib]) else x[, ib]

        fc[total_g1 == 0] <- fc[total_g1 == 0] * 0 + BIG
        fc[total_g2 == 0] <- fc[total_g2 == 0] * 0 - BIG
        fc[total_g1 == 0 & total_g2 == 0] <- 1
    }

    return(list(p.value = p.value, fc = fc))
}
