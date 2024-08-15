#######################################################################
#
# Author: Thang V. Pham, t.pham@vumc.nl
# 2013-06-25
#
# Citation for the beta-binomial test:
#
# Pham TV, Piersma SR, Warmoes M, Jimenez CR (2010) On the beta binomial model for analysis of
# spectral count data in label-free tandem mass spectrometry-based proteomics.
# Bioinformatics, 26(3):363-369.
#
#######################################################################

bb.test <- function(x,
                    tx,
                    group,
                    alternative = c("two.sided", "less", "greater"),
                    n.threads = -1,
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

    # check size
    if (length(group) != ncol(x)) {
        stop("length of 'group' must be equal to the number of columns of 'x'")
    }

    if (ncol(tx) != ncol(x)) {
        stop("The number of columns of 'tx' must be equal to the number of columns of 'x'")
    }

    gf <- factor(group)

    group.values <- levels(gf)

    M <- length(group.values)

    if (length(group.values) < 2) {
        stop("Please compare at least two groups !!!")
    }

    g.size <- double(M)
    g.ind <- double(M)
    g.list <- vector("list", M)


    for (i in 1:M) {
        igroup <- which(gf == group.values[i])
        g.size[i] <- length(igroup)
        if (i > 1) {
            g.ind[i] <- g.ind[i - 1] + g.size[i - 1]
        }
        g.list[[i]] <- igroup
    }


    ## prepare C call

    N <- ncol(x)
    a <- double(K * N)
    ta <- double(K * N)

    ind.beg <- 1

    for (k in 1:K) {
        for (m in 1:M) {
            ind.end <- ind.beg + g.size[m] - 1
            a[ind.beg:ind.end] <- x[k, g.list[[m]]]
            ta[ind.beg:ind.end] <- tx[k, g.list[[m]]]
            ind.beg <- ind.end + 1
        }
    }

    block <- 2 * N + M # magic
    factor <- 1 # sizeof(TYPE) / sizeof(double)

    max.threads <- .set.no.thread(n.threads, nrow(x))

    mem <- double(max.threads * block * factor)

    mem[1] <- 1 # theta.equal = TRUE

    alternative <- match.arg(alternative)

    if (alternative == "two.sided") {
        mem[2] <- 0
    } else {
        if (M != 2) {
            stop("One-sided test for two groups only !!!")
        }

        if (alternative == "less") {
            mem[2] <- -1
        } else {
            mem[2] <- 1
        }
    }

    if (!verbose) {
        max.threads <- -max.threads
    }

    out <- .C("bb",
        as.integer(K),
        as.double(a),
        as.double(ta),
        as.integer(M),
        as.integer(g.size),
        as.integer(g.ind),
        as.double(mem),
        as.integer(max.threads),
        p.value = double(K)
    )

    p.value <- out$p.value

    return(list(p.value = p.value))
}
