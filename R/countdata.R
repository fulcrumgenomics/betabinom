#######################################################################
#
# Author: Thang V. Pham, t.pham@vumc.nl
# 2013-06-25
#
# Citation for the beta-binomial test:
#
# Pham TV, Piersma SR, Warmoes M, Jimenez CR (2010) On the beta binomial model for analysis of spectral count data in label-free tandem mass spectrometry-based proteomics. Bioinformatics, 26(3):363-369.
#
#
# Citation for the inverted beta-binomial test:
#
# Pham TV, Jimenez CR (2012) An accurate paired sample test for count data. Bioinformatics, 28(18):i596-i602.
#
#######################################################################

ibb.test <- function(x,
                     tx,
                     group,
                     alternative = c("two.sided", "less", "greater"),
                     n.threads = -1,
                     BIG = 1e4,
                     verbose = TRUE) {

    #check x & tx
    if (any(tx <= 0))  {
        stop("tx must be positive.\n")
    }

    if (any(x < 0))  {
        stop("x must be non-negative.\n")
    }

    #making matrices
    if (is.vector(x)) {
        x <- matrix(x, nrow = 1)
    }

    K <- nrow(x)

    if (is.vector(tx)) {
        tx <- matrix(rep(tx, K), byrow=TRUE, nrow = K)
    }

    if (nrow(tx) == 1 && K > 1) {
        tx <- matrix(rep(tx, K), byrow=TRUE, nrow = K)
    }

    # check group
    if (length(group) != ncol(x)) {
        stop("length of 'group' must be equal the number of columns of 'x'")
    }

    if (ncol(tx) != ncol(x)) {
        stop("The number of columns of 'tx' must be equal to the number of columns of 'x'")
    }


    gf <- factor(group, levels = unique(group))

    group.values <- levels(gf)
    if (length(group.values) != 2) {
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

    init_block <- 3 * q.size;                 # qz + qlw + sqrt2phi
    block <- 4 * N + N * q.size + 5 * q.size; # magic
    factor <- 1;                              # sizeof(TYPE) / sizeof(double)

    max.threads <- .set.no.thread(n.threads, nrow(x))

    mem <- double((max.threads * block + init_block) * factor)

    mem[1] <- lower.bound;

    alternative <- match.arg(alternative)

    if (alternative == "two.sided") {
        mem[2] <- 0;
    }
    else {
        if (alternative == "less") {
            mem[2] <- -1;
        }
        else {
            mem[2] <- 1;
        }
    }

    a <- double(K*N)
    b <- double(K*N)
    ta <- double(K*N)
    tb <- double(K*N)

    ind_beg <- 1;
    ind_end <- N;

    for (i in 1:K) {

        a[ind_beg:ind_end] <- x[i, ia]
        b[ind_beg:ind_end] <- x[i, ib]
        ta[ind_beg:ind_end] <- tx[i, ia]
        tb[ind_beg:ind_end] <- tx[i, ib]

        ind_beg <- ind_beg + N;
        ind_end <- ind_end + N;
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
              fc = double(K))

    p.value <- out$p.value
    fc <- out$fc

    fc[fc < 1.0] <- -1.0 / fc[fc < 1.0]

    if (nrow(x) > 1) {
        total_g1 <- if (length(ia) > 1) rowSums(x[, ia]) else x[, ia]
        total_g2 <- if (length(ib) > 1) rowSums(x[, ib]) else x[, ib]

        fc[total_g1 == 0] <- fc[total_g1 == 0]*0 + BIG
        fc[total_g2 == 0] <- fc[total_g2 == 0]*0 - BIG
        fc[total_g1 == 0 & total_g2 == 0] <- 1
    }

    return(list(p.value = p.value, fc = fc))
}


bb.test <- function(x,
                    tx,
                    group,
                    alternative = c("two.sided", "less", "greater"),
                    n.threads = -1,
                    verbose = TRUE) {

    #check x & tx
    if (any(tx <= 0))  {
        stop("tx must be positive.\n")
    }

    if (any(x < 0))  {
        stop("x must be non-negative.\n")
    }


    #making matrices

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
    g.list <- vector('list', M)


    for (i in 1:M) {
        igroup <- which(gf == group.values[i])
        g.size[i] <- length(igroup);
        if (i > 1) {
            g.ind[i] <- g.ind[i-1] + g.size[i-1]
        }
        g.list[[i]] <- igroup
    }


    ## prepare C call

    N <- ncol(x)
    a <- double(K*N)
    ta <- double(K*N)

    ind.beg <- 1

    for (k in 1:K) {

        for (m in 1:M) {

            ind.end <- ind.beg + g.size[m] - 1

            a[ind.beg:ind.end] <- x[k, g.list[[m]]]

            ta[ind.beg:ind.end] <- tx[k, g.list[[m]]]

            ind.beg <- ind.end + 1
        }
    }

    block <- 2*N + M; # magic
    factor <- 1;      # sizeof(TYPE) / sizeof(double)

    max.threads <- .set.no.thread(n.threads, nrow(x))

    mem <- double(max.threads * block * factor)

    mem[1] <- 1 # theta.equal = TRUE

    alternative <- match.arg(alternative)

    if (alternative == "two.sided") {
        mem[2] <- 0;
    }
    else {
        if (M != 2) {
            stop("One-sided test for two groups only !!!")
        }

        if (alternative == "less") {
            mem[2] <- -1;
        }
        else {
            mem[2] <- 1;
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
              p.value = double(K))

    p.value <- out$p.value

    return (list(p.value = p.value))

}

normalize <- function(d) {

    total <- apply(d, 2, sum)

    m <- mean(total)

    factor <- total / m

    dnorm <-  d / (matrix(1, nrow(d), 1) %*% factor)

    return(dnorm)

}


fold.change <- function(d1, d2, BIG = 1e4) {

    if (is.vector(d1)) {
        d1 <- matrix(d1, nrow = length(d1))
    }

    if (is.vector(d2)) {
        d2 <- matrix(d2, nrow = length(d2))
    }

    val <- matrix(0, nrow(d1), 1)

    for (i in 1:nrow(d1)) {

        v1 <- sum(d1[i,]) / ncol(d1)

        v2 <- sum(d2[i,]) / ncol(d2)

        if (v1 == 0) {

            if (v2 == 0) {
                val[i] <- 1
            } else {
                val[i] <- BIG
            }
        } else {
            if (v2==0) {
                val[i] <- -BIG
            } else {
                if (v1 > v2) {
                    val[i] <- -v1/v2
                } else {
                    if (v1 < v2) {
                        val[i] <- v2/v1
                    } else {
                        val[i] <- 1
                    }
                }
            }
        }
    }

    return(val)
}


.ncores <- function() {

    out <- .C("bbCores", ncores = as.integer(0));

    return(out$ncores)

}


.set.no.thread <- function(n.threads, max.thread) {

    out <- n.threads

    if (out != 1) {

        ncores <- .ncores()

        if (out > ncores) {
            out <- ncores
        }

        if (out < 1)  {
            out <- ncores + out
            if (out < 1) {
                out <- 1
            }
        }
    }

    if (out > max.thread) {
        out <- max.thread
    }

    return(out)

}
