
fold.change <- function(d1, d2, BIG = 1e4) {
    # fold change will be in range (-inf, -1) [1, inf)
    # and fold change for zeroed rows will be -BIG or BIG

    # convert input vectors and scalars to matrices
    if (is.vector(d1) || (is.atomic(d1) && length(d1) == 1)) {
        d1 <- matrix(d1, nrow = length(d1))
    }
    if (is.vector(d2) || (is.atomic(d2) && length(d2) == 1)) {
        d2 <- matrix(d2, nrow = length(d2))
    }

    # throw error if number of rows doesn't match
    if (nrow(d1) != nrow(d2)) {
        stop(sprintf("number of rows in d1 (%d) and d2 (%d) do not match.\n", nrow(d1), nrow(d2)))
    }

    # TODO: what if number of columns doesn't match?

    # for each row, calculate fold change (in avg value across cols)
    val <- matrix(0, nrow(d1), 1)

    for (i in 1:nrow(d1)) {
        v1 <- sum(d1[i,]) / ncol(d1)
        v2 <- sum(d2[i,]) / ncol(d2)

        if (v1 == 0) {
            if (v2 == 0) { # no data in either, fold.change = 1
                val[i] <- 1
            } else { # no data in v1, avoid dividing by zero
                val[i] <- BIG
            }

        } else {
            if (v2 == 0) { # no data in v2, avoid dividing by zero
                val[i] <- -BIG

            } else { # calculate fold change
                if (v1 > v2) {
                    val[i] <- -v1/v2
                } else {
                    if (v1 < v2) {
                        val[i] <- v2/v1
                    } else { # v2 == v1
                        val[i] <- 1
                    }
                }
            }
        }
    }

    return(val)
}

