
normalize <- function(d) {

    if (any(d < 0))  {
        stop("values must be non-negative.\n")
    }
    total <- apply(d, 2, sum)
    if (any(total <= 0))  {
        stop("totals must be positive.\n")
    }
    m <- mean(total)
    factor <- total / m
    dnorm <-  d / (matrix(1, nrow(d), 1) %*% factor)
    return(dnorm)
}
