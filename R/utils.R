.ncores <- function() {
    out <- .C("bbCores", ncores = as.integer(0))

    return(out$ncores)
}


.set.no.thread <- function(n.threads, max.thread) {
    out <- n.threads

    if (out != 1) {
        ncores <- .ncores()

        if (out > ncores) {
            out <- ncores
        }

        if (out < 1) {
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
