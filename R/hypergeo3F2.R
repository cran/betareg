h3f2 <- function(a, b, z, n, maxiter = 5000, eps = 1e-14) {
    res <- .Call("h3f2",
                 a = as.double(a),
                 b = as.double(b),
                 z = as.double(z),
                 n = as.integer(n),
                 maxiter = as.integer(maxiter),
                 eps = as.double(eps),
                 out = double(n))
    res$out
}

## hypergeo3F2 <- function(u1, u2, u3, l1, l2, z,
##                         n, maxiter = 5000, eps = 1e-14) {
##     res <- .C("hypergeo3f2",
##               u1 = as.double(u1),
##               u2 = as.double(u2),
##               u3 = as.double(u3),
##               l1 = as.double(l1),
##               l2 = as.double(l2),
##               z = as.double(z),
##               n = as.integer(n),
##               maxiter = as.integer(maxiter),
##               eps = as.double(eps),
##               out = double(n))
##     res$out
## }
