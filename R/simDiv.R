#' @title Simulate diversification on an adaptive landscape
#'  
#' @description Simulate diversification on an adaptive landscape
#' 
#' @details The adaptiveness of the landscape is determined by the functions \code{la} and \code{mu}
#' which capture how speciation and extinction change as a function of position in the landscape
#' la, mu, r, d, x0, opt0, bsd, ngen
#' @param la function of position in the landscape giving birth rate at that point
#' @param mu function of position in the landscape giving death rate at that point
#' @param r the rate of evolutionary movement in the landscape
#' @param d the critical distance past which a new species is generated
#' @param x0 initial trait value of the initializing population
#' @param opt0 initial optimal value in the adaptive landscape
#' @param bsd SD of the Brownian drift of the adaptive optimum
#' @param ngen number of generations over which to run the simulation
#' 
#' @return A phylogeny with tip labels and edge length are derived from the input matrix
#' 
# @examples
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
# @seealso 
#' @export

# r <- 10
# d <- 0.1
# ngen <- 200
# bsd <- 0.05
# opt <- 0.5
# 
# la <- function(x) {
#     o <- rep(NA, length(x))
#     o[!is.na(x)] <- 1
#     
#     return(o)
# }
# 
# mu <- function(x, opt) {
#     1 * (dnorm(opt, opt, 0.2) - dnorm(x, opt, 0.2)) + 0
# }

simDiv <- function(la, mu, r, d, x0, opt0, bsd, ngen) {
    x <- rep(NA, 1e+05)
    s <- rep(NA, length(x))
    
    x[1] <- x0
    s[1] <- 1
    
    sim <- matrix(NA, nrow = ngen, ncol = 5)
    colnames(sim) <- c('x', 's', 'orig', 'ext', 'parent')
    sim[1, ] <- c(x[1], s[1], 1, NA, 0)
    
    xmat <- matrix(NA, nrow = length(s), ncol = ngen)
    xmat[, 1] <- x
    
    smat <- matrix(NA, nrow = length(x), ncol = ngen)
    smat[, 1] <- s
    
    j <- 2
    scounter <- 2
    
    for(i in 2:ngen) {
        if(all(is.na(x))) break
        
        opt <- rnorm(1, opt, bsd)
        
        allLa <- la(x)
        allMu <- mu(x, opt)
        
        birth <- sample(c(TRUE, FALSE), 1, prob = c(sum(allLa, na.rm = TRUE), sum(allMu, na.rm = TRUE)))
        
        if(sum(!is.na(x)) == 1) {
            thisOne <- which(!is.na(x))
        } else if(birth) {
            thisOne <- sample(which(!is.na(x)), size = 1, prob = allLa[!is.na(x)])
        } else {
            thisOne <- sample(which(!is.na(x)), size = 1, prob = allMu[!is.na(x)])
        }
        
        if(birth) {
            wiggle <- sample(c(-1, 1), 1) * rexp(1, r)
            
            x[j] <- x[thisOne] + wiggle
            if(abs(x[j] - x[thisOne]) > d) {
                s[j] <- scounter
                scounter <- scounter + 1
            } else {
                s[j] <- s[thisOne]
            }
            
            sim[j, ] <- c(x[j], s[j], i, NA, thisOne)
            
            j <- j + 1
        } else {
            x[thisOne] <- NA
            s[thisOne] <- NA
            
            sim[thisOne, 'ext'] <- i
        }
        
        # xmat[, i] <- x
        # smat[, i] <- s
    }
    
    sim <- as.data.frame(sim[!is.na(sim[, 1]), ])
    sim$parent_s <- c(0, sim$s[sim$parent])
    sim$ext[is.na(sim$ext)] <- ngen
    
    sim4tre <- plyr::ddply(sim, 's', function(x) {
        c(p = min(x$parent_s), o = min(x$orig), e = max(x$ext))
    })
    
    tre <- phyloFromSP(as.matrix(sim4tre))
    return(tre)
}
