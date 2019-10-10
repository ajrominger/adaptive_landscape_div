#' @title Create phylogeny from a 'species-parent' matrix with origin and extinction times
#'  
#' @description Create a phylogeny from a matrix giving species and ancestry information.
#' 
#' @details An 'sp' matrix is a matrix where the first column is the species ID, the 
#' second column is the parent species, the third is the origin time of the species indicated by
#' the first column, and the forth column is that species extinction time
#' 
#' @param sp a 'species-parent' matrix as described in the details
#' 
#' @return A phylogeny with tip labels and edge length are derived from the input matrix
#' 
# @examples
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
# @seealso 
#' @export

phyloFromSP <- function(sp) {
    dep <- max(sp[, 4]) - sp[, 3]
    trim <- max(sp[, 4]) - sp[, 4]
    sp <- sp[, 1:2]
    
    # identify nodes in sp and branch lengths in oe
    spNew <- matrix(NA, ncol = 2, nrow = 0)
    el <- numeric(0)
    
    nodeCount <- -1
    for(i in 1:nrow(sp)) {
        el[spNew[, 1] == sp[i, 2]] <- el[spNew[, 1] == sp[i, 2]] - dep[i]
        spNew[spNew[, 1] == sp[i, 2], 1] <- nodeCount
        
        el <- c(el, rep(dep[i], 2))
        spNew <- rbind(spNew, cbind(sp[i, ], nodeCount))
        
        nodeCount <- nodeCount - 1
    }
    
    # browser()
    
    el[spNew[, 1] > 0] <- el[spNew[, 1] > 0] - trim[spNew[spNew[, 1] > 0, 1]]
    
    sp <- spNew
    rownames(sp) <- colnames(sp) <- NULL
    
    # print(sp)
    
    # get tip labels
    tiplab <- 0:max(sp)
    
    # convert nodes from negative to positive
    sp[sp < 0] <- max(sp) - sp[sp < 0]
    sp <- sp + 1
    
    # make phylo object
    tre <- list(edge = sp[, 2:1], tip.label = tiplab, Nnode = max(tiplab), edge.length = el)
    class(tre) <- 'phylo'
    
    return(tre)
}
