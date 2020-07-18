gs.heatmap <- function(gs.data, limit = NULL, cols = NULL, margins = c(6, 10),
    ...) {
    
    ds = dim(gs.data)
    if (length(ds) != 2) {
        print("gs.data needs to be a matrix-like object!")
        return(invisible(1))
    }
    else if (ds[1] < 2) {
        print("The number of gene sets found in gs.data is 0 or 1, no need to proceed")
        return(invisible(1))
    }
    
    if (is.numeric(limit)) {
        if (length(limit) == 1) 
            limit = c(-abs(limit), abs(limit))
        gs.data[gs.data > limit[2]] = limit[2]
        gs.data[gs.data < limit[1]] = limit[1]
    }
    
    nc = round(max(abs(range(gs.data))) * 100) * 2
    if (is.null(cols)) {
        cols = greenred(nc)
        cols = cols[max(1, round(sum(range(gs.data)) * 100)):min(nc, 
            nc + round(sum(range(gs.data)) * 100))]
    }
    
    heatmap2(gs.data, Colv = F, Rowv = F, dendrogram = "none", 
        col = cols, scale = "none", symkey = FALSE, density.info = "none", 
        trace = "none", margins = margins, keysize = 1, ...)
    
    return(invisible(1))
} 
