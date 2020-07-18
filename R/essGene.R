essGene <- function(gs, exprs, ref = NULL, samp = NULL, 
    gsets = NULL, compare = "paired", use.fold = TRUE, rank.abs = FALSE, 
    use.chi = FALSE, chi.p = 0.05, ...) {
    if (is.null(gsets)) 
        a = gs
    else {
        if (class(gsets[[gs]])[1] == "smc") 
            a = gsets[[gs]]@ids
        else a = gsets[[gs]]
    }
    a = unique(a[a %in% rownames(exprs)])
    expdata = exprs
    
    exprs = gagePrep(exprs, ref = ref, samp = samp, compare = compare, 
        use.fold = use.fold, same.dir = T, ...)
    b = cbind(exprs[a, ])
    nc = ncol(b)
    ind = apply(b, 1, mean)
    if (rank.abs) 
        ind = abs(ind)
    b = cbind(b[order(ind, decreasing = T), ])
    means = apply(exprs, 2, mean, na.rm = T)
    sds = apply(exprs, 2, sd, na.rm = T)
    if (use.chi) {
        chis = apply(((t(b) - means)/sds)^2, 2, sum)
        sel = chis > qchisq(chi.p, nc, lower.tail = F)
    }
    else sel = apply(abs((t(b) - means)/sds), 2, mean) > 1
    if (!is.null(ref) & !is.null(samp)) {
        return(expdata[rownames(b)[sel], c(ref, samp)])
    }
    else return(expdata[rownames(b)[sel], ])
}
 
