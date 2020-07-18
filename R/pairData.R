pairData <- function(exprs, ref.list, samp.list, comp.list = "paired", 
    use.fold = TRUE, ...) {
    if (!is.list(ref.list) | !is.list(samp.list)) 
        stop("ref.list and samp.list both need to be lists")
    refs = unlist(ref.list)
    samps = unlist(samp.list)
    if (any(duplicated(refs))) 
        stop("ref.list list contains overlapping vectors")
    if (any(duplicated(samps))) 
        stop("samp.list list contains overlapping vectors")
    if (length(intersect(refs, samps)) > 0) 
        stop("ref.list and samp.list lists overlap")
    if (any(sapply(ref.list, length) == 0) | length(refs) == 
        0) 
        stop("ref.list list contains empty vector(s)")
    if (any(sapply(samp.list, length) == 0) | length(samps) == 
        0) 
        stop("samp.list list contains empty vector(s)")
    lref = length(ref.list)
    lsamp = length(samp.list)
    if (lref != lsamp) 
        stop("ref.list and samp.list lists need to match in length")
    
    if (is.list(comp.list)) 
        comp.list = unlist(comp.list)
    if (!all(comp.list %in% c("paired", "unpaired"))) 
        stop("comp.list values need to be either 'paired' or 'unpaired'")
    if (length(comp.list) == 1) 
        comp.list = rep(comp.list, lref)
    if (length(comp.list) != lref) 
        stop("comp.list needs to match in length with ref.list and samp.list")
    
    expData = weights = NULL
    for (i in 1:lref) {
        exprs2 = gagePrep(exprs, ref = ref.list[[i]], samp = samp.list[[i]], 
            compare = comp.list[i], ...)
        expData = cbind(expData, exprs2)
        if (use.fold & ("unpaired" %in% comp.list)) {
            #if(use.fold){
            if (comp.list[i] == "unpaired") 
                dr = length(ref.list[[i]])
            if (comp.list[i] == "paired") 
                dr = 1
            weights = c(weights, rep(1/dr, length(samp.list[[i]]) * 
                dr))
        }
    }
    return(list(exprs = expData, weights = weights))
}
 
