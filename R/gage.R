gage <- function(exprs, gsets, ref = NULL, samp = NULL, 
    set.size = c(10, 500), same.dir = TRUE, compare = "paired", 
    rank.test = FALSE, use.fold = TRUE, FDR.adj = TRUE, weights = NULL, 
    full.table = FALSE, saaPrep = gagePrep, saaTest = gs.tTest, 
    saaSum = gageSum, use.stouffer = TRUE, ...) {
    exprs = saaPrep(exprs, ref = ref, samp = samp, same.dir = same.dir, 
        compare = compare, rank.test = rank.test, use.fold = use.fold, 
        weights = weights)
    rawRes = saaTest(exprs, gsets = gsets, set.size = set.size, 
        same.dir = same.dir)
    greater.res = saaSum(rawRes, ref = ref, same.dir = same.dir, 
        compare = compare, use.fold = use.fold, weights = weights, 
        full.table = full.table, use.stouffer = use.stouffer)
    if (same.dir) {
        less.res = saaSum(rawRes, ref = ref, test4up = FALSE, 
            same.dir = same.dir, compare = compare, use.fold = use.fold, 
            weights = weights, full.table = full.table, use.stouffer = use.stouffer)
        return(list(greater = greater.res$p.glob, less = less.res$p.glob, 
            stats = greater.res$results))
    }
    else return(list(greater = greater.res$p.glob, stats = greater.res$results))
}
 
