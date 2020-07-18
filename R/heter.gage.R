heter.gage <- function(exprs, gsets, ref.list, samp.list, 
    comp.list = "paired", use.fold = TRUE, ...) {
    pData = pairData(exprs, ref.list = ref.list, samp.list = samp.list, 
        comp.list = comp.list, use.fold = use.fold, ...)
    hg.p = gage(pData$exprs, gsets, ref = NULL, samp = NULL, 
        weights = pData$weights, ...)
    return(hg.p)
}
 
