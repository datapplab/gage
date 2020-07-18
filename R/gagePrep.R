gagePrep <- function(exprs, ref = NULL, samp = NULL, 
    same.dir = TRUE, compare = "paired", rank.test = FALSE, use.fold = TRUE, 
    weights = NULL, full.table = FALSE, ...) {
    if (class(exprs)[1] == "data.frame") 
        exprs <- as.matrix(exprs)
    if (is.numeric(exprs)) {
        exprs <- cbind(exprs)
        nc <- ncol(exprs)
        if(nc==1 & is.null(colnames(exprs))) colnames(exprs)='data'
        if(is.null(colnames(exprs))) stop("column names need to be specified")
      }
    else stop("exprs needs to be a numeric matrix or vector")
    if (!is.null(ref)) {
        if (!is.numeric(ref)) 
            stop("reference column index's required")
    }
    if (!is.null(ref)) {
        if (options()$verbose) 
            cat("Creating ratios...", "\n")
        if (!is.null(samp)) {
            if (!is.numeric(samp)) 
                stop("sample column index's required")
        }
        else samp = (1:nc)[-ref]
        if(length(ref)==0 | !all(ref %in% 1:nc)) stop("wrong reference column index")
        if(length(samp)==0 | !all(samp %in% 1:nc)) stop("wrong sample column index")
        
        if (use.fold) {
            if (compare == "as.group") {
                exprs = cbind(mean.fc=apply(cbind(exprs[, samp]), 1, 
                  mean) - apply(cbind(exprs[, ref]), 1, mean))
            }
            else if (compare == "paired") {
                if (!is.null(weights) & !all(c(length(ref), length(samp)) == 
                  length(weights))) 
                  stop("please make sure 'weights' comparable and of equal length to 'ref' and 'samp'")
                if (length(samp)%%length(ref) == 0 & length(intersect(samp, 
                  ref)) == 0) {
                  if (length(samp) > 1) 
                    exprs = exprs[, samp] - exprs[, ref]
                  else exprs = matrix(exprs[, samp] - exprs[, 
                    ref], ncol = length(samp), dimnames = list(rownames(exprs), 
                    colnames(exprs)[samp]))
                }
                else stop("please make sure 'ref' and 'samp' are comparable and of equal length or compare='unpaired'")
            }
            else if (compare == "unpaired" & (length(ref) * length(samp) > 
                1)) {
                exprs <- exprs[, rep(samp, each = length(ref))] - 
                  exprs[, rep(ref, length(samp))]
            }
            else {
                ref_mean <- (if (length(ref) > 1) 
                  apply(exprs[, ref], 1, mean, na.rm = TRUE)
                else exprs[, ref])
                cn <- colnames(exprs)[samp]
                exprs <- cbind(exprs[, samp] - ref_mean)
                colnames(exprs) <- cn
            }
        }
        else if (length(ref) > 1 & length(samp) > 1) {
            exprs = cbind(t.stat=apply(exprs, 1, function(x) t.test(x[samp], 
                x[ref], alternative = "two.sided", paired = (compare == 
                  "paired" & length(ref) == length(samp)))$statistic))
        }
        else stop("for t-test, please make sure 'ref' and 'samp' both have length >1")
    }
    else if (nc > 1) {
        if (use.fold) {
            if (compare == "as.group") {
                exprs = cbind(mean.fc=apply(exprs, 1, mean))
            }
            else if (compare == "paired") {
                if (!is.null(weights) & nc != length(weights)) 
                  stop("please make sure 'weights' comparable and equal to exprs column number")
            }
            else stop("improper 'compare' argument  value")
        }
        else {
            exprs = cbind(t.stat=apply(exprs, 1, function(x) t.test(x, 
                alternative = "two.sided", paired = F)$statistic))
            print("one sample t-test per gene")
        }
    }
    
#    exprs = cbind(exprs)
    if(ncol(exprs)==1) colnames(exprs)="exp1"
    if (!same.dir) 
        exprs = abs(exprs)
    if (rank.test) 
        exprs = apply(exprs, 2, rank)
    return(exprs)
}
 
