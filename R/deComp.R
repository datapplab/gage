deComp <- function(..., sampnames, outname = NULL, 
    org = NULL, q.cutoff = 0.1, common = FALSE, return.all = FALSE) {
    input = list(...)
    dims = lapply(input, dim)
    if (!all(sapply(dims, length) == 2)) 
        stop("all inputs need to be matrices")
    if (!all(sapply(dims, identical, dims[[1]]))) 
        stop("all input matrices need to be of the same size")
    if (length(sampnames) != length(input)) 
        stop("all input matrices need to be of the same size")
    
    rn = rownames(input[[1]])
    rn.eq = sapply(input, function(x) {
        return(all(rownames(x) %in% rn))
    })
    if (!all(rn.eq)) 
        stop("all input matrices need to have the same set of rownames")
    
    il = length(input)
    if (il < 2) 
        stop("need to have more than 1 input matrices")
    
    comb = input[[1]]
    for (i in 2:il) {
        comb = cbind(comb, input[[i]][rn, ])
    }
    
    colnames(comb) = paste(rep(sampnames, each = 2), c("stat.mean", 
        "q.value"), sep = "_")
    pp = apply(comb[, c(1:il) * 2], 1, prod)
    nsig = apply(comb[, c(1:il) * 2] > q.cutoff, 1, sum)
    comb = cbind(comb, hits = il - nsig)[order(nsig, pp), ]
    comb = comb[!is.na(comb[, "hits"]), ]
    sel = comb[, "hits"] > 0
    if (!is.null(outname) & sum(sel) > 0) {
        filename = paste(outname, ".comb.txt", sep = "")
        outdata = rbind(comb[sel, ])
        if (sum(sel) == 1) 
            rownames(outdata) = rownames(comb)[sel]
        if (!is.null(org)) {
            annot = NULL
            if (exists("annot.db")) {
                annot = annot.db(rownames(outdata), org = org)
            }
            else if (all(is.numeric(rownames(outdata)))) {
                annot = cbind(rownames(outdata), eg2sym(rownames(outdata)))
            }
            else annot = rownames(outdata)
            write.table(cbind(annot, data.frame(outdata)), file = filename, 
                row.names = F, sep = "\t")
        }
        else {
            cat("Gene Set\t", paste(colnames(comb), collapse = "\t"), 
                "\n", file = filename, sep = "")
            write.table(outdata, file = filename, sep = "\t", 
                col.names = F, append = T)
        }
    }
    else if (sum(sel) == 0) 
        print(paste("No genes or gene sets are signficant at q.value <", 
            q.cutoff))
    
    if (common) 
        return(rownames(comb)[comb[, "hits"] == il])
    
    qs = comb[, (1:il) * 2]
    qs = 1 - (qs > q.cutoff)
    vc = vennCounts(qs)
    colnames(vc) = gsub("_q.value", "", colnames(vc))
    if (!is.null(outname)) {
        return(vc)
    }
    else if (!return.all) {
        return(list(comb = comb[sel, ], vc = vc))
    }
    else {
        return(list(comb = comb, vc = vc))
    }
}
 
