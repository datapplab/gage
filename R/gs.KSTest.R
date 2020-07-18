gs.KSTest <- function(exprs, gsets, set.size = c(10, 
    500), same.dir = TRUE, ...) {
    if (class(gsets)[1] != "list") 
        stop("gsets need to be a list")
    exprs = apply(exprs, 2, rank)
    
    results <- matrix(NA, length(gsets), ncol(exprs))
    setsize <- rep(NA, length(gsets))
    names(setsize) <- rownames(results) <- names(gsets)
    colnames(results) <- colnames(exprs)
    mode(results) <- "numeric"
    nc <- ncol(exprs)
    nr <- nrow(exprs)
    p.results = ps.results = results
    mstat <- setsize
    for (i in 1:length(gsets)) {
        if (class(gsets[[i]])[1] == "smc") {
            clids <- gsets[[i]]@ids
        }
        else {
            clids <- gsets[[i]]
        }
        if (options()$verbose) 
            cat("Testing region ", i, "\n")
        ix <- match(clids, rownames(exprs))
        ix <- ix[!is.na(ix)]
        setsize[i] <- length(ix)
        present <- sum(!is.na(ix))
        if (present < set.size[1]) {
            if (options()$verbose) 
                cat("Skipping region ", i, " because too small-", 
                  present, ",\n")
            next
        }
        if (present > set.size[2]) {
            if (options()$verbose) 
                cat("Skipping region ", i, " because too large-", 
                  present, "\n")
            next
        }
        
        
        texprs <- matrix(exprs[ix, ], ncol = ncol(exprs))
        ks.res <- apply(texprs, 2, function(x) {
            res <- ks.test(x, seq_len(nr)[-x], alternative = "less")
            return(c(res$statistic, res$p.value))
        })
        p.results[i, ] <- ks.res[2, ]

        ks.res2 <- apply(texprs, 2, function(x) {
            res <- ks.test(x, seq_len(nr)[-x], alternative = "greater")
            return(c(res$statistic, res$p.value))
        })
        ps.results[i, ] <- ks.res2[2, ]

	results[i, ] <- mapply(max, ks.res[1, ], ks.res2[1, ])
        mstat[i] = mean(results[i, ])

    }
    if (!same.dir) 
        ps.results = NULL
    rawRes = list(results = results, p.results = p.results, ps.results = ps.results, 
        mstat = mstat, setsizes = setsize)
    return(rawRes)
}
 
