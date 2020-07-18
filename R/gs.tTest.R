gs.tTest <- function(exprs, gsets, set.size = c(10, 
    500), same.dir = TRUE, ...) {
    if (class(gsets)[1] != "list") 
        stop("gsets need to be a list")
    
    s = apply(exprs, 2, sd, na.rm = TRUE)^2
    mu = apply(exprs, 2, mean, na.rm = TRUE)
    results <- matrix(NA, length(gsets), ncol(exprs))
    setsize <- rep(NA, length(gsets))
    names(setsize) <- rownames(results) <- names(gsets)
    colnames(results) <- colnames(exprs)
    mode(results) <- "numeric"
    nc <- ncol(exprs)
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
        
        a = apply(texprs, 2, sd, na.rm = TRUE)^2/length(ix)
        b = s/length(ix)
        df = (a + b)^2/(a^2/(length(ix) - 1) + b^2/(length(ix) - 
            1))
        mod <- (a + b)^(-0.5)
        m <- apply(texprs, 2, mean, na.rm = TRUE) - mu
        stat <- m * mod
        results[i, ] <- as.numeric(stat)
        mstat[i] = mean(as.numeric(stat))
        p.results[i, ] <- pt(as.numeric(stat), df, lower.tail = FALSE)
        ps.results[i, ] <- pt(as.numeric(stat), df, lower.tail = TRUE)
    }
    if (!same.dir) 
        ps.results = NULL
    rawRes = list(results = results, p.results = p.results, ps.results = ps.results, 
        mstat = mstat, setsizes = setsize)
    return(rawRes)
}
 
