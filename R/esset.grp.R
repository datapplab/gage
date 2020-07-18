esset.grp <- function(setp, exprs, gsets, ref = NULL, 
    samp = NULL, test4up = TRUE, same.dir = TRUE, compare = "paired", 
    use.fold = TRUE, cutoff = 0.01, use.q = FALSE, pc = 10^-10, 
    output = TRUE, outname = "esset.grp", make.plot = FALSE, 
    pdf.size = c(7, 7), core.counts = FALSE, get.essets = TRUE, 
    bins = 10, bsize = 1, cex = 0.5, layoutType = "circo", name.str = c(10, 
        100), ...) {
    
    exprs = gagePrep(exprs, ref = ref, samp = samp, same.dir = same.dir, 
        compare = compare, use.fold = use.fold)
    
    if (use.q) {
        ssp = setp[!is.na(setp[, "q.val"]) & setp[, "q.val"] < 
            cutoff, ]
    }
    else ssp = setp[!is.na(setp[, "p.val"]) & setp[, "p.val"] < 
        cutoff, ]
    if (length(ssp) < 2) 
        stop("There are less than 1 significant gene set, try to increase the cutoff (P-value)")
    if (is.null(nrow(ssp))) 
        stop("There are only 1 significant gene set, no redundant gene sets need to be combined")
    ssname = rownames(ssp)
    m = apply(exprs, 1, mean, na.rm = T)
    s = sd(m, na.rm = T)
    m = mean(m, na.rm = T)
    
    core = vector("list", nrow(ssp))
    names(core) = ssname
    sign = test4up
    
    for (i in 1:nrow(ssp)) {
        if (class(gsets[[ssname[i]]]) == "smc") 
            a = gsets[[ssname[i]]]@ids
        else a = gsets[[ssname[i]]]
        a = a[a %in% rownames(exprs)]
        b = apply(cbind(exprs[a, ]), 1, mean) - m
        b = b[order(b, decreasing = sign)]
        core[[i]] = names(b[if (sign | !same.dir) b > s else b < 
            (-s)])
    }
    len.core = sapply(core, length)
    b = apply(cbind(exprs), 1, mean) - m
    b = b[order(b, decreasing = sign)]
    ess.genes = names(b[if (sign | !same.dir) b > 2 * s else b < 
        (-2 * s)])
    
    pmat = matrix(NA, ncol = nrow(ssp), nrow = nrow(ssp))
    olmat = pmat
    ind = cbind(rep(1:nrow(ssp), each = nrow(ssp)), rep(1:nrow(ssp), 
        nrow(ssp)))
    pv = apply(ind, 1, function(x) {
        overlap = length(intersect(core[[x[1]]], core[[x[2]]]))
        if (overlap < 1) 
            return(c(overlap, 1))
        else {
            return(c(overlap, phyper(overlap, len.core[x[1]], 
                length(ess.genes) - len.core[x[1]], len.core[x[2]], 
                lower.tail = F, log.p = FALSE)))
        }
    })
    pmat[ind] = pv[2, ]
    olmat[ind] = pv[1, ]
    
    b <- pmat
    b[b > pc] <- 1
    b[b < pc] <- 0
    b = 1 - b
    diag(b) = 0
    rownames(b) = colnames(b) = 1:nrow(b)
    IDs = ssname
    names(IDs) = rownames(b)
    b[is.na(b)] <- 0
    
    sim.Graph <- as(b, "graphNEL")
    sim.Graph.cc <- connComp(sim.Graph)
    setGroup = lapply(sim.Graph.cc, function(x) ssname[as.numeric(x)])
    eset = sapply(setGroup, function(x) x[1])
    iset = sapply(setGroup, function(x) x[which.max(len.core[x])])
    if (output) {
        outdata = cbind(rbind(ssp[eset, ]), setGroup = sapply(setGroup, 
            function(x) paste(x, sep = "", collapse = "; ")))
        filename = paste(outname, ".esgp.txt", sep = "")
        cat("essentialSets\t", paste(colnames(outdata), collapse = "\t"), 
            "\n", file = filename, sep = "")
        write.table(outdata, file = filename, sep = "\t", col.names = F, 
            append = T)
    }
    
    if (make.plot & exists("make.graph")) {
        ##gograph generation
        colnames(b) = rownames(b) = ssname
        gg <- as(b, "graphNEL")
        goes = setp
        nn = length(nodes(gg))
        gn = substr(nodes(gg), name.str[1], name.str[2])
        
        
        if (length(dim(goes)) == 2) {
            ll = nodes(gg)
            sel = ll %in% rownames(goes)
            if (core.counts) 
                ccs = paste("; ", sapply(core[ll[sel]], length), 
                  sep = "")
            else ccs = NULL
            gn[sel] = paste(gn[sel], " (", round(-log10(goes[ll[sel], 
                "p.val"]), 1), "; ", goes[ll[sel], "set.size"], 
                ccs, ")", sep = "")
            ps = rep(NA, nn)
            names(ps) = nodes(gg)
            ps[sel] = -log10(goes[ll[sel], "p.val"])
            if (!is.null(bins)) 
                ps[ps > (bins * bsize)] = bins * bsize
        }
        else {
            ps = NULL
        }
        
        #edge (or overlap) p-values
        eps <- -log10(pmat[lower.tri(pmat)])
        eidx = eps > -log10(pc)
        eps <- eps[eidx]
        eps <- eps + log10(pc) + 1
        eps[eps > (bins * bsize)] = bins * bsize
        names(eps) <- edgeNames(gg)
        #edge weights (overlap counts)
        ew <- olmat[lower.tri(olmat)]
        ew <- ew[eidx]
        names(ew) <- edgeNames(gg)
        
        
        ##make.graph
        pdf(paste(outname, ".pdf", sep = ""), width = pdf.size[1], 
            height = pdf.size[2])
        grg.xy = make.graph(gg, gn, ps, eps, ew, bins = bins, 
            cex = cex, layoutType = layoutType, ...)
        
        ##mark essential nodes
        if (get.essets && length(dim(goes)) == 2) {
            maxi = ps[eset]
            
            mxy = lapply(grg.xy, function(x) x[nodes(gg) %in% 
                names(maxi)])
            points(mxy, pch = 3, cex = 2 * cex, col = "blue")
            ineset = iset[!iset %in% eset]
            if (length(ineset) > 0) {
                ixy = lapply(grg.xy, function(x) x[nodes(gg) %in% 
                  ineset])
                points(ixy, pch = 0, cex = 2 * cex, col = "blue")
            }
        }
        
        dev.off()
        
    }
    
    return(list(essentialSets = eset, setGroups = setGroup, allSets = ssname, 
        connectedComponent = sim.Graph.cc, overlapCounts = olmat, 
        overlapPvals = pmat, coreGeneSets = core))
}
 
