gageComp <- function(sampnames, dataname, gsname = c("kegg.gs", 
    "go.gs"), use.cols = c("stat.mean", "q.val"), q.cutoff = 0.1, 
    do.plot = TRUE) {
    
    try(load(paste(dataname, ".gage.RData", sep = "")), silent = T)
    
    use.cols = paste("[,c('", use.cols[1], "','", use.cols[2], 
        "')]", sep = "")
    for (gs in gsname) {
        
        outname = paste(dataname, gs, sep = ".")
        dc = parse(text = paste("deComp(", paste(sampnames, paste(".", 
            gs, ".2d.p", sep = ""), "$greater", use.cols, sep = "", 
            collapse = ", "), ", sampnames=sampnames, outname=outname, q.cutoff=q.cutoff", 
            ")", sep = ""))
        assign(paste(gs, ".comp", sep = ""), eval(dc))
        
        outname = paste(dataname, gs, "up", sep = ".")
        dc = parse(text = paste("deComp(", paste(sampnames, paste(".", 
            gs, ".p", sep = ""), "$greater", use.cols, sep = "", 
            collapse = ", "), ", sampnames=sampnames, outname=outname, q.cutoff=q.cutoff", 
            ")", sep = ""))
        assign(paste(gs, ".up.comp", sep = ""), eval(dc))
        
        outname = paste(dataname, gs, "dn", sep = ".")
        dc = parse(text = paste("deComp(", paste(sampnames, paste(".", 
            gs, ".p", sep = ""), "$less", use.cols, sep = "", 
            collapse = ", "), ", sampnames=sampnames, outname=outname, q.cutoff=q.cutoff", 
            ")", sep = ""))
        assign(paste(gs, ".dn.comp", sep = ""), eval(dc))
        
    }
    
    if (do.plot & length(sampnames) < 4) {
        pdf(paste(dataname, "gage.comp.pdf", sep = "."))
        for (gs in gsname) {
            vennDiagram2(eval(as.name(paste(gs, ".comp", sep = ""))), 
                include = "both")
            vennDiagram2(eval(as.name(paste(gs, ".up.comp", sep = ""))), 
                eval(as.name(paste(gs, ".dn.comp", sep = ""))))
        }
        dev.off()
    }
    else if (length(sampnames) > 3) 
        warning("No venn diagram generated for comparison with >3 parties")
    
    return(invisible(1))
}
 
