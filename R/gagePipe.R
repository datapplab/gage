gagePipe <- function(arraydata, dataname = "arraydata", 
    trim.at = TRUE, sampnames, gsdata = NULL, gsname = c("kegg.gs", 
        "go.gs"), ref.list, samp.list, weight.list = NULL, comp.list = "paired", 
    q.cutoff = 0.1, heatmap = TRUE, pdf.size = c(7, 7), p.limit = c(0.5, 
        5.5), stat.limit = 5, ...) {
    
    if (length(arraydata) == 1 & is.character(arraydata)) {
        dataname = gsub("[.]R[dD]ata", "", basename(arraydata))
        load(arraydata)
        arraydata = eval(as.name(dataname))
    }
    if (trim.at & any(grep("_at", rownames(arraydata)[1:10]))) 
        rownames(arraydata) = gsub("_at", "", rownames(arraydata))
    
    if (is.null(dataname)) 
        stop("dataname needs to be specified")
    if (!is.null(gsdata)) 
        load(gsdata)
    for (gs in gsname) {
        if (!exists(gs)) 
            stop(paste(gs, "does not exsit"))
    }
    if (!is.list(samp.list)) 
        stop("samp.list needs to be a list")
    if (length(sampnames) != length(samp.list)) 
        stop("samp.list needs to match sampnames in length")
    if (is.list(ref.list)) {
        if (length(sampnames) != length(ref.list)) {
            if (length(ref.list) != 1) {
                stop("ref.list needs to match sampnames in length")
            }
            else ref.list = ref.list[[1]]
        }
    }
    
    #weight may be NULL, hence can not be unlisted
    if (is.list(weight.list)) {
        if (length(sampnames) != length(weight.list)) {
            if (length(weight.list) != 1) {
                stop("weight.list needs to match sampnames in length")
            }
            else weight.list = weight.list[[1]]
        }
    }
    
    #comp.list is different from weight.list
    if (is.list(comp.list)) 
        comp.list = unlist(comp.list)
    if (length(comp.list) == 1) 
        comp.list = rep(comp.list, length(sampnames))
    if (length(comp.list) != length(sampnames)) 
        stop("comp.list needs to match sampnames in length")
    
    for (i in 1:length(sampnames)) {
        
        ref = if (is.list(ref.list)) 
            ref.list[[i]]
        else ref.list
        weights = if (is.list(weight.list)) 
            weight.list[[i]]
        else weight.list
        comp = comp.list[i]
        for (gs in gsname) {
            assign(paste(sampnames[i], paste(gs, ".p", sep = ""), 
                sep = "."), gage(arraydata, gsets = eval(as.name(gs)), 
                ref = ref, samp = samp.list[[i]], weights = weights, 
                compare = comp, same.dir = T, ...))
            gage.p = eval(as.name(paste(sampnames[i], paste(gs, 
                ".p", sep = ""), sep = ".")))
            rn1 = rownames(gage.p$greater)[gage.p$greater[, "q.val"] < 
                q.cutoff & !is.na(gage.p$greater[, "q.val"])]
            if (length(rn1) > 0) {
                gage.p.sel1 = rbind(gage.p$greater[rn1, ])
                rownames(gage.p.sel1) = rn1
            }
            else gage.p.sel1 = NULL
            rn2 = rownames(gage.p$less)[gage.p$less[, "q.val"] < 
                q.cutoff & !is.na(gage.p$less[, "q.val"])]
            if (length(rn2) > 0) {
                gage.p.sel2 = rbind(gage.p$less[rn2, ])
                rownames(gage.p.sel2) = rn2
            }
            else gage.p.sel2 = NULL
            outdata = rbind(gage.p.sel1, gage.p.sel2)
            if (!is.null(outdata)) {
                filename = paste(dataname, paste(sampnames[i], 
                  paste(gs, ".p", sep = ""), sep = "."), "txt", 
                  sep = ".")
                cat("Gene Set\t", paste(colnames(eval(as.name(paste(sampnames[i], 
                  paste(gs, ".p", sep = ""), sep = ".")))$greater), 
                  collapse = "\t"), "\n", file = filename, sep = "")
                write.table(outdata, file = filename, sep = "\t", 
                  col.names = F, append = T)
                if (heatmap & nrow(outdata) > 1 & ncol(outdata) > 6) {
                  pdfname = paste(dataname, paste(sampnames[i], 
                    gs, sep = "."), "heatmap.pdf", sep = ".")
                  pdf(pdfname, width = pdf.size[1], height = pdf.size[2])
                  if (length(rn1) > 1) {
                    gs.heatmap(-log10(gage.p.sel1[, -c(1:5)]), 
                      limit = p.limit, main = "UP test: -log10(p-value)", 
                      ...)
                  }
                  else print(paste("No heatmap produced for up-regulated", 
                    gs, "gene sets, only 1 or none signficant."))
                  if (length(rn2) > 1) {
                    gs.heatmap(-log10(gage.p.sel2[, -c(1:5)]), 
                      limit = p.limit, main = "Down test: -log10(p-value)", 
                      ...)
                  }
                  else print(paste("No heatmap produced for down-regulated", 
                    gs, "gene sets, only 1 or none signficant."))
                  gs.heatmap(gage.p$stats[c(rn1, rn2), -1], limit = stat.limit, 
                    main = "Test statistics", ...)
                  dev.off()
                }
                else if (heatmap & ncol(outdata) > 6) 
                  print(paste("No heatmap produced for up- or down-regulated", 
                    gs, "gene sets, only 1 or none signficant."))
                else if (heatmap)
                  print(paste("No heatmap produced for up- or down-regulated", 
                    gs, "gene sets, due to single-column data"))
                
            }
            else print(paste("No", gs, "gene sets are signficant in one-direction!"))
            
            assign(paste(sampnames[i], paste(gs, ".2d.p", sep = ""), 
                sep = "."), gage(arraydata, gsets = eval(as.name(gs)), 
                ref = ref, samp = samp.list[[i]], weights = weights, 
                compare = comp, same.dir = F, ...))
            gage.p = eval(as.name(paste(sampnames[i], paste(gs, 
                ".2d.p", sep = ""), sep = ".")))
            rn1 = rownames(gage.p$greater)[gage.p$greater[, "q.val"] < 
                q.cutoff & !is.na(gage.p$greater[, "q.val"])]
            if (length(rn1) > 0) {
                gage.p.sel1 = rbind(gage.p$greater[rn1, ])
                rownames(gage.p.sel1) = rn1
                filename = paste(dataname, paste(sampnames[i], 
                  paste(gs, ".2d.p", sep = ""), sep = "."), "txt", 
                  sep = ".")
                cat("Gene Set\t", paste(colnames(eval(as.name(paste(sampnames[i], 
                  paste(gs, ".2d.p", sep = ""), sep = ".")))$greater), 
                  collapse = "\t"), "\n", file = filename, sep = "")
                write.table(gage.p.sel1, file = filename, sep = "\t", 
                  col.names = F, append = T)
                
                if (heatmap & length(rn1) > 1 & ncol(gage.p.sel1) > 6) {
                  pdfname = paste(dataname, paste(sampnames[i], 
                    gs, sep = "."), "2d.heatmap.pdf", sep = ".")
                  pdf(pdfname, width = pdf.size[1], height = pdf.size[2])
                  gs.heatmap(-log10(gage.p.sel1[, -c(1:5)]), 
                    limit = p.limit, main = "Two-way test: -log10(p-value)", 
                    ...)
                  gs.heatmap(gage.p$stats[rn1, -1], limit = stat.limit, 
                    main = "Test statistics", ...)
                  dev.off()
                }
                else if (heatmap & ncol(gage.p.sel1) > 6) 
                  print(paste("No heatmap produced for two-way perturbed", 
                    gs, "gene sets, only 1 signficant."))
                else if (heatmap)
                  print(paste("No heatmap produced for two-way perturbed", 
                    gs, "gene sets, due to single-column data"))
                
            }
            else print(paste("No", gs, "gene sets are signficant in two-direction!"))
            
        }
        
    }
    rm(gage.p)
    
    save(list = objects(pattern = "[.]p$"), file = paste(dataname, 
        ".gage.RData", sep = ""))
    return(invisible(1))
}
 
