sigGeneSet <- function(setp, cutoff = 0.1, dualSig = (0:2)[2], 
    qpval = c("q.val", "p.val")[1], heatmap = TRUE, outname = "array", 
    pdf.size = c(7, 7), p.limit = c(0.5, 5.5), stat.limit = 5, 
    ...) {
    if (is.list(setp) & length(setp) == 3) {
        gs.name = rownames(setp$greater)
        setp$less = setp$less[gs.name, ]
        greater = setp$greater[, "p.val"] < setp$less[, "p.val"]
        sel.greater.0 = !is.na(setp$greater[, qpval]) & setp$greater[, 
            qpval] < cutoff
        sel.less.0 = !is.na(setp$less[, qpval]) & setp$less[, 
            qpval] < cutoff
        if (dualSig == 0) {
            sel.greater = sel.greater.0 & !sel.less.0
            sel.less = sel.less.0 & !sel.greater.0
        }
        else if (dualSig == 1) {
            sel.greater = sel.greater.0 & greater
            sel.less = sel.less.0 & !greater
        }
        else if (dualSig == 2) {
            sel.greater = sel.greater.0
            sel.less = sel.less.0
        }
        else {
            print("incorrect value for dualSig arguement")
        }
        
        ord = order(setp$less[sel.less, "p.val"])
        setp.sig = setp
        setp.sig$greater = rbind(setp$greater[sel.greater, ])
        setp.sig$less = rbind(setp$less[gs.name[sel.less][ord], 
            ])
        sig.gs = c(gs.name[sel.greater], gs.name[sel.less][ord])
        setp.sig$stats = rbind(setp.sig$stats[sig.gs, ])
        rownames(setp.sig$greater) = gs.name[sel.greater]
        rownames(setp.sig$less) = gs.name[sel.less][ord]
        rownames(setp.sig$stats) = sig.gs
        
        if (heatmap & length(sig.gs) > 1) {
            pdf(paste(outname, ".gs.heatmap.pdf", sep = ""), 
                width = pdf.size[1], height = pdf.size[2])
            if (sum(sel.greater) > 1) {
                gs.heatmap(-log10(setp.sig$greater[, -c(1:5)]), 
                  limit = p.limit, main = "GAGE Up-test: -log10(p-value)", 
                  ...)
            }
            else print("No heatmap produced for up-regulated gene sets, only 1 or none signficant.")
            if (sum(sel.less) > 1) {
                gs.heatmap(-log10(setp.sig$less[, -c(1:5)]), 
                  limit = p.limit, main = "  GAGE Down-test: -log10(p-value)", 
                  ...)
            }
            else print("No heatmap produced for down-regulated gene sets, only 1 or none signficant.")
            gs.heatmap(setp.sig$stats[, -1], limit = stat.limit, 
                main = "GAGE test statistics", ...)
            dev.off()
        }
        else if (heatmap) 
            print("No heatmap produced for up- or down-regulated gene sets, only 1 or none signficant.")
        
        print(paste("there are", sum(sel.greater), "signficantly up-regulated gene sets"))
        print(paste("there are", sum(sel.less), "signficantly down-regulated gene sets"))
        return(setp.sig)
    }
    else if (is.list(setp) & length(setp) == 2) {
        
        gs.name = rownames(setp$greater)
        sel = !is.na(setp$greater[, qpval]) & setp$greater[, 
            qpval] < cutoff
        
        setp.sig = setp
        setp.sig$greater = rbind(setp$greater[sel, ])
        setp.sig$stats = rbind(setp.sig$stats[sel, ])
        rownames(setp.sig$greater) = gs.name[sel]
        rownames(setp.sig$stats) = gs.name[sel]
        
        if (heatmap & sum(sel) > 1) {
            pdf(paste(outname, ".gs.2d.heatmap.pdf", sep = ""), 
                width = pdf.size[1], height = pdf.size[2])
            gs.heatmap(-log10(setp.sig$greater[, -c(1:5)]), limit = p.limit, 
                main = "GAGE Two-way test: -log10(p-value)", ...)
            gs.heatmap(setp.sig$stats[, -1], limit = stat.limit, 
                main = "GAGE test statistics", ...)
            dev.off()
        }
        else if (heatmap) 
            print("No heatmap produced for two-way perturbed gene sets, only 1 or none signficant.")
        print(paste("there are", sum(sel), "signficantly two-direction perturbed gene sets"))
        return(setp.sig)
    }
    else print("setp is an invalid input!")
}
 
