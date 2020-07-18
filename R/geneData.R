geneData <- function(genes, exprs, ref = NULL, samp = NULL, 
    outname = "array", txt = TRUE, heatmap = FALSE, scatterplot = FALSE, 
    samp.mean = FALSE, pdf.size = c(7, 7), cols = NULL, scale = "row", 
    limit = NULL, label.groups = TRUE, ...) {
    if (!is.null(ref) & is.null(samp)) 
        samp = (1:ncol(exprs))[-ref]
    if (is.null(ref) & !is.null(samp)) 
        ref = (1:ncol(exprs))[-samp]
    if (!is.null(ref) & !is.null(samp)) 
        icol = c(ref, samp)
    else icol = 1:ncol(exprs)
    
    genes = cbind(genes)
    sel = rownames(exprs) %in% genes[, 1]
    if (sum(sel) < 2) {
        print("The number of genes found in exprs is 0 or 1, no need to proceed")
        return(invisible(1))
    }
    gData = cbind(exprs[sel, ]) #exprs[sel, ]
    if (ncol(genes) > 1) 
        rownames(gData) = genes[match(rownames(gData), genes[, 
            1]), 2]
    if (txt) {
        filename = paste(outname, ".geneData.txt", sep = "")
        cat("Gene\t", paste(colnames(gData[, icol]), collapse = "\t"), 
            "\n", file = filename, sep = "")
        write.table(gData[, icol], file = filename, sep = "\t", 
            col.names = F, append = T)
    }
    
    if (heatmap & length(icol)>1 & !is.null(ref) & !is.null(samp)) {
        if (scale == "row") {
            gData.h = rownorm(gData[, icol])
        }
        else if (scale == "column") {
            gData.h = rownorm(t(gData[, icol]))
        }
        else gData.h = gData[, icol]
        
        if (is.numeric(limit)) {
            gData.h[gData.h > abs(limit)] = abs(limit)
            gData.h[gData.h < -abs(limit)] = -abs(limit)
        }
        
        nc = round(max(abs(range(gData.h))) * 100) * 2
        if (is.null(cols)) {
            cols = greenred(nc)
            cols = cols[max(1, round(sum(range(gData.h)) * 100)):min(nc, 
                nc + round(sum(range(gData.h)) * 100))]
        }
        
        pdf(paste(outname, ".geneData.heatmap.pdf", sep = ""), 
            width = pdf.size[1], height = pdf.size[2])
        if (label.groups & !is.null(ref) & !is.null(samp)) {
            col.cols = colorpanel(2, low = "black", high = "yellow")
            col.cols = rep(col.cols, c(length(ref), length(samp)))
            heatmap2(gData.h, col = cols, scale = "none", symkey = FALSE, 
                density.info = "none", trace = "none", ColSideColors = col.cols, 
                keysize = 1, ...)
        }
        else {
            heatmap2(gData.h, col = cols, scale = "none", symkey = FALSE, 
                density.info = "none", trace = "none", keysize = 1, ...)
        }
        dev.off()
    }
    
    if (scatterplot & !is.null(ref) & !is.null(samp)) {
        pdf(paste(outname, ".geneData.pdf", sep = ""), width = pdf.size[1], 
            height = pdf.size[2])
        sc1 = 1.5
        op = par(lwd = 2)
        if (samp.mean) {
            x = apply(gData[, ref], 1, mean, na.rm = T)
            y = apply(gData[, samp], 1, mean, na.rm = T)
            xlim = ylim = range(x, y)
            plot(x, y, type = "p", pch = 19, xlab = "Control Mean", 
                ylab = "Experiment Mean", xlim = xlim, ylim = ylim, 
                cex = sc1, cex.axis = sc1, cex.lab = sc1)
            abline(0, 1)
        }
        else {
            if (length(ref) > 1 & length(samp) > 1) {
                xlim = ylim = range(gData[, ref[1:2]], gData[, 
                  samp[1:2]])
            }
            else xlim = ylim = range(gData[, ref[1]], gData[, 
                samp[1]])
            plot(gData[, ref[1]], gData[, samp[1]], type = "n", 
                pch = 19, col = "gray", xlab = "Control", ylab = "Experiment", 
                xlim = xlim, ylim = ylim, cex = sc1, cex.axis = sc1, 
                cex.lab = sc1)
            abline(0, 1)
            points(gData[, ref[1]], gData[, samp[1]], pch = 19, 
                cex = sc1)
            if (length(ref) > 1 & length(samp) > 1) {
                points(gData[, ref[2]], gData[, samp[2]], pch = 24, 
                  col = "red", cex = sc1)
                legend("topleft", c("Sample 1", "Sample 2"), pch = c(19, 
                  24), col = c("black", "red"), bty = "n", pt.cex = sc1)
            }
        }
        par(op)
        dev.off()
    }
    
    return(invisible(1))
}
 
