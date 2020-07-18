vennCounts <- function(x, include = "both") {
    x <- as.matrix(x)
    include <- match.arg(include, c("both", "up", "down"))
    x <- sign(switch(include, both = abs(x), up = x > 0, down = x < 
        0))
    nprobes <- nrow(x)
    ncontrasts <- ncol(x)
    names <- colnames(x)
    if (is.null(names)) 
        names <- paste("Group", 1:ncontrasts)
    noutcomes <- 2^ncontrasts
    outcomes <- matrix(0, noutcomes, ncontrasts)
    colnames(outcomes) <- names
    for (j in 1:ncontrasts) outcomes[, j] <- rep(0:1, times = 2^(j - 
        1), each = 2^(ncontrasts - j))
    xlist <- list()
    for (i in 1:ncontrasts) xlist[[i]] <- factor(x[, ncontrasts - 
        i + 1], levels = c(0, 1))
    counts <- as.vector(table(xlist))
    structure(cbind(outcomes, Counts = counts), class = "VennCounts")
}
 
