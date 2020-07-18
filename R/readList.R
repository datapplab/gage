readList <- function(file) {
    f <- readLines(file)
    lst = sapply(f, function(x) unlist(strsplit(x, "\t", fixed = TRUE)))
    names(lst) = sapply(lst, function(x) x[1])
    lst = lapply(lst, function(x) x[-(1:2)])
    return(lst)
}
 
