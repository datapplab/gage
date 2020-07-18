rownorm <- function(x) {
    dat.z <- apply(x, 1, function(y) (y - mean(y, na.rm = T))/sd(y, 
        na.rm = T))
    dat.z <- t(dat.z)
    dat.z
}
 
