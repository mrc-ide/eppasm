# Natural age to index
a2i <- function(x, min=15, max=80) which(min:max %in% x)

# Non number to value
na2num <- function(x, y) {x[is.na(x)] <- y; return(x)}