
pick_maxlpd <- function(set){set <- set[!sapply(set, is.null)]; set[[which.max(sapply(set, function(x) tail(x$stat,1)[1]))]]}
get_imisiter <- function(x) nrow(x$stat)
get_logmargpost <- function(x) tail(x$stat, 1)[1]
get_nunique <- function(x) tail(x$stat, 1)[2]
