# Plotting bits
# -----------------------------------------------------------------------------
showPch <- function() {
    dev.new()
    plot(0:25, pch=0:25, col=1:25)
    text(0:25, labels=0:25, pos=3, xpd=T)
    text(0:25, labels=0:25, pos=3, xpd=T)
}

solarized <- c("#002b36", "#dc322f", "#b58900", "#268bd2", "#859900", "#6c71c4", "#d33682", "#2aa198")
harvard <- c(cod.gray="#0b0b09", vivid.burgundy="#961b36", medium.champagne="#f6eeab", tulip.tree="#eba938", monte.carlo="#7dc4ba", red.damask="#d96043", wattle="#d4d849", wheat="#f7deb2", light.blue="#b2dbde", fern.frond="#64821c", swans.down="#d8e9dc", chocolate="#d46619", allports="#206b87", red.damask="#d76340", boston.blue="#4683a8", aluminum="#989897", silver="#c2c2c1")
gruvstd <- c("#282828", "#CD3B27", "#98971B", "#D7992A", "#458588", "#B16286", "#689D6A", "#A89984")
gruvlt <- c("#928374", "#ED4631", "#B8BB26", "#F8BD32", "#83A598", "#D3869B", "#8EC07C", "#EBDBB2")
gruvGrey <- c("#1D2021", "#282828", "#3C3836", "#504945", "#665C54", "#7C6F64", "#928374", "#D65D2A")
gruvGrey2 <- c("#FFFFFF", "#32302F", "#A89984", "#BDAE93", "#D5C4A1", "#EBDBB2", "#FBF1C7", "#F17E30")

Kgrid <- function(bg = "white", cols = "gray93" ) {
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg,
         border = NA)
    xaxp <- par("xaxp"); yaxp <- par("yaxp")
    Vvec <- seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3])
    Hvec <- seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3])
    vvec <- seq(xaxp[1] + diff(Vvec)[1]/2, xaxp[2], by = abs(diff(Vvec)[1]))
    hvec <- seq(yaxp[1] + diff(Hvec)[1]/2, yaxp[2], by = abs(diff(Hvec)[1]))
    abline(v=Vvec, h = Hvec, lty=1, col = cols, lwd = 1)
    abline(v=vvec, h = hvec, lty=1, col = cols, lwd = .5)
}
vline <- function(v = 0, ...) abline(v = v, lty = 3, ...)
hline <- function(h = 0, ...) abline(h = h, lty = 3, ...)
put <- function(n.row, n.col, mar = NULL,...) {
  if (is.null(mar)) par(mfrow = c(n.row, n.col), mar = c(5,4,2,1)+0.1, ...)
  else par(mfrow = c(n.row, n.col), ...)
}
blankplot <- function(..., autoax = TRUE) {
  plot(..., type = "n", axes = FALSE)
  Kgrid()
  if (autoax) {
    Kaxis(1); Kaxis(2, las = 2) 
  }
}
lineplot <- function(..., autoax = TRUE) {
  plot(..., type = "n", axes = FALSE)
  Kgrid()
  if (autoax) {
    Kaxis(1); Kaxis(2, las = 2)
  }
  lines(...)
}
Kaxis <- function(side = 1, col='gray93', colticks='dimgray', ...) {
    axis(side, col=col, col.ticks=colticks, ...)
}
# Line plot of FOI, infections, death
plot.eppmix <- function(mod, which="FOI",...){
  palette(solarized)
  out <- attr(mod, which)
  if (length(dim(out)) !=3 ) stop("not implemented yet")
  put(dim(out)[2]/2, 2)
  mains <- c(sapply(c("Male", "Female"), paste, 1:(dim(out)[2]/2) ))
  invisible(sapply(1:dim(out)[2], function(y) {
    blankplot(1:52, numeric(52), xlab="Years", ylab=which, ylim=c(0, max(out)),
              main=mains[y])
    sapply(1:66, function(x) lines(out[x, y, ], col=x))
  }))
}
