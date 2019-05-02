#' @importFrom grDevices adjustcolor dev.new palette
#' @importFrom graphics abline axis barplot boxplot layout legend lines 
#' @importFrom graphics mtext par plot points polygon rect segments text
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

#' @importFrom grDevices col2rgb rgb
AddAlpha <- function (plotclr, alpha = 0.5, verbose = 0) {
    tmp <- col2rgb(plotclr, alpha = alpha)
    tmp[4, ] = round(alpha * 255)
    for (i in 1:ncol(tmp)) {
        plotclr[i] = rgb(tmp[1, i], tmp[2, i], tmp[3, i], tmp[4, 
            i], maxColorValue = 255)
    }
    return(plotclr)
}

Kgrid <- function(bg = NA, cols = "gray93" ) {
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
  else par(mfrow = c(n.row, n.col), mar=mar, ...)
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
pointplot <- function(..., autoax = TRUE) {
  plot(..., type = "n", axes = FALSE)
  Kgrid()
  if (autoax) {
    Kaxis(1); Kaxis(2, las = 2)
  }
  points(...)
}
Kaxis <- function(side = 1, col='gray93', colticks='dimgray', ...) {
    axis(side, col=col, col.ticks=colticks, ...)
}

Kolygon <- function(x, y, ylow=NULL, col='gray70', alpha=0.4, border='white',...) {
  if (is.null(border)) border <- col
  if (missing(x)) x <- 1:length(y)
  if (missing(y)) {
    y <- x
    x <- seq(length(y))
  }
  xx <- c(x, rev(x))
  if (is.null(ylow)) yy <- c(y, rep(0, length(y)))
    else yy <- c(y, rev(ylow))
  polygon(xx, yy, col = AddAlpha(col, alpha), border = border, ...)
}

#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices  colorRampPalette
genCols <- function(n, pallete = "Set1") {
  getCols <- colorRampPalette(brewer.pal(8, pallete))
  return(getCols(n))
}

area_plot <- function(..., add=FALSE, autoax = TRUE) {
  if (!add) {
    plot(..., type = "n", axes = FALSE)
    Kgrid()
  }
  Kolygon(...)
  if (autoax) {
    Kaxis(1); Kaxis(2, las = 2) 
  }
}

# Line plot of FOI, infections, death
# -----------------------------------------------------------------------------
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

# plot demographic.
# -----------------------------------------------------------------------------
plot.dempp <- function(mod, start_year=1970, min_age=15, bin_year=5, bin_age=5, max_age=80, colors=solarized, ...){
  mod <- mod$data[,,1,]
  if (interactive()) dev.new(width=9, height=6)
  put(1, 2, c(5,4,2,2))
  palette(colors)
  mains  <- c("Male", "Female")
     yl  <- dim(mod)[3]
  y_name <- start_year:(start_year+yl-1)
  n_age  <- dim(mod)[1]
  cols   <- rep(1:floor(n_age/5), each=bin_age, length.out = n_age)
  at_age <- which(!duplicated(cols))
  at_yr  <- which(!duplicated(1:yl%/%bin_year))
  sapply(1:dim(mod)[2], function(y) {
    blankplot(1:yl, numeric(52), xlab="Years", ylab="Pop. size",
              ylim = c(0, max(mod)),
              main=mains[y], autoax=FALSE)
    Kaxis(1, at=at_yr, labels=y_name[at_yr])
    Kaxis(2)
    sapply(1:n_age, function(x) lines(mod[x, y, ], col=cols[x],...))
    text(yl+1, mod[at_age,y,yl], (min_age:max_age)[at_age], cex=.7)
  })
  invisible()
}

# Plot mod prevalence
#' @importFrom grDevices dev.new devAskNewPage
plot_prev_r6 <- function(mod, byAge=TRUE, byAgeGroup=FALSE, byYear=FALSE,
                      years=1970:2021, removeZero=TRUE, cols=c(4,5), stats=FALSE,
                      add=FALSE, stackbar = FALSE, colset = "Pastel1", 
                      separate=FALSE, debut=FALSE, ...) {
    palette(solarized)
    options(font.main=1)
    list2env(mod$ss, environment())
    sumByAGs <- function(x) x
    if (byAgeGroup)
      sumByAGs <- function(x) apply(x, 2, fastmatch::ctapply, ag.idx, sum)
    years <- as.character(years)
    inmod <- mod$data

    dimnames(inmod)[[4]] <- as.character(proj_start:(proj_start+PROJ_YEARS-1))
    ages <- AGE_START:(AGE_START+pAG-1)
    xlabs <- ages

    if (byAgeGroup) 
      xlabs <- c(ages[db_agr], paste(ages[agfirst.idx][-db_agr],
                                     ages[aglast.idx][-db_agr], sep='-'))

    if (!stackbar) {
      m.prev <- sumByAGs(inmod[,1,2,]) / 
                sumByAGs(rowSums(aperm(inmod[,1,,], c(1,3,2)),,2))
      f.prev <- sumByAGs(inmod[,2,2,]) / 
                sumByAGs(rowSums(aperm(inmod[,2,,], c(1,3,2)),,2))
    }

    if (stackbar) {
        m.prev <- sweep(inmod[,1,2,], 2, colSums(inmod[,1,2,]), '/')
        f.prev <- sweep(inmod[,2,2,], 2, colSums(inmod[,2,2,]), '/')
        m.prev[is.na(m.prev)] <- 0
        f.prev[is.na(f.prev)] <- 0
    }

    if (debut) {
        byAge <- FALSE
        indb <- mod$pop_db
        dimnames(indb)[[4]] <- as.character(proj_start:(proj_start+PROJ_YEARS-1))
        # prevalence among total population
        m.prev <- sumByAGs(inmod[db_aid,1,2,years]) / 
                sumByAGs(rowSums(aperm(inmod[db_aid,1,,years], c(1,3,2)),,2))
        f.prev <- sumByAGs(inmod[db_aid,2,2,years]) / 
                sumByAGs(rowSums(aperm(inmod[db_aid,2,,years], c(1,3,2)),,2))
        # prevalence among sexually population
        a_pop  <- inmod[db_aid,,,years]
        a_pop  <- a_pop - indb[,,,years]
        m.prev_a <- sumByAGs(a_pop[,1,2,]) / 
                    sumByAGs(rowSums(aperm(a_pop[,1,,], c(1,3,2)),,2))
        f.prev_a <- sumByAGs(a_pop[,2,2,]) / 
                    sumByAGs(rowSums(aperm(a_pop[,2,,], c(1,3,2)),,2))
        # incidence rate among sexually active population vs. total population
        # zero for inactive
        # % sexually active by age
        pc_active <- 1 - rowSums(aperm(indb[,,,years],c(1,2,4,3)),,3) / 
                         rowSums(aperm(inmod[db_aid,,,years],c(1,2,4,3)),,3)
    }

    # Beer's coffiecients

    if (removeZero) {
        start <- which(colSums(m.prev) > 0)[1]
        end   <- dim(m.prev)[2]
        identical(start, which(colSums(f.prev) > 0)[1])
        m.prev <- m.prev[, start:end]
        f.prev <- f.prev[, start:end]
        if (debut) {
          m.prev_a <- m.prev_a[, start:end]
          f.prev_a <- f.prev_a[, start:end]
        }
    }

    if (debut) {
      ylim <- c(0, max(m.prev, f.prev, m.prev_a, f.prev_a))
      def.par <- par(no.readonly = TRUE)
      layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), heights=c(.6,.4))
      sexdebutplot(m.prev, m.prev_a, ages, years, ylim, db_aid, main="Male") 
      sexdebutplot(f.prev, f.prev_a, ages, years, ylim, db_aid, main="Female") 
      lineplot(pc_active[,1,1], col=years[1], lwd=2, ylab="% active", xlab="", cex=.7, autoax=F, ylim=c(0,1))
      Kaxis(1, at=db_aid, labels=ages[db_aid])
      Kaxis(2)
      lines(pc_active[,1,2], col=years[2], lty=2, lwd=2)
      lines(pc_active[,1,3], col=years[3], lty=3, lwd=2)
      lines(pc_active[,1,4], col=years[4], lty=4, lwd=2)
      lineplot(pc_active[,2,1], col=years[1], lwd=2, ylab="% active", xlab="", cex=.7, autoax=F, ylim=c(0,1))
      Kaxis(1, at=db_aid, labels=ages[db_aid])
      Kaxis(2)
      lines(pc_active[,2,2], col=years[2], lty=2, lwd=2)
      lines(pc_active[,2,3], col=years[3], lty=3, lwd=2)
      lines(pc_active[,2,4], col=years[4], lty=4, lwd=2)
      par(def.par)
    }

    if (stackbar) {
        byAge <- FALSE
        xlabs <- proj_start:(proj_start+PROJ_YEARS-1)
        xlabs <- xlabs[start:end]
        cls <- genCols(66, colset)
        if (separate) {
          put(1,1,c(4,2,1,2))
          stackbarplot(m.prev, xlabs, ages, cls, main="Male")
          if (interactive()) devAskNewPage(ask = 1)
          stackbarplot(f.prev, xlabs, ages, cls, main="Female")
        } else {
          put(1,2,c(4,2,1,2))
          stackbarplot(m.prev, xlabs, ages, cls, main="Male")
          stackbarplot(f.prev, xlabs, ages, cls, main="Female")
        }
        if (stats)
          return(list(f.prev, m.prev))
    }

    if (byAge) {
        ncol   <- dim(t(m.prev))[2]
        m.prev <- boxplot(t(m.prev), plot=FALSE)
        f.prev <- boxplot(t(f.prev), plot=FALSE)
        ylim <- c(0, max(m.prev$stats, f.prev$stats))
        if (!add) {
            blankplot(m.prev$stats[3, ], ylim=ylim, autoax=F, xlab="Age", ylab="Prevalence (IQR)", ...)
            Kaxis(1, at=1:ncol, labels=xlabs); Kaxis(2)
        }
        Kolygon(y=m.prev$conf[2,], ylow=m.prev$conf[1,], col=cols[1])
        Kolygon(y=f.prev$conf[2,], ylow=f.prev$conf[1,], col=cols[2])
        lines(m.prev$stats[3,], col=cols[1], lwd=2)
        lines(f.prev$stats[3,], col=cols[2], lwd=2)
        legend("topright", lwd=2, col=cols, legend= c("Male", "Female"), bty='n')
    }
    if (stats) return(list(m.prev, f.prev))
    on.exit(devAskNewPage(ask = 0))
}

stackbarplot <- function(x, xlabs, ages, cls=cls, space=0, axes=F, xaxs='i',
                         ats=4, border=AddAlpha('gray95'), ...) {
  tm <- barplot(x, col=cls, space=space, border=border, axes=axes, xaxs=xaxs,
                font.main=1, ...)
  ncol <- length(tm)
  atYear <- seq(1, ncol, ats)
  atAge  <- seq(1, length(ages), ats)
  text(tm[atYear], par("yaxp")[1]-.05, labels=xlabs[atYear], srt=90,
       xpd=T, cex=.7)
  text(tm[1]-2, cumsum(x[,1])[atAge], labels=ages[atAge],
       cex=.7, xpd=T, col='dimgray')
  text(tm[ncol]+2, cumsum(x[,ncol])[atAge], labels=ages[atAge],
       cex=.7, xpd=T, col='dimgray')
}

sexdebutplot <- function(m.prev, m.prev_a, ages, years, ylim, db_aid,...) {
    blankplot(m.prev[, 1], ylim=ylim, autoax=F, xlab="Age", ylab="Prevalence",
              font.main=1, ...)
    Kaxis(1, at=db_aid, labels=ages[db_aid])
    Kaxis(2)
    for (i in years) {
        lines(m.prev[, i], col = i, lty=1, lwd=2)
        lines(m.prev_a[, i], col = i, lty=3, lwd=2)
    }
    legend("topleft", lwd=2, col=years, cex=.7, legend=years, bty='n')
    legend("bottomright", lwd=2, lty=c(1,3), cex=.7,
           legend=c("Total", "Sexual active"), bty='n')
}