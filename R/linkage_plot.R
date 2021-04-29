linkage_plot <- function (mapthis, outfile, mapthese = NULL, at.axis = NULL, 
          autoconnadj = TRUE, cex.axis = par("cex.axis"), cex.lgtitle = par("cex.main"), 
          cex.main = par("cex.main"), col.axis = par("col.axis"), 
          col.lgtitle = par("col.main"), col.main = par("col.main"), 
          conndf = NULL, denmap = FALSE, dupnbr = FALSE, font.axis = par("font.axis"), 
          font.lgtitle = par("font.main"), font.main = par("font.main"), 
          header = TRUE, labdist = 0.3, labels.axis = TRUE, lcex = par("cex"), 
          lcol = par("col"), lfont = par("font"), lgperrow = NULL, 
          lgtitles = NULL, lgw = 0.25, lg.col = NULL, lg.lwd = par("lwd"), 
          lty.axis = "solid", lwd.axis = 1, lwd.ticks.axis = lwd.axis, 
          main = NULL, markerformatlist = NULL, maxnbrcolsfordups = 3, 
          pdf.bg = "transparent", pdf.family = "Helvetica", pdf.fg = "black", 
          pdf.width = NULL, pdf.height = NULL, pdf.pointsize = 12, 
          pdf.title = "LinkageMapView R output", posonleft = NULL, 
          prtlgtitles = TRUE, qtldf = NULL, revthese = NULL, rcex = par("cex"), 
          rcol = par("col"), rfont = par("font"), roundpos = 1, rsegcol = TRUE, 
          ruler = FALSE, sectcoldf = NULL, segcol = NULL, qtlscanone = NULL, 
          showonly = NULL, units = "cM", ylab = units) 
{
  pgx <- 0.5
  pgy <- 0.5
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (!is.null(roundpos)) {
    if (!is.numeric(roundpos)) {
      stop("roundpos - number of places after decimal point - must be numeric")
    }
  }
  font.axis <- convertfont("font.axis", font.axis)
  font.lgtitle <- convertfont("font.lgtitle", font.lgtitle)
  font.main <- convertfont("font.main", font.main)
  lfont <- convertfont("lfont", lfont)
  rfont <- convertfont("rfont", rfont)
  if (!is.null(markerformatlist$font)) {
    markerformatlist$font <- convertfont("markerformatlist$font", 
                                         markerformatlist$font)
  }
  if ("cross" %in% class(mapthis)) {
    if (is.null(mapthese)) {
      mapthese <- qtl::chrnames(mapthis)
    }
    lgin <- readlgcross(mapthis, mapthese)
  }
  else if ("character" %in% class(mapthis)) {
    if (is.null(mapthese)) {
      lgin <- readlgtext(mapthis, header = header)
      mapthese <- unique(lgin$group)
    }
    else {
      lgin <- readlgtext(mapthis, mapthese, header = header)
    }
  }
  else if ("data.frame" %in% class(mapthis)) {
    if (is.null(mapthese)) {
      lgin <- readlgdf(as.data.frame(mapthis))
      mapthese <- unique(lgin$group)
    }
    else {
      lgin <- readlgdf(mapthis, mapthese)
    }
  }
  else {
    stop("first parameter, mapthis, must be a data frame, a filename, or an r/qtl cross object")
  }
  if (!is.null(revthese)) {
    if (!all(revthese %in% mapthese)) {
      stop("Requested chrnames to reverse not in those requested to plot")
    }
  }
  if (!is.null(lgperrow)) {
    if (!(is.numeric(lgperrow) && floor(lgperrow) == lgperrow)) {
      stop("lgperrow must be NULL or an integer")
    }
  }
  else {
    lgperrow <- length(mapthese)
  }
  nbrrows <- ceiling(length(mapthese)/lgperrow)
  if (denmap) {
    autoconnadj <- FALSE
    rsegcol <- FALSE
    ruler <- TRUE
    showonly <- NULL
    markerformatlist <- NULL
    conndf <- NULL
  }
  if (!is.null(segcol)) {
    if (!(segcol %in% colnames(lgin))) {
      stop(c("segcol column ", segcol, " not found in mapthis"))
    }
    else {
      rsegcol <- FALSE
    }
  }
  if (!is.null(showonly)) {
    if (!all(showonly %in% lgin$locus)) {
      stop(c("Requested showonly locus not found:", showonly[!(showonly %in% 
                                                                 lgin$locus)]))
    }
    else {
      showonly <- unique(showonly)
    }
  }
  else {
    if (any(lgin$locus == "") || any(is.na(lgin$locus))) {
      notnull <- lgin$locus[which(lgin$locus != "")]
      showonly <- notnull[which(!is.na(notnull))]
    }
  }
  if (!is.null(qtldf)) {
    fas <- sapply(qtldf, is.factor)
    qtldf[fas] <- lapply(qtldf[fas], as.character)
  }
  if (!is.null(qtlscanone)) {
    if (!("data.frame" %in% class(qtlscanone) & "scanone" %in% 
          class(qtlscanone))) {
      stop(c("qtlscanone should be a data.frame for r/qtl"))
    }
    else {
      qtldf <- usescanone(qtlscanone, qtldf, mapthese, 
                          pdf.fg, maxdec = roundpos)
    }
  }
  if (is.null(sectcoldf) && denmap == TRUE) {
    sectcoldf <- lmvdencolor(lgin)
  }
  if (!is.null(sectcoldf)) {
    fas <- sapply(sectcoldf, is.factor)
    sectcoldf[fas] <- lapply(sectcoldf[fas], as.character)
  }
  if (!is.null(posonleft)) {
    if (!(all(posonleft %in% c(T, F)) && length(posonleft) == 
          length(mapthese))) {
      stop("Position on vector must be same length as # of linkage groups to map and only TRUE or FALSE")
    }
  }
  else {
    posonleft <- rep(TRUE, length.out = length(mapthese))
  }
  if (!(is.null(ruler)) && !(length(ruler) == 1)) {
    stop("ruler must be a single TRUE or FALSE")
  }
  if (!(ruler == T) && !(ruler == F)) {
    stop("ruler may only be TRUE or FALSE")
  }
  if (!(is.null(dupnbr)) && !(length(dupnbr) == 1)) {
    stop("dupnbr must be a single TRUE or FALSE")
  }
  if (!(dupnbr == T) && !(dupnbr == F)) {
    stop("dupnbr may only be TRUE or FALSE")
  }
  if (dupnbr == T) {
    maxnbrcolsfordups <- 2
  }
  if (!(prtlgtitles == T) && !(prtlgtitles == F)) {
    stop("prtlgtitles may only be TRUE or FALSE")
  }
  if (ruler) {
    leftmar <- 1 + ceiling(cex.axis)
    if (!is.null(units)) {
      leftmar <- leftmar + 1
    }
  }
  else {
    leftmar <- 1
  }
  if (prtlgtitles) {
    topmar <- ceiling(cex.lgtitle) + 1
  }
  else {
    topmar <- 1
  }
  if (!is.null(main)) {
    topoma <- 1 + ceiling(cex.main)
  }
  else {
    topoma <- 1
  }
  pdf.options(bg = pdf.bg, title = pdf.title, family = pdf.family, 
              pointsize = pdf.pointsize, fg = pdf.fg)
  on.exit(pdf.options(reset = TRUE), add = TRUE)
  pdf(outfile, width = 30, height = 30)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, leftmar, topmar, 0), oma = c(1, 0, topoma, 
                                              1), new = FALSE)
  lg <- list()
  dim <- list()
  lrcol <- list()
  llcol <- list()
  lrfont <- list()
  llfont <- list()
  lrcex <- list()
  llcex <- list()
  qtl <- list()
  solist <- list()
  sectcol <- list()
  for (i in 1:length(mapthese)) {
    lg[[i]] <- getlg(lgin, mapthese[i], dupnbr, roundpos)
    editlgdf(lg[[i]])
    if (lg[[i]]$group[1] %in% revthese) {
      lg[[i]]$locus <- rev(lg[[i]]$locus)
      lg[[i]]$position <- revpos(lg[[i]]$position, roundpos)
      if (!is.null(segcol)) {
        lg[[i]][[eval(segcol)]] <- rev(lg[[i]][[eval(segcol)]])
      }
      if (!is.null(sectcoldf)) {
        sectcol[[i]] <- rev(subset(sectcoldf, sectcoldf$chr == 
                                     lg[[i]][1, 1]))
      }
    }
    else {
      lg[[i]]$position <- round(lg[[i]]$position, roundpos)
      if (!is.null(sectcoldf)) {
        sectcol[[i]] <- subset(sectcoldf, sectcoldf$chr == 
                                 lg[[i]][1, 1])
      }
    }
    if (!is.null(qtldf)) {
      qtl[[i]] <- subset(qtldf, qtldf$chr == lg[[i]][1, 
                                                     1])
      for (qtlix in 1:nrow(qtl[[i]])) {
        if (nrow(subset(qtldf, qtldf$chr == lg[[i]][1, 
                                                    1])) != 0) {
          if (qtl[[i]]$chr[qtlix] %in% revthese) {
            so <- qtl[[i]]$so[qtlix]
            si <- qtl[[i]]$si[qtlix]
            ei <- qtl[[i]]$ei[qtlix]
            eo <- qtl[[i]]$eo[qtlix]
            qtl[[i]]$eo[qtlix] <- max(lg[[i]]$position) - 
              so - min(lg[[i]]$position)
            qtl[[i]]$ei[qtlix] <- max(lg[[i]]$position) - 
              si - min(lg[[i]]$position)
            qtl[[i]]$si[qtlix] <- max(lg[[i]]$position) - 
              ei - min(lg[[i]]$position)
            qtl[[i]]$so[qtlix] <- max(lg[[i]]$position) - 
              eo - min(lg[[i]]$position)
          }
        }
      }
    }
    if (!is.null(sectcoldf)) {
      sectcol[[i]] <- subset(sectcoldf, sectcoldf$chr == 
                               lg[[i]][1, 1])
    }
  }
  miny <- vector(length = nbrrows)
  maxy <- vector(length = nbrrows)
  for (nr in 1:nbrrows) {
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese), nr * lgperrow))
    miny[nr] <- min(lg[[fromlg]]$position)
    maxy[nr] <- max(lg[[fromlg]]$position)
    for (i in fromlg:tolg) {
      if (max(lg[[i]]$position) > maxy[nr]) {
        maxy[nr] <- max(lg[[i]]$position)
      }
      if (min(lg[[i]]$position) < miny[nr]) {
        miny[nr] <- min(lg[[i]]$position)
      }
    }
  }
  for (i in 1:length(mapthese)) {
    lrcol[[i]] <- rcol
    lrcol[[i]] = rep(lrcol[[i]], length.out = length(lg[[i]]$locus))
    lrfont[[i]] <- rfont
    lrfont[[i]] = rep(lrfont[[i]], length.out = length(lg[[i]]$locus))
    lrcex[[i]] <- rcex
    lrcex[[i]] <- rep(lrcex[[i]], length.out = length(lg[[i]]$locus))
    llcol[[i]] <- lcol
    llcol[[i]] = rep(llcol[[i]], length.out = length(lg[[i]]$position))
    llfont[[i]] <- lfont
    llfont[[i]] = rep(llfont[[i]], length.out = length(lg[[i]]$position))
    llcex[[i]] <- lcex
    llcex[[i]] <- rep(llcex[[i]], length.out = length(lg[[i]]$position))
    if (!is.null(markerformatlist)) {
      for (ml in 1:length(markerformatlist)) {
        if (!is.null(markerformatlist[[ml]]$col) && 
            !is.na(markerformatlist[[ml]]$col)) {
          lrcol[[i]][which(lg[[i]]$locus %in% markerformatlist[[ml]]$locus)] <- markerformatlist[[ml]]$col
        }
        if (!is.null(markerformatlist[[ml]]$cex) && 
            !is.na(markerformatlist[[ml]]$cex)) {
          lrcex[[i]][which(lg[[i]]$locus %in% markerformatlist[[ml]]$locus)] <- markerformatlist[[ml]]$cex
        }
        if (!is.null(markerformatlist[[ml]]$font) && 
            !is.na(markerformatlist[[ml]]$font)) {
          lrfont[[i]][which(lg[[i]]$locus %in% markerformatlist[[ml]]$locus)] <- markerformatlist[[ml]]$font
        }
      }
    }
    if (!is.null(qtldf)) {
      qtldfone <- qtl[[i]]
      if (nrow(qtl[[i]]) == 0) {
        qtldfone <- NULL
      }
    }
    else {
      qtldfone <- NULL
    }
    dim[[i]] <- reqdim(lg[[i]], c(miny[nr], maxy[nr]), denmap = denmap, 
                       maxnbrcolsfordups = maxnbrcolsfordups, pdf.width = 30, 
                       labdist = labdist, lcol = llcol[[i]], lfont = llfont[[i]], 
                       lcex = llcex[[i]], rcol = lrcol[[i]], rfont = lrfont[[i]], 
                       rcex = lrcex[[i]], cex.lgtitle = cex.lgtitle, qtldf = qtldfone, 
                       ruler = ruler, prtlgtitles = prtlgtitles, showonly = showonly)
  }
  pgxlg <- vector(length = length(mapthese))
  width <- vector(length = length(mapthese))
  totwidth <- rep(0, nbrrows)
  totheight <- rep(0, nbrrows)
  for (nr in 1:nbrrows) {
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese), nr * lgperrow))
    dim[[fromlg]]$reqwidth <- dim[[fromlg]]$reqwidth + 0.3
    for (i in fromlg:tolg) {
      totwidth[nr] <- dim[[i]]$reqwidth + totwidth[nr]
      if (dim[[i]]$reqheight > totheight[nr]) {
        totheight[nr] <- dim[[i]]$reqheight
      }
    }
  }
  allrowwidth <- max(totwidth)
  if (denmap & !is.null(sectcoldf$dens)) {
    allrowheight <- sum(totheight) + 1.5
    relheight <- vector(length = (nbrrows + 1))
    relheight[(nbrrows + 1)] <- 1.5/allrowheight
  }
  else {
    allrowheight <- sum(totheight)
    relheight <- vector(length = nbrrows)
  }
  for (nr in 1:nbrrows) {
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese), nr * lgperrow))
    maxrowheight <- 0
    for (i in fromlg:tolg) {
      if (dim[[i]]$reqheight > maxrowheight) {
        maxrowheight <- dim[[i]]$reqheight
      }
    }
    relheight[nr] <- maxrowheight/allrowheight
  }
  allrowwidth <- allrowwidth + par("mai")[2] + par("mai")[4] + 
    par("omi")[2] + par("omi")[4]
  allrowheight <- allrowheight + par("omi")[2] + par("omi")[4]
  if (is.null(pdf.width)) {
    pdf.width <- ceiling(allrowwidth)
  }
  if (is.null(pdf.height)) {
    pdf.height <- ceiling(allrowheight)
  }
  message(c("Required pdf.width = ", allrowwidth))
  message(c("Required pdf.height = ", allrowheight))
  message(c("Using pdf.width = ", pdf.width))
  message(c("Using pdf.height = ", pdf.height))
  dev.off()
  pdf.options(bg = pdf.bg, title = pdf.title, family = pdf.family, 
              pointsize = pdf.pointsize, fg = pdf.fg)
  pdf(outfile, width = pdf.width, height = pdf.height)
  par(mar = c(0, leftmar, topmar, 0), oma = c(1, 0.2, topoma, 
                                              0.2), new = FALSE)
  if (denmap & !is.null(sectcoldf$dens)) {
    layout(c(seq(1, nbrrows + 1)), heights = relheight)
  }
  else {
    layout(c(seq(1, nbrrows)), heights = relheight)
  }
  for (nr in 1:nbrrows) {
    plot(0.5, 0.5, xlim = c(0, 1), ylim = c(maxy[nr], miny[nr]), 
         type = "n", cex = 1, xaxt = "n", yaxt = "n", xlab = "", 
         ylab = "", xaxs = "i", bty = "n")
    pin <- par("pin")[1]
    if (ruler) {
      axis(side = 2, col.axis = col.axis, cex.axis = cex.axis, 
           font.axis = font.axis, at = at.axis, labels = labels.axis, 
           lty = lty.axis, lwd = lwd.axis, lwd.ticks = lwd.ticks.axis)
      mtext(ylab, side = 2, line = leftmar - 1)
    }
    if (nr == nbrrows) {
      mtext("Rendered by LinkageMapView", side = 1, cex = 0.5, 
            outer = TRUE, col = pdf.fg, adj = 0)
    }
    yrlabwidth <- list()
    adjyr <- list()
    adjyl <- list()
    dups <- list()
    widthused <- 0
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese), nr * lgperrow))
    for (i in fromlg:tolg) {
      pctwidth <- dim[[i]]$reqwidth/totwidth[nr]
      width <- pctwidth * pin
      if (posonleft[i]) {
        if (!ruler) {
          llablen <- dim[[i]]$maxlenllab
          lspace <- strwidth("M", units = "inches") * 
            max(lcex)
          llabdist <- labdist
        }
        else {
          if (i > 1 && !posonleft[i - 1]) {
            llablen <- 0
            lspace <- max(strwidth(lg[[i]][1, 1]) * 
                            cex.lgtitle, strwidth(lg[[i - 1]][1, 1]) * 
                            cex.lgtitle)/2 + strwidth("M", units = "inches") * 
              cex.lgtitle
            llabdist <- 0
          }
          else {
            llablen <- 0
            lspace <- 0
            llabdist <- 0
          }
        }
      }
      else {
        llablen <- dim[[i]]$maxlenrlab
        lspace <- strwidth("M", units = "inches") * 
          max(rcex)
        llabdist <- labdist
      }
      pgxlg[i] <- (sum(llablen, lgw/2, llabdist, lspace))/pin + 
        widthused/pin
      if (i == (nr - 1) * lgperrow + 1) {
        pgxlg[i] <- pgxlg[i] + 0.3/pdf.width
      }
      widthused <- widthused + width
      if (!is.null(qtldf)) {
        qtldfone <- qtl[[i]]
        if (nrow(qtl[[i]]) == 0) {
          qtldfone <- NULL
        }
      }
      else {
        qtldfone <- NULL
      }
      if (!is.null(sectcoldf)) {
        sectcoldfone <- sectcol[[i]]
        if (nrow(sectcol[[i]]) == 0) {
          sectcoldfone <- NULL
        }
      }
      else {
        sectcoldfone <- NULL
      }
      if (!is.null(lgtitles)) {
        lgtitleone <- lgtitles[i]
      }
      else {
        lgtitleone <- NULL
      }
      if (i == 1) {
        main <- main
      }
      else {
        main = NULL
      }
      dolist <- drawone(lg[[i]], dim[[i]], totwidth[nr], 
                        c(miny[nr], maxy[nr]), denmap = denmap, maxnbrcolsfordups = maxnbrcolsfordups, 
                        pdf.width = pin, pdf.fg = pdf.fg, lgw = lgw, 
                        lg.col = lg.col, lg.lwd = lg.lwd, pgx = pgxlg[i], 
                        labdist = labdist, lcol = llcol[[i]], lfont = llfont[[i]], 
                        lcex = llcex[[i]], rcol = lrcol[[i]], rfont = lrfont[[i]], 
                        rcex = lrcex[[i]], rsegcol = rsegcol, main = main, 
                        cex.main = cex.main, font.main = font.main, 
                        col.main = col.main, cex.lgtitle = cex.lgtitle, 
                        font.lgtitle = font.lgtitle, col.lgtitle = col.lgtitle, 
                        qtldf = qtldfone, posonleft = posonleft[i], 
                        ruler = ruler, prtlgtitles = prtlgtitles, lgtitles = lgtitleone, 
                        segcol = segcol, showonly = showonly, sectcoldf = sectcoldfone)
      yrlabwidth[[i]] <- dolist$yrlabwidth
      adjyr[[i]] <- dolist$adjyr
      adjyl[[i]] <- dolist$adjyl
      dups[[i]] <- dolist$dups
      solist[[i]] <- dolist$solist
    }
    if (autoconnadj == TRUE) {
      autoconndf <- autoconn(lg[fromlg:tolg], pdf.fg, 
                             lgperrow)
      if (!is.null(conndf)) {
        fas <- sapply(conndf, is.factor)
        conndf[fas] <- lapply(conndf[fas], as.character)
        allconndf <- rbind(conndf, autoconndf)
        allconndf <- allconndf[!duplicated(allconndf[, 
                                                     1:4]), ]
      }
      else {
        allconndf <- autoconndf
      }
    }
    else {
      if (!is.null(conndf)) {
        fas <- sapply(conndf, is.factor)
        conndf[fas] <- lapply(conndf[fas], as.character)
        allconndf <- conndf
      }
      else {
        allconndf <- data.frame()
      }
    }
    if (nrow(allconndf) > 0) {
      for (i in 1:nrow(allconndf)) {
        fromi <- match(allconndf$fromchr[i], mapthese)
        if (is.na(fromi)) {
          stop(c("Connect marker from position not found for chr = ", 
                 allconndf$fromchr[i], " and locus = ", allconndf$fromlocus[i]))
        }
        toi <- match(allconndf$tochr[i], mapthese)
        if (is.na(toi)) {
          stop(c("Connect marker to position not found for chr = ", 
                 allconndf$tochr[i], " and locus = ", allconndf$tolocus[i]))
        }
        if (toi < fromi) {
          temp <- toi
          toi <- fromi
          fromi <- temp
        }
        if (ceiling(fromi/lgperrow) != ceiling(toi/lgperrow)) {
          stop(c("Connect marker from chr ", allconndf$fromchr[i], 
                 " must be on the same row as to chr ", allconndf$tochr[i], 
                 " and you have lgperrow = ", lgperrow))
        }
        fpos <- lg[[fromi]]$position[lg[[fromi]]$locus == 
                                       allconndf$fromlocus[i]]
        if (identical(fpos, numeric(0))) {
          stop(c("Connect marker from position not found for chr = ", 
                 allconndf$fromchr[i], " and locus = ", allconndf$fromlocus[i]))
        }
        tpos <- lg[[toi]]$position[lg[[toi]]$locus == 
                                     allconndf$tolocus[i]]
        if (identical(tpos, numeric(0))) {
          stop(c("Connect marker to position not found for chr = ", 
                 allconndf$tochr[i], " and locus = ", allconndf$tolocus[i]))
        }
        if (!is.null(showonly)) {
          fy <- match(fpos, solist[[fromi]]$newllab)
          ty <- match(tpos, solist[[toi]]$newllab)
        }
        else {
          fy <- match(fpos, lg[[fromi]]$position)
          ty <- match(tpos, lg[[toi]]$position)
        }
        if (posonleft[fromi]) {
          connfxpos <- ((yrlabwidth[[fromi]][fy] + lgw/2 + 
                           labdist/2)/pin) + pgxlg[fromi]
          if (!is.null(showonly)) {
            connfypos <- adjyr[[fromi]][match(fpos, 
                                              solist[[fromi]]$newllab[dups[[fromi]]$rkeep])]
          }
          else {
            connfypos <- adjyr[[fromi]][match(fpos, 
                                              lg[[fromi]]$position[dups[[fromi]]$rkeep])]
          }
        }
        else {
          if (!ruler) {
            connfxpos <- ((dim[[fromi]]$maxlenllab + 
                             lgw/2 + labdist + strwidth("M", units = "inches"))/pin) + 
              pgxlg[fromi]
            if (!is.null(showonly)) {
              connfypos <- adjyl[[fromi]][match(fpos, 
                                                solist[[fromi]]$newllab[dups[[fromi]]$lkeep])]
            }
            else {
              connfypos <- adjyl[[fromi]][match(fpos, 
                                                lg[[fromi]]$position[dups[[fromi]]$lkeep])]
            }
          }
          else {
            connfxpos <- ((lgw/2)/pin) + pgxlg[fromi]
            connfypos <- fpos
          }
        }
        if (!posonleft[toi]) {
          conntxpos <- -((yrlabwidth[[toi]][ty] + lgw/2 + 
                            labdist/2)/pin) + pgxlg[toi]
          if (!is.null(showonly)) {
            conntypos <- adjyr[[toi]][match(tpos, solist[[toi]]$newllab[dups[[toi]]$rkeep])]
          }
          else {
            conntypos <- adjyr[[toi]][match(tpos, lg[[toi]]$position[dups[[toi]]$rkeep])]
          }
        }
        else {
          if (!ruler) {
            conntxpos <- -((dim[[toi]]$maxlenllab + 
                              lgw/2 + labdist + strwidth("M", units = "inches"))/pin) + 
              pgxlg[toi]
            if (!is.null(showonly)) {
              conntypos <- adjyl[[toi]][match(tpos, 
                                              solist[[toi]]$newllab[dups[[toi]]$lkeep])]
            }
            else {
              conntypos <- adjyl[[toi]][match(tpos, 
                                              lg[[toi]]$position[dups[[toi]]$lkeep])]
            }
          }
          else {
            conntxpos <- -((lgw/2)/pin) + pgxlg[toi]
            conntypos <- tpos
          }
        }
        if (is.null(allconndf$col[i])) {
          allconndf$col[i] = pdf.fg
        }
        segments(connfxpos, connfypos, conntxpos, conntypos, 
                 col = allconndf$col[i])
      }
    }
  }
  if (denmap & !is.null(sectcoldf$dens)) {
    leg <- sectcoldf[order(sectcoldf$dens, na.last = NA), 
                     ]
    if (max(leg$dens) > 100) {
      leg$dens <- round(leg$dens, digits = 0)
    }
    else {
      leg$dens <- round(leg$dens, digits = 1)
    }
    uleg <- leg[!duplicated(leg$dens), ]
    if ((pdf.width/0.25) < nrow(uleg)) {
      nbrbuckets <- round(pdf.width/0.25, digits = 0)
      bplotdens <- uleg$dens[seq(1, length(uleg$dens), 
                                 length(uleg$dens)/nbrbuckets)]
      bplotcol <- uleg$col[seq(1, length(uleg$col), length(uleg$col)/nbrbuckets)]
    }
    else {
      bplotdens <- uleg$dens
      bplotcol <- uleg$col
    }
    if (!bplotdens[length(bplotdens)] == uleg$dens[length(uleg$dens)]) {
      bplotdens <- append(bplotdens, uleg$dens[length(uleg$dens)])
      bplotcol <- append(bplotcol, uleg$col[length(uleg$col)])
    }
    par(mar = c(5, leftmar, 1, 1))
    barplot(rep(1, length(bplotcol)), col = bplotcol, space = 0, 
            axes = F, xlab = sprintf("Density (Loci/%s)",units), 
            names = 1/bplotdens, cex.names = 0.75, 
            las = 2, cex.lab = 0.75)
  }
}
