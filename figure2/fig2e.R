# Fig. 2e: WGD-centred networks of the retained alterations, pan-cancer and per cancer type. Boxes = retained copy-number
# alterations of one class and direction, black dots = mutated genes; green = co-occurrent with WGD, purple = mutually
# exclusive with WGD (CO right of the hub, ME left); * = retained with FDR <= 0.10 in one cohort and P < 0.05 in the other
# (tier2 in the table; tier1 = FDR <= 0.10 in both cohorts).
# The numbers under each title are WGD+ / WGD− tumours in MSK-IMPACT | TCGA. Layout in millimetres (coord_equal).
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
DAT <- "figure2/data"; OUT <- "results/figure2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
CO <- "#3A9E6F"; ME <- "#6A4C9C"; INK <- "#1A1A1A"
PT <- 5.2; LH <- PT * 1.12 / 72 * 25.4; CW <- PT * 0.50 / 72 * 25.4      # font pt -> line height / mean char width (mm)
PAD <- 1.0; GAP <- 0.9; HUB <- 3.2; DOT <- 0.9                           # box padding, slot spacing, hub-to-box gap, dot offset (mm)
hdr <- c(broad_gain = "Arm gain", broad_loss = "Arm loss", broad_loh = "Arm LOH", focal_amp = "Focal amp",
         focal_del = "Focal del", focal_loh = "Focal LOH")
sp <- read.delim(file.path(DAT, "fig2_scope_pairs.tsv"), check.names = FALSE)   # the retained associations, one row per alteration x scope
disp <- function(d) { x <- ifelse(is.na(d$display_MSK) | d$display_MSK == "", d$display_TCGA, d$display_MSK)   # focal LOH: cytoband (gene)
  ifelse(d$class == "focal_loh" & !is.na(x) & x != "", x, d$unit) }
sp$unit <- disp(sp)
r <- sp; r$dir <- r$direction_MSK
r$unit <- ifelse(r$tier == "tier2", paste0(r$unit, "*"), r$unit)
cnt <- read.delim(file.path(DAT, "fig2e_scope_counts.tsv"))                 # n(WGD+) / n(WGD-) per cohort and scope
ttl <- function(sc) { m <- cnt[cnt$scope == sc & cnt$cohort == "MSK-IMPACT", ]; t2 <- cnt[cnt$scope == sc & cnt$cohort == "TCGA", ]
  sprintf("%s\n%s / %s  |  %s / %s", ifelse(sc == "pancancer", "Pan-cancer", sc),
          format(m$wgd_pos, big.mark = ","), format(m$n - m$wgd_pos, big.mark = ","),
          format(t2$wgd_pos, big.mark = ","), format(t2$n - t2$wgd_pos, big.mark = ",")) }
band_group <- function(u) {                                              # "1p36 (CDC42)" "1p36 (ERRFI1)*" -> "1p36 (CDC42, ERRFI1*)"
  m <- regmatches(u, regexec("^(\\S+) \\((.+)\\)(\\*?)$", u)); band <- sapply(m, `[`, 2); gene <- paste0(sapply(m, `[`, 3), sapply(m, `[`, 4))
  ok <- !is.na(band); out <- u[!ok]
  if (any(ok)) { sp <- split(gene[ok], band[ok]); out <- c(out, sprintf("%s (%s)", names(sp), sapply(sp, function(x) paste(sort(x), collapse = ", ")))) }
  arm_sort(out) }
arm_sort <- function(u) {                                                 # 1p 1q 2p ... 22q, then gene symbols alphabetically
  ch <- suppressWarnings(as.integer(sub("^(\\d+).*$", "\\1", u))); u[order(is.na(ch), ch, u)] }
scopes <- c("pancancer", sort(setdiff(unique(cnt$scope), "pancancer"))); stopifnot(length(scopes) == 13)   # the layout below is for 12 cancer types

nodes_of <- function(sc, W) {                                            # W = panel width (mm); returns nodes with mm coords
  x <- r[r$scope == sc, ]; wrap <- floor((W / 2 - HUB - 2 * PAD - 0.5) / CW)
  cna <- x[x$class != "snv", ]
  b <- if (nrow(cna)) do.call(rbind, lapply(split(cna, list(cna$class, cna$dir), drop = TRUE), function(g) {
    if (g$class[1] == "focal_loh") {                                     # one cytoband per line: "1p36 (CDC42, ERRFI1, MAD2L2)"
      items <- unlist(lapply(band_group(g$unit), function(s) strwrap(s, wrap, exdent = 2)))
      txt <- c(paste0(hdr[g$class[1]], ":"), items)
    } else txt <- strwrap(paste0(hdr[g$class[1]], ": ", paste(arm_sort(g$unit), collapse = " ")), wrap)
    data.frame(type = "box", dir = g$dir[1], n = nrow(g), cls = g$class[1], label = paste(txt, collapse = "\n"),
               lines = length(txt), wchar = max(nchar(txt))) })) else NULL
  s <- x[x$class == "snv", ]
  g <- if (nrow(s)) data.frame(type = "gene", dir = s$dir, n = 1, cls = "snv", label = s$unit, lines = 1, wchar = nchar(s$unit)) else NULL
  n <- rbind(b, g); if (is.null(n) || !nrow(n)) return(NULL)
  n$h <- ifelse(n$type == "box", n$lines * LH + 2 * PAD + GAP, LH + 0.3)
  out <- lapply(c("CO", "ME"), function(dd) {
    m <- n[n$dir == dd, ]; if (!nrow(m)) return(NULL)
    m <- m[order(m$type, -m$n), ]                                        # boxes (largest set first), then genes
    mid <- cumsum(m$h) - m$h / 2; m$y <- -(mid - sum(m$h) / 2); m$x <- ifelse(dd == "CO", HUB, -HUB); m
  })
  n <- do.call(rbind, out); n$scope <- sc; n
}
draw <- function(sc, W, Hh, title = ttl(sc)) {                                  # Hh = half height (mm)
  n <- nodes_of(sc, W)
  if (is.null(n)) {                                                        # no retained alteration: title, n and the bare hub
    return(ggplot() + annotate("point", 0, 0, size = 4.2, shape = 21, fill = INK, colour = "white", stroke = 0.5) +
             annotate("text", 0, 0, label = "WGD", colour = "white", family = FONT, size = 4.6 / .pt, fontface = "bold") +
             annotate("text", 0, -4, label = "no replicated WGD-associated event", colour = "grey40", family = FONT, size = PT / .pt) +
             coord_equal(xlim = c(-W / 2, W / 2), ylim = c(-Hh, Hh), expand = FALSE, clip = "off") + labs(title = title) +
             theme_void(base_family = FONT) +
             theme(plot.title = element_text(size = 7.0, hjust = 0.5, colour = "grey10", lineheight = 0.95, margin = margin(b = 1.5)),
                   plot.margin = margin(1, 2, 1, 2)))
  }
  bx <- n[n$type == "box", ]; gn <- n[n$type == "gene", ]
  ggplot() +
    geom_segment(data = n, aes(x = 0, y = 0, xend = x, yend = y, colour = dir), linewidth = 0.45) +
    geom_label(data = bx, aes(x, y, label = label, colour = dir, hjust = ifelse(x > 0, 0, 1)), fill = "white", family = FONT,
               size = PT / .pt, linewidth = 0.3, label.padding = unit(PAD, "mm"), label.r = unit(0.5, "mm"), lineheight = 0.95,
               show.legend = FALSE) +
    geom_point(data = gn, aes(x, y), size = 1.0, colour = INK) +
    geom_text(data = gn, aes(x + sign(x) * DOT, y, label = label, hjust = ifelse(x > 0, 0, 1)), colour = INK, family = FONT,
              size = PT / .pt, fontface = "italic") +
    annotate("point", 0, 0, size = 4.2, shape = 21, fill = INK, colour = "white", stroke = 0.5) +
    annotate("text", 0, 0, label = "WGD", colour = "white", family = FONT, size = 4.6 / .pt, fontface = "bold") +
    scale_colour_manual(values = c(CO = CO, ME = ME), guide = "none") +
    coord_equal(xlim = c(-W / 2, W / 2), ylim = c(-Hh, Hh), expand = FALSE, clip = "off") +
    labs(title = title) + theme_void(base_family = FONT) +
    theme(plot.title = element_text(size = 7.0, hjust = 0.5, colour = "grey10", lineheight = 0.95, margin = margin(b = 1.5)),
          plot.margin = margin(1, 2, 1, 2))
}
half_h <- function(sc, W) { n <- nodes_of(sc, W); if (is.null(n)) 5 else max(abs(n$y) + n$h / 2) + 1 }
W_T <- 44; W_P <- 58; TT <- 11                                              # panel widths, title allowance (mm)
hh <- sapply(scopes[-1], half_h, W = W_T); hh_p <- half_h("pancancer", W_P)
small <- names(sort(hh))[1:3]; big <- names(sort(hh, decreasing = TRUE)); big <- big[!big %in% small]
rows <- split(big, rep(1:3, each = 3))                                       # right block: 3 rows x 3, tallest row first
row_h <- sapply(rows, function(r) max(hh[r]))
right <- wrap_plots(lapply(seq_along(rows), function(i) wrap_plots(lapply(rows[[i]], function(s) draw(s, W_T, row_h[i])), ncol = 3)),
                    ncol = 1, heights = 2 * row_h + TT)
left_h <- c(hh_p, hh[small])
left <- wrap_plots(c(list(draw("pancancer", W_P, hh_p)), lapply(small, function(s) draw(s, W_P, hh[s]))),
                   ncol = 1, heights = 2 * left_h + TT)
H_FIG <- max(sum(2 * left_h + TT), sum(2 * row_h + TT)) + 4
p <- (left | right) + plot_layout(widths = c(W_P, 3 * W_T))
ggsave(file.path(OUT, "fig2e.png"), p, width = W_P + 3 * W_T + 4, height = H_FIG, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig2e.pdf"), p, width = W_P + 3 * W_T + 4, height = H_FIG, units = "mm", device = cairo_pdf, bg = "white")
