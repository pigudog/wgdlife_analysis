# Fig. 5l: abnormal promoter hypermethylation of interferon-signalling, antigen-presentation (MHC-I / APM, MHC-II) and
# innate-sensing genes in WGD+ tumours: per gene, the lineage-stratified Cochran-Mantel-Haenszel odds ratio (95% CI) of a
# promoter methylation above the 95th percentile of the same lineage's WGD− tumours, WGD+ versus WGD−, in two columns
# (interferon pathway left, antigen presentation right); genes with fewer than two promoter probes are not drawn.
# Filled = BH FDR < 0.05 across the genes shown. Log-scale axis with the dashed null and a right-hand OR (95% CI) column.
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(patchwork) })
DAT <- "figure5/data"; OUT <- "results/figure5"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
BASE <- 7; TICK <- 6; LW <- 0.45
RED  <- "#C0504D"; GREY <- "grey55"

R <- fread(file.path(DAT, "fig5l_promoter_hypermethylation.tsv"))
MODS <- c("IFN signalling", "MHC-I / APM", "Innate sensing", "MHC-II")
R <- R[module %in% MODS]
R[, module := factor(module, MODS)]
R <- R[measurable == TRUE]
R[, sig := FDR < 0.05]
R[, txt := sprintf("%.2f (%.2f–%.2f)", OR, lo, hi)]
# two columns: left = IFN signalling over innate sensing, right = MHC-I / APM over MHC-II; within a module by OR
col_rows <- function(mods) {
  y <- 0; out <- list(); hd <- list()
  for (m in rev(mods)) {
    g <- R[module == m][order(OR)]
    g[, y := y + seq_len(.N)]; y <- max(g$y) + 1
    out[[m]] <- g; hd[[m]] <- data.table(module = m, y = y); y <- y + 1.3
  }
  list(G = rbindlist(out), H = rbindlist(hd))
}
XL <- c(0.62, 3.3); XT <- 3.75; XR <- 7.4
forest <- function(cr, ytop, lab_x = TRUE) {
  G <- cr$G; H <- cr$H
  ggplot(G, aes(y = y)) +
    geom_vline(xintercept = 1, linetype = "22", linewidth = LW, colour = "grey50") +
    annotate("segment", x = XL[1], xend = XL[2], y = 0.3, yend = 0.3, linewidth = LW, colour = "black") +
    geom_errorbar(aes(xmin = lo, xmax = hi, colour = sig), width = 0, linewidth = 0.5, orientation = "y") +
    geom_point(aes(x = OR, colour = sig, fill = sig), shape = 21, size = 1.3, stroke = 0.45) +
    geom_text(aes(x = 0.58, label = gene), hjust = 1, family = FONT, size = 2.15, colour = "black") +
    geom_text(aes(x = XT, label = txt), hjust = 0, family = FONT, size = 1.95, colour = "grey15") +
    geom_text(data = H, aes(x = 0.58, y = y, label = module), hjust = 1, family = FONT, size = 2.15, fontface = "bold") +
    annotate("text", x = XT, y = ytop - 0.6, label = "OR (95% CI)", hjust = 0, family = FONT, size = 1.95, colour = "grey40") +
    annotate("text", x = 1.03, y = ytop - 0.6, label = "more in WGD+ →", hjust = 0, family = FONT, size = 1.95, colour = "grey40") +
    scale_colour_manual(values = c(`TRUE` = RED, `FALSE` = GREY), guide = "none", aesthetics = "colour",
                        breaks = c(TRUE, FALSE)) +
    scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = "white"), guide = "none") +
    scale_x_continuous(transform = "log", breaks = c(0.8, 1, 1.5, 2, 3), labels = c("0.8", "1", "1.5", "2", "3"),
                       name = if (lab_x) "OR, abnormal promoter hypermethylation" else NULL) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(XL[1], XR), ylim = c(0.3, ytop), clip = "off") +
    labs(y = NULL) + theme_classic(base_family = FONT, base_size = BASE) +
    theme(axis.line = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.x = element_line(linewidth = LW), axis.text.x = element_text(size = TICK, colour = "black"),
          axis.title.x = element_text(size = TICK, margin = margin(t = 1), hjust = 0.25), plot.margin = margin(2, 2, 1, 40))
}
L <- col_rows(c("IFN signalling", "Innate sensing")); Rr <- col_rows(c("MHC-I / APM", "MHC-II"))
# top-align the two columns: shift the shorter column up so both first headers sit on the same line
sh <- max(Rr$H$y) - max(L$H$y); if (sh > 0) { L$G[, y := y + sh]; L$H[, y := y + sh] } else { Rr$G[, y := y - sh]; Rr$H[, y := y - sh] }
ytop <- max(c(L$H$y, Rr$H$y)) + 1.6
pl <- forest(L, ytop)
pr <- forest(Rr, ytop)
p <- (pl | pr) + plot_annotation(caption = "filled: BH FDR < 0.05",
                                theme = theme(plot.caption = element_text(family = FONT, size = TICK, colour = "grey40", hjust = 0)))
W <- 120; Hmm <- 66
ggsave(file.path(OUT, "fig5l.png"), p, width = W, height = Hmm, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig5l.pdf"), p, width = W, height = Hmm, units = "mm", device = cairo_pdf, bg = "white")
