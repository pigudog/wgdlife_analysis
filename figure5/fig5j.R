# Fig. 5j: where WGD-associated hypermethylation lands. Per carcinoma lineage, the share of CpGs that lie on PRC2-target
# CGIs among all CpGs tested in that lineage and among the CpGs inside its WGD hyper-DMRs (the within-lineage DMR calls),
# one line per lineage, boxes summarising the 12 lineages; paired Wilcoxon signed-rank test across lineages.
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(data.table) })
DAT <- "figure5/data"; OUT <- "results/figure5"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
PAL <- c(NSCLC = "#4E79A7", Breast = "#F28E2B", Colorectal = "#59A14F", Prostate = "#E15759",
         Pancreatic = "#B07AA1", Bladder = "#76B7B2", Endometrial = "#EDC948", Ovarian = "#FF9DA7",
         Esophagogastric = "#9C755F", Hepatobiliary = "#8CD17D", Renal = "#86BCB6",
         `Head & Neck` = "#BAB0AC")                       # the lineage palette
COL <- c("all CpGs tested" = "#BDBDBD", "hyper-DMR" = "#C0504D")  # hyper-DMR = the WGD+ red of Fig. 5i
INK <- "#222222"; AX <- 0.35
YLIM <- c(0, 64)                      # data 8-52 % + headroom for the two stat lines
SYMFAM <- cairoSymbolFont(FONT, usePUA = FALSE)

tab <- fread(file.path(DAT, "fig5j_prc2_cgi_share.tsv"))
st <- read.delim(file.path(DAT, "fig5j_panel_stats.tsv"))      # the paired Wilcoxon P, the number of lineages with the share higher in hyper-DMRs, the DMR cohort size
# per lineage: the share of tested CpGs on PRC2-target CGIs among all tested CpGs (the array background) and among the
# CpGs inside the lineage's hyper-DMRs
tab[, `:=`(bg = pct_all_tested, dmr = pct_in_hyper_dmr)]
stopifnot(nrow(tab) == 12)
long <- rbind(tab[, .(lineage, grp = "all CpGs tested", y = bg)],
              tab[, .(lineage, grp = "hyper-DMR",  y = dmr)])
long[, grp := factor(grp, levels = names(COL))]
long[, x := as.numeric(grp)]
seg <- tab[, .(lineage, x = 1, xend = 2, y = bg, yend = dmr)]
lab <- tab[, .(lineage, x = 2, y = dmr)]
p <- ggplot() +
  geom_boxplot(data = long, aes(x = x, y = y, group = grp, fill = grp),
               width = 0.34, linewidth = AX, colour = INK, outlier.shape = NA, alpha = 0.55) +
  geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend, colour = lineage),
               linewidth = 0.42, alpha = 0.85) +
  geom_point(data = long, aes(x = x, y = y, fill = lineage),
             shape = 21, size = 2.1, stroke = 0.35, colour = INK) +
  annotate("text", x = 0.62, y = YLIM[2] - 1.5, hjust = 0, family = FONT, size = 2.4,
           colour = INK, parse = TRUE,
           label = sprintf("paste(\"paired Wilcoxon \", italic(P), \" = %s\")",
                           formatC(st$wilcoxon_p, format = "g", digits = 2))) +
  annotate("text", x = 0.62, y = YLIM[2] - 5.5, hjust = 0, family = FONT, size = 2.4,
           colour = INK, label = sprintf("higher in %d/%d lineages", st$n_lineages_up, nrow(tab))) +
  annotate("text", x = 0.62, y = 1.8, hjust = 0, family = FONT, size = 2.0, colour = "#777777",
           label = sprintf("Hetzel 2025 PRC2 CGIs · %s tumours · %d cohort blocks", format(st$n_tumours, big.mark = ","), st$n_blocks)) +
  scale_fill_manual(values = c(PAL, COL), guide = "none") +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_y_continuous(limits = YLIM, breaks = seq(0, 60, 20), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6, base_family = FONT) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = AX, colour = INK),
        axis.ticks = element_line(linewidth = AX, colour = INK),
        axis.ticks.length = unit(1.1, "pt"),
        axis.text = element_text(size = 6, colour = INK),
        axis.text.x = element_text(size = 6.5, colour = INK),
        axis.title = element_text(size = 6.5, colour = INK),
        plot.margin = margin(2, 2, 1, 1)) +
  ggrepel::geom_text_repel(
    data = lab, aes(x = x, y = y, label = lineage, colour = lineage),
    family = FONT, size = 2.1, direction = "y", hjust = 0, nudge_x = 0.30,
    segment.size = 0.25, segment.colour = "#B0B0B0", box.padding = 0.08,
    min.segment.length = 0, max.overlaps = Inf, seed = 1) +
  scale_x_continuous(breaks = 1:2, labels = names(COL), limits = c(0.55, 3.15), expand = c(0, 0)) +
  labs(y = "CpGs on PRC2-target CGIs (%)")

ggsave(file.path(OUT, "fig5j.png"), p, width = 62, height = 72, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig5j.pdf"), p, width = 62, height = 72, units = "mm", device = cairo_pdf, bg = "white", symbolfamily = SYMFAM)
