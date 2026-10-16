# Fig. 5i: global methylation loss with WGD. Per carcinoma lineage, the median methylation over the common-PMD solo-WCGW
# probes in WGD− and WGD+ tumours (the purity-adjusted predictions of the per-lineage model, so the gap between a lineage's
# two points is its adjusted WGD coefficient), one line per lineage, boxes summarising the 12 lineages; paired Wilcoxon
# signed-rank test across lineages.
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(data.table) })
DAT <- "figure5/data"; OUT <- "results/figure5"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
PAL <- c(NSCLC = "#4E79A7", Breast = "#F28E2B", Colorectal = "#59A14F", Prostate = "#E15759",
         Pancreatic = "#B07AA1", Bladder = "#76B7B2", Endometrial = "#EDC948", Ovarian = "#FF9DA7",
         Esophagogastric = "#9C755F", Hepatobiliary = "#8CD17D", Renal = "#86BCB6",
         `Head & Neck` = "#BAB0AC")                       # the lineage palette
COL <- c("WGD−" = "#8AA9C9", "WGD+" = "#C0504D")
INK <- "#222222"; AX <- 0.35
SYMFAM <- cairoSymbolFont(FONT, usePUA = FALSE)

tab <- fread(file.path(DAT, "fig5i_pmd_methylation.tsv"))
st <- read.delim(file.path(DAT, "fig5i_panel_stats.tsv"))      # the paired Wilcoxon P, the number of lineages with WGD+ lower, the tumour and block totals
tab[, `:=`(neg = adj_neg, pos = adj_pos)]

long <- rbind(tab[, .(lineage, grp = "WGD−", y = neg)],
              tab[, .(lineage, grp = "WGD+",      y = pos)])
long[, grp := factor(grp, levels = c("WGD−", "WGD+"))]
long[, x := as.numeric(grp)]
seg <- tab[, .(lineage, x = 1, xend = 2, y = neg, yend = pos)]
lab <- tab[, .(lineage, x = 2, y = pos)]

p <- ggplot() +
  geom_boxplot(data = long, aes(x = x, y = y, group = grp, fill = grp),
               width = 0.34, linewidth = AX, colour = INK, outlier.shape = NA, alpha = 0.55) +
  geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend, colour = lineage),
               linewidth = 0.42, alpha = 0.85) +
  geom_point(data = long, aes(x = x, y = y, fill = lineage),
             shape = 21, size = 2.1, stroke = 0.35, colour = INK) +
  ggrepel::geom_text_repel(
    data = lab, aes(x = x, y = y, label = lineage, colour = lineage),
    family = FONT, size = 2.1, direction = "y", hjust = 0, nudge_x = 0.30,
    segment.size = 0.25, segment.colour = "#B0B0B0", box.padding = 0.08,
    min.segment.length = 0, max.overlaps = Inf, seed = 1) +
  annotate("text", x = 0.62, y = max(long$y) + 4.0, hjust = 0, family = FONT, size = 2.4,
           colour = INK, parse = TRUE,
           label = sprintf("paste(\"paired Wilcoxon \", italic(P), \" = %.3f\")", st$wilcoxon_p)) +
  annotate("text", x = 0.62, y = max(long$y) + 1.6, hjust = 0, family = FONT, size = 2.4,
           colour = INK, label = sprintf("WGD+ lower in %d/%d lineages", st$n_lineages_down, nrow(tab))) +
  annotate("text", x = 0.62, y = min(long$y) - 1.2, hjust = 0, family = FONT, size = 2.0,
           colour = "#777777",
           label = sprintf("%s tumours · %s cohort blocks", format(st$n_tumours, big.mark = ","), st$n_blocks)) +
  scale_fill_manual(values = c(PAL, COL), guide = "none") +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_x_continuous(breaks = 1:2, labels = c("WGD−", "WGD+"),
                     limits = c(0.55, 3.15), expand = c(0, 0)) +
  scale_y_continuous(name = "common-PMD solo-WCGW methylation (%), purity-adjusted",
                     expand = expansion(mult = c(0.04, 0.10))) +
  labs(x = NULL) +
  theme_minimal(base_size = 6, base_family = FONT) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = AX, colour = INK),
        axis.ticks = element_line(linewidth = AX, colour = INK),
        axis.ticks.length = unit(1.1, "pt"),
        axis.text = element_text(size = 6, colour = INK),
        axis.text.x = element_text(size = 6.5, colour = INK),
        axis.title = element_text(size = 6.5, colour = INK),
        plot.margin = margin(2, 2, 1, 1))

ggsave(file.path(OUT, "fig5i.png"), p, width = 62, height = 72, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig5i.pdf"), p, width = 62, height = 72, units = "mm", device = cairo_pdf, bg = "white", symbolfamily = SYMFAM)
