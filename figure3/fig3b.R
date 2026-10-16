# Fig. 3b: onset of the HRD-associated mutational process (SBS3) relative to WGD in PCAWG WGD+ tumours with SBS3 >= 10% of
# SNVs. One point per tumour: median onset minus the WGD time in mutation time, grouped by cancer type (types with fewer than
# five tumours pooled as Other); blue = the bootstrap interval of the onset lies before WGD, grey = it spans WGD, gold = after.
# The boxes are drawn from the per-type box statistics of the second table.
suppressPackageStartupMessages({ library(ggplot2) })
DAT <- "figure3/data"; OUT <- "results/figure3"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
H <- read.delim(file.path(DAT, "fig3b_hrd_onset.tsv"), stringsAsFactors = FALSE)       # one row per tumour: onset minus WGD time, side of its interval
B <- read.delim(file.path(DAT, "fig3b_cancer_types.tsv"), stringsAsFactors = FALSE)    # per cancer type, in panel order: n and the box statistics
lev <- paste0(B$cancer_type, "\n", B$n)                                                  # axis label: the type over its n
B$tlab <- factor(lev, levels = lev); H$tlab <- factor(lev[match(H$cancer_type, B$cancer_type)], levels = lev)
g <- ggplot(H, aes(x = tlab, y = onset_minus_wgd)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#6E8FB3", alpha = 0.07) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4, linetype = "22") +
  geom_boxplot(data = B, aes(x = tlab, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax), stat = "identity",
               inherit.aes = FALSE, width = 0.62, fill = "white", colour = "grey25", linewidth = 0.35) +
  geom_point(aes(colour = side), size = 0.7, alpha = 0.85, position = position_jitter(width = 0.16, height = 0, seed = 1)) +
  scale_colour_manual(values = c(pre = "#6E8FB3", spans = "grey55", post = "#C9A66B"), breaks = c("pre", "spans", "post"),
                      labels = c(pre = "before WGD", spans = "spans WGD", post = "after WGD"), name = NULL) +
  scale_y_continuous(breaks = c(-1, -0.5, 0), labels = c("-1", "-0.5", "0"), name = "HRD onset − WGD (mutation time)") +
  coord_cartesian(ylim = c(-1.2, 0.25)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.5, alpha = 1))) +
  theme_classic(base_family = FONT, base_size = 9) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 5.8, colour = "grey10", lineheight = 0.85),
        axis.line = element_line(linewidth = 0.45), axis.ticks = element_line(linewidth = 0.45), axis.text.y = element_text(size = 6.5),
        axis.title.y = element_text(size = 6.6), legend.position = "top", legend.text = element_text(size = 5.5),
        legend.key.size = unit(2.4, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.5, "mm"), plot.margin = margin(2, 3, 2, 2))
ggsave(file.path(OUT, "fig3b.png"), g, width = 58, height = 58, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig3b.pdf"), g, width = 58, height = 58, units = "mm", device = cairo_pdf, bg = "white")
