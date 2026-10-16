# Fig. 3f: timing of chromosome-arm gains relative to WGD across arm x cancer-type groups of WGD+ PCAWG carcinomas. Two
# caterpillars: F = probability that the first gain preceded WGD, L = probability that the last gain followed WGD, one point
# per group with its 95% CI, at the position of the group in the value-sorted caterpillar (all from the table); filled = arms
# associated with WGD in Fig. 2. Line and band: random-effects mean
# and +/- tau (between-group s.d. after removing sampling error), read from the result table.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
DAT <- "figure3/data"; OUT <- "results/figure3"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
C <- read.delim(file.path(DAT, "fig3f_gain_timing_cells.tsv"), stringsAsFactors = FALSE)
re <- read.delim(file.path(DAT, "fig3f_random_effects.tsv"), stringsAsFactors = FALSE)   # random-effects mean and tau of F and L
st <- lapply(setNames(c("F", "L"), c("F", "L")), function(k) c(mu = re$mu[re$quantity == k], tau = re$tau[re$quantity == k]))
cat_plot <- function(y, lo, hi, i, s, lab, col, xlab = NULL) {
  d <- data.frame(y = y, lo = lo, hi = hi, i = i, fig2 = C$fig2_gain); d <- d[order(d$i), ]   # drawn left to right
  ggplot(d) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = s["mu"] - s["tau"], ymax = s["mu"] + s["tau"], fill = col, alpha = 0.18) +
    geom_hline(yintercept = s["mu"], colour = col, linewidth = 0.5) +
    geom_hline(yintercept = 0.5, colour = "grey55", linewidth = 0.35, linetype = "22") +
    geom_linerange(aes(x = i, ymin = lo, ymax = hi), colour = "grey60", linewidth = 0.3) +
    geom_point(aes(x = i, y = y, fill = fig2), shape = 21, size = 1.15, stroke = 0.3, colour = "grey25") +
    scale_fill_manual(values = c(`TRUE` = "grey15", `FALSE` = "white"), labels = c(`TRUE` = "WGD-associated (Fig. 2)", `FALSE` = "other arms"), name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0, 0), name = lab) +
    scale_x_continuous(expand = expansion(add = 1.5), name = xlab) +
    theme_classic(base_family = FONT, base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          axis.line.y = element_line(linewidth = 0.45), axis.ticks.y = element_line(linewidth = 0.45), axis.text.y = element_text(size = 6.5),
          axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 6.5), legend.position = "none", plot.margin = margin(2, 3, 2, 1))
}
pF <- cat_plot(C$F, C$F_lo, C$F_hi, C$x_F, st$F, "first gain before WGD", "#6E8FB3") + theme(plot.margin = margin(2, 3, 4, 1))
pL <- cat_plot(C$L, C$L_lo, C$L_hi, C$x_L, st$L, "last gain after WGD", "#C9A66B", xlab = sprintf("%d arm × cancer-type cells, sorted", nrow(C))) +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5), legend.key.size = unit(2.4, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.5, "mm"))
g <- pF / pL + plot_layout(heights = c(1, 1.12))
ggsave(file.path(OUT, "fig3f.png"), g, width = 64, height = 68, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig3f.pdf"), g, width = 64, height = 68, units = "mm", device = cairo_pdf, bg = "white")
