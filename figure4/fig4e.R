# Fig. 4e: |difference in WGD-association effect size| between metastasis-site and primary-site samples per locus (observed)
# against the admixed controls (expected); red dashed line = the 99th percentile of the controls (the threshold of the
# category rule); P = two-sided Wilcoxon test. Threshold and P are read from the result table.
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
DAT <- "figure4/data"; OUT <- "results/figure4"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
theme_f4 <- function(bs = 7) theme_classic(base_family = FONT, base_size = bs) +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35), axis.text = element_text(size = bs - 0.8, colour = "grey10"),
        axis.title = element_text(size = bs - 0.3), legend.text = element_text(size = bs - 1.3), legend.title = element_text(size = bs - 1),
        legend.key.size = unit(2.6, "mm"), legend.background = element_blank(), plot.margin = margin(3, 4, 2, 2))
fmtp <- function(p) if (p < 1e-3) sprintf("italic(P) == %s", sub("e-0?", " %*% 10^-", formatC(p, format = "e", digits = 0))) else sprintf("italic(P) == %.2g", p)
SYMFAM <- cairoSymbolFont(FONT, usePUA = FALSE)
d <- fread(file.path(DAT, "fig4_pm_categories.tsv")); nl <- fread(file.path(DAT, "fig4_pm_null.tsv"))
st <- fread(file.path(DAT, "fig4e_stats.tsv"))                                 # the 99th-percentile threshold and the Wilcoxon P
B <- rbind(data.table(g = "Observed", v = abs(d$dist_obs)), data.table(g = "Expected", v = abs(nl$dist))); B[, g := factor(g, c("Expected", "Observed"))]
g1 <- ggplot(B, aes(v, g)) + geom_point(position = position_jitter(height = 0.18, width = 0, seed = 1), size = 0.25, colour = "grey15", alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, width = 0.45, fill = NA, linewidth = 0.3, colour = "grey30") +
  geom_vline(xintercept = st$threshold_99, colour = "#C0392B", linetype = "22", linewidth = 0.35) +
  annotate("text", x = max(B$v), y = 1.5, label = fmtp(st$wilcoxon_p), parse = TRUE, hjust = 1, family = FONT, size = 1.9) +
  scale_x_sqrt(breaks = c(0, 5, 10, 25, 50, 100)) + labs(x = "|ES difference|, metastasis vs primary", y = NULL) + theme_f4() +
  theme(axis.text.y = element_text(size = 5.8), plot.margin = margin(2, 4, 1, 2))
ggsave(file.path(OUT, "fig4e.png"), g1, width = 66, height = 26, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig4e.pdf"), g1, width = 66, height = 26, units = "mm", device = cairo_pdf, bg = "white", symbolfamily = SYMFAM)
