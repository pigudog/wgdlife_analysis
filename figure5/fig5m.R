# Fig. 5m: one point per gene (the measurable interferon, antigen-presentation and innate-sensing genes of the methylation
# screen): the pan-cancer WGD effect on promoter methylation (Δ M value, WGD+ vs WGD−; x) against the expression change of
# the tumours whose promoter is abnormally methylated relative to the others (log2, adjusted for local copy number, purity
# and study; y), coloured by gene group. Least-squares fit with its 95% band; Spearman correlation across the genes and the
# number of genes in the lower-right quadrant (methylation gained in WGD+ and expression down).
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(data.table) })
DAT <- "figure5/data"; OUT <- "results/figure5"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
INK <- "#222222"; AX <- 0.35
GCOL <- c("IFN signalling" = "#F28E2B", "MHC-I / APM" = "#4E79A7", "Innate sensing" = "#59A14F", "MHC-II" = "#B07AA1")
darken <- function(col, f = 0.72) { m <- grDevices::col2rgb(col) * f / 255; grDevices::rgb(t(m)) }
SYMFAM <- cairoSymbolFont(FONT, usePUA = FALSE)

P <- fread(file.path(DAT, "fig5m_methylation_effect.tsv"))
P <- merge(P, fread(file.path(DAT, "fig5m_silencing_effect.tsv"))[, .(gene, est_abn)], by = "gene")
P[, est_RNA := est_abn]                    # y = within-tumour silencing
P[, group := factor(group, names(GCOL))]
st <- read.delim(file.path(DAT, "fig5m_panel_stats.tsv"))      # Spearman rho and P across the genes, genes in the lower-right quadrant
FIT <- fread(file.path(DAT, "fig5m_fit_line.tsv"))   # the least-squares line and its 95% confidence band across the genes
# P in plotmath: 9.2 x 10^-5 when small, two significant digits otherwise
p_txt <- if (st$spearman_p < 0.001) {
  e <- floor(log10(st$spearman_p)); sprintf("%.1f %%*%% 10^%d", st$spearman_p / 10^e, e)
} else sprintf("%.2g", st$spearman_p)

p <- ggplot(P, aes(est_methylation, est_RNA)) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.35, colour = "grey55") +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.35, colour = "grey55") +
  # the least-squares fit across the genes with its 95% band
  geom_ribbon(data = FIT, aes(x = x, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey85", alpha = 0.6, colour = NA) +
  geom_line(data = FIT, aes(x = x, y = y), inherit.aes = FALSE, colour = "grey25", linewidth = 0.45) +
  geom_point(aes(fill = group), shape = 21, size = 1.9, stroke = 0.3, colour = "grey15") +
  geom_text_repel(aes(label = gene, colour = group), family = FONT, size = 1.9, seed = 3, max.overlaps = Inf,
                  box.padding = 0.32, point.padding = 0.18, min.segment.length = 0, segment.size = 0.2,
                  force = 2.2, force_pull = 0.6, max.time = 3, max.iter = 2e4, show.legend = FALSE) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.04, vjust = -2.6, family = FONT, size = 2.1, colour = INK,
           parse = TRUE, label = sprintf("Spearman~rho == %.2f*','~italic(P) == %s", st$spearman_rho, p_txt)) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.04, vjust = -0.9, family = FONT, size = 2.1, colour = INK,
           label = sprintf("gained in WGD+ and silenced: %d/%d", st$n_lower_right, nrow(P))) +
  scale_fill_manual(values = GCOL, name = NULL) +
  scale_colour_manual(values = vapply(GCOL, darken, ""), guide = "none") +
  labs(x = "Δ promoter M value (WGD+ vs WGD−)",
       y = "Expression change when promoter\nmethylated (log2)") +
  theme_classic(base_family = FONT, base_size = 6) +
  theme(axis.line = element_line(linewidth = AX, colour = INK), axis.ticks = element_line(linewidth = AX, colour = INK),
        axis.text = element_text(size = 6, colour = INK), axis.title = element_text(size = 6.5, colour = INK),
        legend.position = "bottom", legend.text = element_text(size = 6), legend.key.size = unit(2.6, "mm"),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(2, 3, 1, 1))

ggsave(file.path(OUT, "fig5m.png"), p, width = 72, height = 66, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig5m.pdf"), p, width = 72, height = 66, units = "mm", device = cairo_pdf, bg = "white", symbolfamily = SYMFAM)
