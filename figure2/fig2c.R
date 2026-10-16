# Fig. 2c: direction of the retained WGD-associated alterations: co-occurrence (CO, green) or mutual exclusivity (ME, purple).
suppressPackageStartupMessages(library(ggplot2))
DAT <- "figure2/data"; OUT <- "results/figure2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
pc <- read.delim(file.path(DAT, "fig2c_direction.tsv"))
pc$direction <- factor(pc$direction, levels = c("CO", "ME"))
lab <- setNames(sprintf("%d %s", pc$n, c(CO = "co-occurrent with WGD", ME = "mutually exclusive with WGD")[as.character(pc$direction)]), pc$direction)
p <- ggplot(pc, aes(x = 1, y = n, fill = direction)) +
  geom_col(width = 1, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_stack(vjust = 0.5), colour = "white", family = FONT, size = 3.1, fontface = "bold") +
  scale_fill_manual(values = c(CO = "#3A9E6F", ME = "#6A4C9C"), labels = lab, name = NULL) +
  coord_polar(theta = "y", start = 0) + xlim(0.5, 1.5) + theme_void(base_family = FONT) +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 7.8), legend.key.size = unit(3.2, "mm"),
        plot.margin = margin(2, 2, 2, 2))
ggsave(file.path(OUT, "fig2c.png"), p, width = 52, height = 62, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig2c.pdf"), p, width = 52, height = 62, units = "mm", device = cairo_pdf, bg = "white")
