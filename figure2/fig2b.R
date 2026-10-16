# Fig. 2b: overlap of the WGD-associated alterations retained pan-cancer and within individual cancer types.
# Retained = FDR <= 0.10 in both discovery cohorts with the same direction, or in one cohort with P < 0.05 in the other.
suppressPackageStartupMessages(library(ggplot2))
DAT <- "figure2/data"; OUT <- "results/figure2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
v <- read.delim(file.path(DAT, "fig2b_retained_overlap.tsv"))
circ <- function(x0, r = 1, n = 300) data.frame(x = x0 + r * cos(seq(0, 2 * pi, length.out = n)), y = r * sin(seq(0, 2 * pi, length.out = n)))
L <- circ(-0.55); R <- circ(0.55)
lens <- rbind(L[(L$x - 0.55)^2 + L$y^2 <= 1, ], R[(R$x + 0.55)^2 + R$y^2 <= 1, ]); lens <- lens[order(atan2(lens$y, lens$x)), ]
p <- ggplot() +
  geom_polygon(data = L, aes(x, y), fill = "#A6A6A6", colour = NA) +
  geom_polygon(data = lens, aes(x, y), fill = "#2C6EA6", colour = NA) +
  geom_path(data = rbind(L, L[1, ]), aes(x, y), colour = "grey30", linewidth = 0.45) +
  geom_path(data = rbind(R, R[1, ]), aes(x, y), colour = "grey30", linewidth = 0.45) +
  annotate("text", x = -0.95, y = 0, label = v$pancancer_only, family = FONT, size = 3.4, colour = "grey10") +
  annotate("text", x = 0, y = 0, label = v$retained_both, family = FONT, size = 3.4, colour = "white", fontface = "bold") +
  annotate("text", x = 0.95, y = 0, label = v$type_only, family = FONT, size = 3.4, colour = "grey10") +
  annotate("text", x = -0.85, y = -1.3, label = sprintf("Pan-cancer\n(%d)", v$retained_pancancer), family = FONT, size = 3.0, lineheight = 0.95) +
  annotate("text", x = 0.85, y = -1.3, label = sprintf("Cancer-type-\nspecific (%d)", v$retained_type_specific), family = FONT, size = 3.0, lineheight = 0.95) +
  coord_equal(xlim = c(-1.7, 1.7), ylim = c(-1.65, 1.1)) + theme_void(base_family = FONT)
ggsave(file.path(OUT, "fig2b.png"), p, width = 62, height = 50, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig2b.pdf"), p, width = 62, height = 50, units = "mm", device = cairo_pdf, bg = "white")
