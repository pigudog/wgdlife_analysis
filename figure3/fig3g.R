# Fig. 3g: convergence of the event spectrum after WGD. Top: schematic (no data) of what is measured: each lineage = a line
# (its centroid) with a band (its within-lineage spread) around the dotted pan-cancer centre; the double arrows are the plotted
# distance before and after WGD. Bottom left (MSK-IMPACT): per cancer type, distance of the type's event-spectrum centroid from the
# pan-cancer centroid scaled by the within-type spread, before vs after WGD (grey = cancer type, black = median across types).
# Bottom right: cancer-type explanatory power of the event spectrum (PERMANOVA R^2), MSK-IMPACT (circle) and TCGA (triangle).
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
DAT <- "figure3/data"; OUT <- "results/figure3"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"; WGD_COL <- "grey45"
base_theme <- function() theme_classic(base_family = FONT, base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.45), axis.ticks = element_line(linewidth = 0.45),
        axis.text = element_text(size = 6.3, colour = "grey10"), axis.title = element_text(size = 6.4), plot.margin = margin(2, 2, 1, 1))

# top: schematic
xs <- seq(0, 1, length.out = 200); wgd_x <- 0.42; CEN <- 0.5; HW <- 0.075
traj <- function(off) { s <- 1 / (1 + exp(-15 * (xs - (wgd_x + 0.14)))); CEN + off * (1 - 0.75 * s) }
offs <- c(A = 0.34, B = 0.10, C = -0.28)
lin <- do.call(rbind, lapply(names(offs), function(l) data.frame(lineage = l, x = xs, y = traj(offs[[l]]))))
greys <- c(A = "grey10", B = "grey40", C = "grey62")
ya <- function(l, x) lin$y[lin$lineage == l][which.min(abs(xs - x))]
x_pre <- 0.2; x_post <- 0.9
pt <- ggplot(lin, aes(x, y)) +
  geom_ribbon(aes(ymin = y - HW, ymax = y + HW, fill = lineage, group = lineage), alpha = 0.18, colour = NA) +
  annotate("segment", x = 0, xend = 1, y = CEN, yend = CEN, colour = "grey55", linewidth = 0.4, linetype = "13") +
  annotate("segment", x = wgd_x, xend = wgd_x, y = 0.02, yend = 0.93, linetype = "22", colour = WGD_COL, linewidth = 0.45) +
  geom_path(aes(colour = lineage, group = lineage), linewidth = 0.8, lineend = "round") +
  annotate("segment", x = c(x_pre, x_post), xend = c(x_pre, x_post), y = CEN, yend = c(ya("A", x_pre), ya("A", x_post)) - 0.012,
           colour = "grey20", linewidth = 0.35, arrow = arrow(ends = "both", length = unit(0.8, "mm"), type = "closed")) +
  annotate("text", x = -0.015, y = CEN + offs, label = names(offs), colour = greys, family = FONT, fontface = "bold", size = 2.2, hjust = 1) +
  annotate("text", x = wgd_x, y = 0.98, label = "WGD", colour = WGD_COL, family = FONT, fontface = "bold", size = 2.1, vjust = 0) +
  annotate("text", x = c(wgd_x / 2, (1 + wgd_x) / 2), y = 0.06, label = c("pre-WGD", "post-WGD"), colour = c("#6E8FB3", "#C9A66B"),
           family = FONT, fontface = "bold", size = 2.0) +
  scale_colour_manual(values = greys, guide = "none") + scale_fill_manual(values = greys, guide = "none") +
  scale_x_continuous(expand = c(0.01, 0.01)) + coord_cartesian(xlim = c(-0.05, 1), ylim = c(0.04, 1.02), clip = "off") +
  theme_void(base_family = FONT) + theme(plot.margin = margin(3, 3, 1, 6))

# bottom left: lineage distance from the pan-cancer centre
bt <- read.delim(file.path(DAT, "fig3g_convergence_by_type.tsv"))
D <- rbind(data.frame(ct = bt$cancer_type, x = 1, v = bt$ratio_between_within_pre), data.frame(ct = bt$cancer_type, x = 2, v = bt$ratio_between_within_post))
md <- read.delim(file.path(DAT, "fig3g_medians.tsv"))                                  # the black line: median across cancer types
Md <- data.frame(x = c(1, 2), v = md$median_distinctiveness[match(c("pre", "post"), md$epoch)])
xsc <- scale_x_continuous(breaks = c(1, 2), labels = c("pre", "post"), limits = c(0.75, 2.25))
slope_theme <- function() base_theme() + theme(axis.title.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())
p1 <- ggplot() +
  geom_line(data = D, aes(x, v, group = ct), colour = "grey70", linewidth = 0.3) +
  geom_line(data = Md, aes(x, v), colour = "grey10", linewidth = 0.7) + geom_point(data = Md, aes(x, v), colour = "grey10", size = 1.4) +
  xsc + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
  labs(y = "lineage distinctiveness") + slope_theme()

# bottom right: cancer-type R^2
cm <- read.delim(file.path(DAT, "fig3g_cancer_type_r2.tsv"))                          # one row per platform: R^2 before and after WGD
r2 <- rbind(data.frame(pf = toupper(cm$platform), x = 1, v = 100 * cm$pre), data.frame(pf = toupper(cm$platform), x = 2, v = 100 * cm$post))
p2 <- ggplot(r2, aes(x, v, group = pf)) +
  geom_line(colour = "grey25", linewidth = 0.55) +
  geom_point(aes(shape = pf), size = 1.5, fill = "white", colour = "grey15", stroke = 0.45) +
  scale_shape_manual(values = c(MSK = 21, TCGA = 24), name = NULL) + xsc +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(y = "cancer-type R² (%)") + slope_theme() +
  theme(legend.position = c(0.72, 0.9), legend.text = element_text(size = 5.3), legend.key.size = unit(2, "mm"),
        legend.background = element_blank(), legend.spacing.y = unit(0, "mm"))
p <- pt / (p1 | p2) + plot_layout(heights = c(0.42, 0.58))
ggsave(file.path(OUT, "fig3g.png"), p, width = 70, height = 68, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig3g.pdf"), p, width = 70, height = 68, units = "mm", device = cairo_pdf, bg = "white")
