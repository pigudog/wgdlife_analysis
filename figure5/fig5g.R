# Fig. 5g: the WGD effect on eight immune-related expression scores in four systems ordered by the immune pressure the cancer
# cells have been under: 2D cell lines and organoids (no adaptive immune system), xenografts in NSG mice (day 14) and tumours
# in immunocompetent A/J mice (days 8 and 14, shaded). Points = estimate with its 95% CI, red = nominal P < 0.05, the 2 v 2
# day-8 intervals in light grey. The culture and mouse columns are on different x scales, each column carrying its own ticks.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid); library(gtable) })
DAT <- "figure5/data"; OUT <- "results/figure5"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
TXT <- 7; TICK <- 6; LW <- 0.35
SIG  <- "#C0504D"   # red        — nominally significant
NS   <- "#9E9E9E"   # grey       — not significant
FAINT<- "#D9D9D9"   # light grey — the 2 v 2 day-8 intervals

G <- fread(file.path(DAT, "fig5g_immune_pressure.tsv"))
A <- fread(file.path(DAT, "fig5g_arm_labels.tsv"))                         # strip labels: the system over its sample sizes, in panel order
ARMLAB <- setNames(paste0(A$label_line1, "\n", A$label_line2), A$arm)
G[, armf := factor(ARMLAB[arm], levels = unname(ARMLAB))]
G[, sig := factor(fifelse(p < 0.05, "nominal P < 0.05", "not significant"),
                  levels = c("nominal P < 0.05", "not significant"))]
G[, wide := arm == "Mouse, A/J d8"]                          # the 2 v 2 time control
G[, y := max(order_ep) + 1 - order_ep]
G[, culture := arm %in% c("2D cell line", "Organoid")]
stopifnot(!any(is.na(G$est)), nrow(G) == 40L)
ROW <- unique(G[, .(y, label)])[order(-y)]

## the two x ranges, forced per column group with an invisible geom_blank
XC <- c(-0.30, 0.55); XM <- c(-1.42, 1.02)
lims <- rbind(G[culture == TRUE,  .(armf = unique(armf)), ][, .(armf, lo = XC[1], hi = XC[2])],
              G[culture == FALSE, .(armf = unique(armf)), ][, .(armf, lo = XM[1], hi = XM[2])])
lims <- melt(lims, id.vars = "armf", value.name = "est")[, y := 1]

g <- ggplot(G, aes(est, y)) +
  ## immunocompetent arms: shaded panel background
  geom_rect(data = G[arm %in% c("Mouse, A/J d8", "Mouse, A/J d14"), .(armf = unique(armf))],
            inherit.aes = FALSE, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            fill = "grey94") +
  geom_blank(data = lims) +
  geom_vline(xintercept = 0, linewidth = LW, colour = "grey45") +
  ## 95% CI — the 2 v 2 day-8 arm at a lighter weight, then everything else
  geom_linerange(data = G[wide == TRUE], aes(xmin = lo, xmax = hi),
                 linewidth = 0.28, colour = FAINT) +
  geom_linerange(data = G[wide == FALSE], aes(xmin = lo, xmax = hi, colour = sig),
                 linewidth = 0.45) +
  geom_point(aes(fill = sig), shape = 21, size = 1.35, stroke = 0.3, colour = "grey20") +
  scale_colour_manual(values = c(`nominal P < 0.05` = SIG, `not significant` = NS),
                      guide = "none") +
  scale_fill_manual(values = c(`nominal P < 0.05` = SIG, `not significant` = NS), name = NULL) +
  scale_y_continuous(breaks = ROW$y, labels = ROW$label,
                     limits = c(0.4, max(G$y) + 0.6), expand = expansion(0)) +
  ## wide (mouse) panels get -1 / 0 ; narrow (culture) panels get their own finer ticks
  scale_x_continuous(expand = expansion(0),
                     breaks = function(l) if (diff(l) > 1.5) c(-1, 0, 1) else c(-0.2, 0.2)) +
  facet_grid(. ~ armf, scales = "free_x") +
  labs(x = "WGD effect (SD)") +
  theme_classic(base_family = FONT, base_size = TICK) +
  theme(axis.line = element_line(linewidth = LW),
        axis.ticks = element_line(linewidth = LW), axis.ticks.length = unit(0.6, "mm"),
        axis.text = element_text(size = TICK, colour = "grey10"),
        axis.text.y = element_text(hjust = 1),
        axis.title.x = element_text(size = TXT, margin = margin(t = 1)),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = TICK, colour = "grey10", lineheight = 1.05,
                                  margin = margin(b = 1.2)),
        panel.spacing.x = unit(1.4, "mm"),
        legend.text = element_text(size = TICK), legend.title = element_blank(),
        legend.key.size = unit(2.2, "mm"), legend.key.spacing.x = unit(1.0, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.4, "mm"),
        legend.position = "bottom", legend.justification = "center",
        plot.margin = margin(9, 2, 1, 1.5))

## ---- regime row above the strips + the scale note under the axis ----------------------
add_furniture <- function(plot) {
  pdf(NULL); on.exit(dev.off(), add = TRUE)      # size the grobs without leaving Rplots.pdf; the text extents baked into the
  gt <- suppressWarnings(ggplotGrob(plot))      # gtable are the null device's (it warns that Arial is not in its font database)
  pan <- gt$layout[grepl("^panel", gt$layout$name), ]; pan <- pan[order(pan$l), ]
  stopifnot(nrow(pan) == 5L)
  ## the two culture columns are on a ~3x finer x scale, so give them 1.22x the width —
  ## otherwise their tick labels collide and the estimates still cannot be read
  gt$widths[pan$l[1:2]] <- gt$widths[pan$l[1:2]] * 1.22
  ## (a) the two regimes, spanning whole groups of columns
  top <- min(gt$layout$t[grepl("^strip", gt$layout$name)])
  gt <- gtable_add_rows(gt, unit(3.1, "mm"), pos = top - 1)
  for (k in list(list(cols = 1:3, lab = "no adaptive immune system"),
                 list(cols = 4:5, lab = "immunocompetent host"))) {
    gr <- grobTree(
      linesGrob(x = unit(c(0.02, 0.98), "npc"), y = unit(c(0.06, 0.06), "npc"),
                gp = gpar(col = "grey55", lwd = LW / 0.75)),
      textGrob(k$lab, y = unit(0.30, "npc"), vjust = 0,
               gp = gpar(fontfamily = FONT, fontsize = TICK, col = "grey25")))
    gt <- gtable_add_grob(gt, gr, t = top, l = pan$l[min(k$cols)], r = pan$l[max(k$cols)],
                          name = paste0("regime", k$cols[1]), clip = "off")
  }
  ## (b) say on the panel that the culture and mouse columns are on different x scales
  ax <- max(gt$layout$b[grepl("^xlab", gt$layout$name)])
  gt <- gtable_add_rows(gt, unit(3.0, "mm"), pos = ax)
  note <- textGrob("culture and mouse columns are on different x scales — read the ticks",
                   y = unit(0.62, "npc"), vjust = 1,
                   gp = gpar(fontfamily = FONT, fontsize = TICK, col = "grey35"))
  gt <- gtable_add_grob(gt, note, t = ax + 1, l = pan$l[1], r = pan$l[5],
                        name = "scalenote", clip = "off")
  gt
}
gt <- add_furniture(g)

W <- 90; H <- 74
ggsave(file.path(OUT, "fig5g.png"), gt, width = W, height = H, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig5g.pdf"), gt, width = W, height = H, units = "mm", device = cairo_pdf, bg = "white")
