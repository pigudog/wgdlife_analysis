# Fig. 5k: the Hallmark programmes at which WGD+ tumours gain (right, red) or lose (left, blue) methylation: over-representation
# of each of the 50 Hallmark gene sets in the per-lineage WGD hyper-DMRs and hypo-DMRs, combined across the 12 lineages
# (Stouffer Z of the one-sided over-representation P, signed + for methylation gained and − for lost). Every programme at
# BH FDR < 0.05 on its side is drawn; bar colour = signed Z; the immune / IFN programmes of Fig. 5a in bold accent.
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
`%|%` <- function(a, b) ifelse(is.na(a), b, a)
DAT <- "figure5/data"; OUT <- "results/figure5"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
ACC <- "#B5533C"; INK <- "#222222"; AX <- 0.35
BLUE <- "#8AA9C9"; RED <- "#C0504D"; MID <- "#F3F1EF"     # the WGD− / WGD+ pair of Fig. 5i
# tint toward white by fraction f -- the light end of each side
mix  <- function(col, f) { m <- grDevices::col2rgb(col) / 255; grDevices::rgb(t(m * (1 - f) + f)) }
BLUE_LT <- mix(BLUE, 0.45); RED_LT <- mix(RED, 0.45)
LAB  <- 2.15                                     # mm, ~6.1 pt

H <- fread(file.path(DAT, "fig5k_hallmark_methylation.tsv"))        # the programmes at FDR < 0.05, with their side
H[, x := fifelse(side == "gained", Z, -Z)]
SHORT <- c("Allograft Rejection" = "Allograft rejection", "Interferon Gamma Response" = "IFN-γ response",
           "Interferon Alpha Response" = "IFN-α response", "Pancreas Beta Cells" = "Pancreas β cells",
           "KRAS Signaling Dn" = "KRAS signaling (down)", "KRAS Signaling Up" = "KRAS signaling (up)",
           "Epithelial Mesenchymal Transition" = "EMT", "Apical Junction" = "Apical junction")
H[, lab := fifelse(set %in% names(SHORT), SHORT[set], set)]
H[, lab := sub("^Estrogen Response (Early|Late)$", "Estrogen response \\L\\1", lab, perl = TRUE)]
stopifnot(!anyDuplicated(H$set))                   # a programme never sits on both sides
setorder(H, -x)
H[, y := as.numeric(rev(seq_len(.N)))]           # numeric: the block gap below adds 0.45

# narrow layout: programme names sit in a left column (immune in bold accent), bars diverge from zero to their right
H[, lab := c("Estrogen response late" = "Estrogen resp. late",
             "Estrogen response early" = "Estrogen resp. early",
             # descriptive names: MSigDB's "KRAS signaling UP / DN" are the genes induced / repressed when oncogenic
             # KRAS is switched on, not a KRAS-activity readout here
             "KRAS signaling (down)" = "KRAS-repressed genes",
             "KRAS signaling (up)" = "KRAS-induced genes")[lab] %|% lab]
XL <- c(-10.4, 9.0)
GAP <- 0.45                                      # extra space between the gained and the lost block
H[side == "gained", y := y + GAP]
YTOP <- max(H$y) + 0.5                           # top edge of the bar area
p <- ggplot(H) +
  # zero line over the bar area only (not through the header); bars without outlines
  annotate("segment", x = 0, xend = 0, y = 0.4, yend = YTOP, linewidth = AX, colour = INK) +
  geom_col(aes(x = x, y = y, fill = x), orientation = "y", width = 0.7, colour = NA) +
  geom_text(aes(x = XL[1] - 0.6, y = y, label = lab, colour = immune,
                fontface = ifelse(immune, "bold", "plain")),
            hjust = 1, family = FONT, size = LAB) +
  annotate("text", x = -0.6, y = YTOP + 0.55, hjust = 1, family = FONT, size = LAB,
           colour = "grey45", label = "lost") +
  annotate("text", x = 0.6, y = YTOP + 0.55, hjust = 0, family = FONT, size = LAB,
           colour = "grey45", label = "gained") +
  # each side runs from its half-saturated tint at |Z| = 2.5 (just below the weakest bar, 2.9) to the
  # full colour at the strongest bar; colour still orders bars by strength
  scale_fill_gradientn(colours = c(BLUE, BLUE_LT, MID, RED_LT, RED),
                       values = scales::rescale(c(-1, -2.5 / max(abs(H$x)), 0,
                                                  2.5 / max(abs(H$x)), 1)),
                       limits = c(-1, 1) * max(abs(H$x)), guide = "none") +
  scale_colour_manual(values = c(`TRUE` = ACC, `FALSE` = INK), guide = "none") +
  scale_x_continuous(breaks = c(-8, 0, 8), name = "enrichment Z") +
  # ranges set on the coordinate system, not the scales, so the name column left of XL[1] is kept
  coord_cartesian(xlim = XL, ylim = c(0.4, YTOP + 1.1), expand = FALSE, clip = "off") +
  labs(y = NULL) +
  theme_minimal(base_size = 6, base_family = FONT) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = AX, colour = INK),
        axis.ticks.x = element_line(linewidth = AX, colour = INK),
        axis.ticks.length = unit(1.1, "pt"),
        axis.text.x = element_text(size = 6, colour = INK),
        axis.title.x = element_text(size = 6.5, colour = INK),
        plot.margin = margin(2, 2, 1, 66))       # pt; the left column holds the programme names

ggsave(file.path(OUT, "fig5k.png"), p, width = 45, height = 55, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig5k.pdf"), p, width = 45, height = 55, units = "mm", device = cairo_pdf, bg = "white")
