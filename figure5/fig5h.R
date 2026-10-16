# Fig. 5h: immune-escape events in WGD+ versus WGD− tumours, pooled (All) and per carcinoma lineage (lineages with at least
# 100 tumours in the smaller WGD arm). Column 1: HLA class I LOH observed / expected from each tumour's genome-wide LOH, as
# log2 of the WGD+ over the WGD− O/E; columns 2-7: truncating mutations in six CRISPR immune-escape gene sets, log2 of the
# WGD+ over the WGD− prevalence. Hatched = fewer than 10 altered patients in an arm; stars = BH q < 0.05 / 0.01 / 0.001
# (never on hatched cells).
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggtext); library(grid); library(cowplot) })
cowplot::set_null_device("agg")   # the default pdf() null device cannot resolve Arial
DAT <- "figure5/data"; OUT <- "results/figure5"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
MINUS <- "−"                 # U+2212; U+207B (superscript minus) is absent from Arial
LW    <- 0.35 / 2.13              # 0.35 pt rules in ggplot linewidth units
LIM   <- 1.5                      # colour scale saturates here (|log2| 1.5 = 2.83-fold)

d <- rbind(fread(file.path(DAT, "fig5h_escape_by_lineage.tsv")),
           fread(file.path(DAT, "fig5h_escape_pooled.tsv")))
tst <- fread(file.path(DAT, "fig5h_escape_tests.tsv"))

ROWS <- c("HLA class I LOH", "T cell killing", "MHC-I surface", "in vivo immune pressure",
          "γδ T cell killing", "NK cell killing", "macrophage")
# Display names. The HLA row is an O/E against each tumour's own genome-wide LOH, so its label says so. The six SNV rows
# are all "genes whose knockout lets the cell escape <pressure>", so they are named by the escape: a gene whose knockout
# LOWERS surface MHC-I is a "positive regulator of MHC-I", one whose knockout confers survival under T-cell attack gives
# "resistance to T-cell killing". The in vivo axis is the only one whose read-out comes from a real immune environment
# (knockout libraries grown in immunocompetent mice), so it is named for immune control.
RDISP <- c("HLA-I LOH, O/E",
           "Resistance to T-cell killing", "Positive regulators of MHC-I",
           "Resistance to immune control, *in vivo*", "Resistance to γδ T-cell killing",
           "Resistance to NK-cell killing", "Resistance to macrophages")
names(RDISP) <- ROWS
# columns: pooled first, then the lineages by WGD+ cohort size, largest first
LIN <- c("All", "NSCLC", "Breast", "Esophagogastric", "Ovarian", "Colorectal", "Bladder",
         "Head&Neck", "Endometrial", "Hepatobiliary", "Renal")
d <- d[lineage %in% LIN]
stopifnot(setequal(d$row, ROWS), setequal(d$lineage, LIN))

d[, row := factor(row, levels = rev(ROWS))]
d[, lineage := factor(lineage, levels = LIN)]
d[, block := factor(fifelse(event_class == "loh", "loh", "snv"), levels = c("loh", "snv"),
                    labels = c("/ genome<br>LOH", "prev."))]
d[, colgrp := factor(fifelse(lineage == "All", "pooled", "bylin"),
                     levels = c("pooled", "bylin"))]

# ------------------------------------------------------------------- fill: log2 ratio -> hex
ramp <- colorRamp(c("#2F4B6E", "#FFFFFF", "#C0504D"), space = "Lab")   # muted navy / white / red
hex  <- function(z) {
  z <- pmax(pmin(z, LIM), -LIM)
  apply(ramp((z + LIM) / (2 * LIM)), 1, function(v) rgb(v[1], v[2], v[3], maxColorValue = 255))
}
wash <- function(col, a) {                    # blend toward white: thin cells recede
  m <- col2rgb(col); apply(m, 2, function(v)
    rgb(v[1] + (255 - v[1]) * a, v[2] + (255 - v[2]) * a, v[3] + (255 - v[3]) * a,
        maxColorValue = 255))
}
GREY <- "#EDEDED"                              # ratio undefined (an arm at 0 altered patients)
d[, fill := GREY]                              # fifelse would evaluate hex() on the NA ratios
d[defined == TRUE, fill := hex(log2_ratio)]
d[thin == TRUE & defined == TRUE, fill := wash(fill, 0.45)]
lum <- function(col) { m <- col2rgb(col); as.numeric(0.299 * m[1, ] + 0.587 * m[2, ] + 0.114 * m[3, ]) }
d[, ink := fifelse(lum(fill) < 130, "white", "grey8")]

# ------------------------------------------------------------ hatching over the thin cells
# a unit tile is [-0.5, .5]^2 around (xi, yi); the stroke v = u + b is clipped to u in [0, 1].
# yi is the row's position INSIDE its facet — facet_grid gives each block its own free y scale.
d[, yfac := as.integer(factor(as.character(row),
                              levels = rev(ROWS[if (.BY$block == levels(block)[1]) 1 else 2:7]))),
  by = block]
d[, xfac := as.integer(lineage) - fifelse(colgrp == "pooled", 0L, 1L)]   # x resets per facet
hatch <- d[thin == TRUE][, {
  xi <- xfac[1]; yi <- yfac[1]
  b  <- c(-0.72, -0.36, 0, 0.36, 0.72)
  u0 <- pmax(0, -b); u1 <- pmin(1, 1 - b)
  k <- 0.93                                    # inset, so no stroke touches the white border
  .(x = xi + k * (u0 - 0.5), xend = xi + k * (u1 - 0.5),
    y = yi + k * (u0 + b - 0.5), yend = yi + k * (u1 + b - 0.5))
}, by = .(row, lineage, block, colgrp)]

# significance: BH q of the tested (pooled) cells, drawn as stars; thin/zero cells never starred
q <- tst[, .(row, lineage, fdr_bh)]
d <- merge(d, q, by = c("row", "lineage"), all.x = TRUE, sort = FALSE)
d[, star := fifelse(is.na(fdr_bh) | thin == TRUE | is.na(fill), NA_character_,
             fifelse(fdr_bh < 0.001, "***", fifelse(fdr_bh < 0.01, "**",
             fifelse(fdr_bh < 0.05, "*", NA_character_))))]
d[, row := factor(as.character(row), levels = rev(ROWS))]   # merge() drops the factor levels
d[, lineage := factor(as.character(lineage), levels = LIN)]
setorder(d, -row, lineage)

# --------------------------------------------------------------------------------- the plot
CB <- guide_colourbar(ticks = FALSE, frame.colour = NA,
                      barwidth = unit(16, "mm"), barheight = unit(1.9, "mm"))
pl <- ggplot(d, aes(row, lineage)) +      # transposed: gene sets across, lineages down
  geom_tile(aes(fill = fill), colour = "white", linewidth = 0.5) +
  geom_text(aes(label = star, colour = ink), size = 7 / .pt, family = FONT,
            vjust = 0.72, na.rm = TRUE) +
  scale_colour_identity() +
  scale_fill_identity(guide = "none") +
  scale_x_discrete(expand = c(0, 0), labels = function(x) RDISP[x]) +
  scale_y_discrete(expand = c(0, 0), limits = rev) +
  facet_grid(colgrp ~ block, scales = "free", space = "free") +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = FONT, base_size = 7) +
  theme(
    axis.text.x = element_markdown(size = 6, colour = "grey10", angle = 45, hjust = 1,
                                   vjust = 1, margin = margin(t = 1)),
    axis.text.y = element_text(size = 6, colour = "grey10", hjust = 1,
                               margin = margin(r = 1.6)),
    strip.text.y.right = element_blank(),
    strip.text.x = element_blank(),
    strip.placement = "outside",
    panel.spacing.y = unit(1.7, "mm"), panel.spacing.x = unit(1.8, "mm"),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(1.5, 1, 1.5, 1.5))

# ------------------------------------------------- the colourbar, drawn from a dummy layer
bar <- ggplot(data.table(z = c(-LIM, LIM)), aes(1, 1, fill = z)) + geom_tile() +
  scale_fill_gradientn(
    colours = rgb(ramp(seq(0, 1, length.out = 64)), maxColorValue = 255),
    limits = c(-LIM, LIM), breaks = c(-1.5, -1, 0, 1, 1.5),
    labels = c(paste0(MINUS, "1.5"), paste0(MINUS, "1"), "0", "1", "1.5"),
    name = NULL, guide = CB) +
  theme_void(base_family = FONT) +
  theme(legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 6, colour = "grey10", family = FONT),
        legend.margin = margin(0, 0, 0, 16), legend.key.spacing = unit(0, "mm"))

W <- 58; H <- 72
draw <- function() {
  # extracted on the ACTIVE device: on the default pdf() one Arial is not resolvable
  leg <- cowplot::get_plot_component(bar, "guide-box-top", return_all = TRUE)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit.c(unit(4.6, "mm"),
                                                                   unit(1, "null")))))
  pushViewport(viewport(layout.pos.row = 1)); grid.draw(leg); popViewport()
  pushViewport(viewport(layout.pos.row = 2)); print(pl, newpage = FALSE); popViewport()
}

ragg::agg_png(file.path(OUT, "fig5h.png"), width = W, height = H, units = "mm", res = 600, background = "white"); draw(); invisible(dev.off())
cairo_pdf(file.path(OUT, "fig5h.pdf"), width = W / 25.4, height = H / 25.4, bg = "white"); draw(); invisible(dev.off())
