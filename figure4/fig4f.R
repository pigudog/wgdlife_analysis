# Fig. 4f: mean WGD-association effect size per locus in primary-site (x) and metastasis-site (y) samples. Categories are
# nested inside the distance gate |ES_M - ES_P| / sqrt(2) > 99th percentile of the admixed controls: specific = significant in
# >= 6 of 10 runs on one side and <= 2 on the other (the other side testable in >= 6 runs), enriched = the rest above the gate
# (side by magnitude), switch = both sides >= 6 with opposite sign; shared = below the gate. Red = metastasis (filled specific,
# open enriched), blue = primary, gold = switch, grey = shared. Labels: the loci farthest from the diagonal, six outside the
# zoom window and four inside it, plus the loci named in the text; the yellow square is the zoomed region (inset).
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(patchwork); library(ggrepel) })
DAT <- "figure4/data"; OUT <- "results/figure4"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
theme_f4 <- function(bs = 7) theme_classic(base_family = FONT, base_size = bs) +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35), axis.text = element_text(size = bs - 0.8, colour = "grey10"),
        axis.title = element_text(size = bs - 0.3), legend.text = element_text(size = bs - 1.3), legend.title = element_text(size = bs - 1),
        legend.key.size = unit(2.6, "mm"), legend.background = element_blank(), plot.margin = margin(3, 4, 2, 2))
MET <- "#B0403F"
d <- fread(file.path(DAT, "fig4_pm_categories.tsv"))
LV <- c("Metastasis-specific", "Primary-specific", "Metastasis-enriched", "Primary-enriched", "M–P switch", "Shared")
cl <- c(Mspec = LV[1], Pspec = LV[2], Menr = LV[3], Penr = LV[4], switch = LV[5], none = LV[6])
d[, cls := factor(cl[cat], LV)]
d[, lab := sub(" \\(.*\\)$", "", gene)]; d[modality == "focal_loh", lab := paste(sub("^.*\\((.*)\\)$", "\\1", gene), "LOH")]; d[modality %in% c("broad_gain", "broad_loh"), lab := paste(lab, ifelse(modality == "broad_gain", "gain", "LOH"))]
d[modality == "focal_amp", lab := paste(lab, "amp")]; d[modality == "focal_del", lab := paste(lab, "del")]
Z <- 42                                                   # zoom window |ES| <= Z (inset)
d[, spec := cat %in% c("Mspec", "Pspec", "switch")]
cand <- d[cls != "Shared"][order(-abs(dist_obs), -spec)]  # ranked by |distance|; specific only breaks ties
cand[, inside := abs(es_P) <= Z & abs(es_M) <= Z]
FORCE <- c("CDKN2A del", "9q LOH", "4p LOH")              # named in the text: always labelled, wherever they fall
out <- unique(rbind(cand[inside == FALSE][1:6], cand[inside == FALSE & lab %in% FORCE]), by = "lab")   # six outside the window (+ forced)
inn <- unique(rbind(cand[inside == TRUE][1:4], cand[inside == TRUE & lab %in% FORCE]), by = "lab")     # four inside it (+ forced)
# hand-chosen label centres (data coordinates) in empty space; outer labels sit outside the window and off the inset,
# inset labels in the inset's empty corners; loci without a position fall back to a nudge away from the diagonal
OPOS <- data.table(lab = c("4p LOH", "7q gain", "9p LOH", "8p LOH", "CDKN2A del", "H1-2 amp", "18p gain"),
                   lx = c(22, -98, -105, -58, 78, -55, 43), ly = c(-122, -80, 18, 63, 118, -112, 138))
IPOS <- data.table(lab = c("9q LOH", "14q LOH", "TBX3 LOH", "HNF1A LOH"), lx = c(-28, -24, 22, 22), ly = c(14, 37.5, -29, -38.5))
place <- function(x, P, dx, dy) { x <- merge(x, P, by = "lab", all.x = TRUE, sort = FALSE)
  if (anyNA(x$lx)) warning("no hand position for ", paste(x[is.na(lx), lab], collapse = ", "), " - using the default nudge")
  x[is.na(lx), `:=`(lx = es_P + ifelse(es_M > es_P, -dx, dx), ly = es_M + ifelse(es_M > es_P, dy, -dy))]; x }
out <- place(out, OPOS, 60, 30); inn <- place(inn, IPOS, 13, 9)
out[, reg := ifelse(ly > Z, "top", ifelse(lx < -Z, "left", "bottom"))]   # region-bounded repel limits, all excluding the window
lab_repel <- function(dd, sz, xl, yl, bp = 0.15) geom_text_repel(data = dd, aes(label = lab), nudge_x = dd$lx - dd$es_P, nudge_y = dd$ly - dd$es_M,
  family = FONT, size = sz, fontface = "italic", colour = "grey15", min.segment.length = 0, segment.size = 0.2, segment.colour = "grey35",
  box.padding = bp, point.padding = 0.15, force = 0, force_pull = 0, seed = 4, max.overlaps = Inf, xlim = xl, ylim = yl, bg.color = "white", bg.r = 0.15)
LEG <- "Stronger WGD association\n(CO or ME) in:"
PB <- "#3F6A9E"
cols <- setNames(c(MET, PB, MET, PB, "#E69F00", "grey75"), LV); shp <- setNames(c(16, 16, 1, 1, 16, 16), LV)
lim <- c(min(c(d$es_P, d$es_M)) - 25, max(c(d$es_P, d$es_M)) + 10); CO <- "#3C8D5A"; ME <- "#7B5EA7"
pts <- function(g, sz = 1.05, ssz = 0.55) g + geom_point(data = d[cls == "Shared"], colour = "grey78", size = ssz) +
  geom_point(data = d[cls != "Shared"], aes(colour = cls, shape = cls), size = sz, stroke = 0.45) +
  scale_colour_manual(values = cols, name = LEG, drop = TRUE) + scale_shape_manual(values = shp, name = LEG, drop = TRUE)
g2 <- ggplot(d, aes(es_P, es_M)) + geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.3, linetype = "22") +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.25) +
  annotate("rect", xmin = -Z, xmax = Z, ymin = -Z, ymax = Z, fill = "#F3D34A", alpha = 0.55, colour = NA)
g2 <- pts(g2) +
  lab_repel(out[reg == "top"], 1.75, c(lim[1] + 4, lim[2] - 2), c(Z + 8, lim[2] - 2)) +
  lab_repel(out[reg == "left"], 1.75, c(lim[1] + 4, -Z - 6), c(lim[1] + 6, lim[2] - 2)) +
  lab_repel(out[reg == "bottom"], 1.75, c(lim[1] + 4, 48), c(lim[1] + 6, -Z - 6)) +
  annotate("segment", x = lim[1], xend = lim[2], y = lim[1], yend = lim[1], arrow = arrow(length = unit(1, "mm"), ends = "both", type = "closed"), colour = "grey45", linewidth = 0.25) +
  annotate("segment", x = lim[1], xend = lim[1], y = lim[1], yend = lim[2], arrow = arrow(length = unit(1, "mm"), ends = "both", type = "closed"), colour = "grey45", linewidth = 0.25) +
  annotate("text", x = c(lim[1] + 16, lim[2] - 4), y = lim[1] + 14, label = c("ME", "CO"), colour = c(ME, CO), hjust = c(0, 1), vjust = 0, family = FONT, size = 2) +
  annotate("text", x = lim[1] + 16, y = lim[2] - 4, label = "CO", colour = CO, hjust = 0, vjust = 1, family = FONT, size = 2) +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) + labs(x = "Primary, mean ES", y = "Metastasis, mean ES") + theme_f4() +
  theme(legend.position = c(0.27, 0.78), legend.key.size = unit(1.8, "mm"), legend.text = element_text(size = 5.0),
        legend.title = element_text(size = 5.0, lineheight = 0.9), legend.spacing.y = unit(0, "mm"))
gz <- ggplot(d[abs(es_P) <= Z & abs(es_M) <= Z], aes(es_P, es_M)) + geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.3, linetype = "22")
gz <- pts(gz, 0.8, 0.4) + lab_repel(inn, 1.6, c(-Z + 1, Z - 1), c(-Z + 1, Z - 1), bp = 0.1) +
  coord_equal(xlim = c(-Z, Z), ylim = c(-Z, Z), expand = FALSE) + theme_void(base_family = FONT) +
  theme(legend.position = "none", panel.background = element_rect(fill = "#FBF1C4", colour = "grey40", linewidth = 0.3))
g2 <- g2 + inset_element(gz, left = 0.56, bottom = 0.11, right = 0.99, top = 0.52, align_to = "panel")
ggsave(file.path(OUT, "fig4f.png"), g2, width = 68, height = 66, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig4f.pdf"), g2, width = 68, height = 66, units = "mm", device = cairo_pdf, bg = "white")
