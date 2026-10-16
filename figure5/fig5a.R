# Fig. 5a: the 50 MSigDB Hallmark programmes ranked by their WGD+ vs WGD− normalised enrichment score (NES) in the TCGA
# transcriptome and the CPTAC proteome and phosphoproteome (cancer type as the only covariate; the layer's sample size is
# printed under its name). Red = immune / IFN programmes (the five with the best mean rank across the three layers are named),
# dark = proliferation programmes, light grey = other; filled = FDR < 0.05, open = n.s. The key is drawn as its own strip
# under the panel.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
DAT <- "figure5/data"; OUT <- "results/figure5"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
ACC <- "#B5533C"; DARK <- "grey30"; LIGHT <- "grey75"
TXT <- 7; TICK <- 6; LAB <- 2.15                          # 2.15 mm = 6.1 pt
H <- read.delim(file.path(DAT, "fig5a_hallmark_enrichment.tsv")); S <- read.delim(file.path(DAT, "fig5a_layer_sizes.tsv"))
IMM <- c("Interferon Gamma Response", "Interferon Alpha Response", "Allograft Rejection", "Inflammatory Response",
         "IL-6/JAK/STAT3 Signaling", "Complement", "TNF-alpha Signaling via NF-kB", "IL-2/STAT5 Signaling")
PROL <- c("E2F Targets", "G2-M Checkpoint", "Myc Targets V1", "Myc Targets V2", "Mitotic Spindle")

LV <- c("RNA", "protein", "phospho"); LB <- c("Transcriptome", "Proteome", "Phosphoproteome")
A <- H[H$layer %in% LV, ]
A$layer <- factor(A$layer, levels = LV, labels = LB)
CONS <- c("TCGA", "CPTAC", "CPTAC")          # consortium per layer, shown under the layer name
A$class <- ifelse(A$pathway %in% IMM, "Immune", ifelse(A$pathway %in% PROL, "Proliferation", "Other"))
A$rk <- A$rank_down; A$sig <- A$padj < 0.05        # y = the NES rank within the layer (1 = most down); filled = FDR < 0.05
# the five labelled immune programmes: the five of the eight with the best mean rank across the three layers
A$lab <- c("Allograft Rejection" = "Allograft rejection",
           "Interferon Gamma Response" = "IFN-γ response",
           "Inflammatory Response" = "Inflammatory\nresponse",
           "IL-6/JAK/STAT3 Signaling" = "IL-6 / JAK / STAT3",
           "Complement" = "Complement")[A$pathway]
A$lab[is.na(A$lab)] <- ""
LL <- do.call(rbind, lapply(split(A[A$lab != "", ], A$layer[A$lab != ""]), function(x) {
  x <- x[order(x$rk), ]; x$ly <- min(x$rk) + (seq_len(nrow(x)) - 1) * 7; x$xl <- 1.0; x }))
nl <- setNames(sprintf("%s\n%s %s", LB, CONS, formatC(S$n[match(LV, S$layer)], big.mark = ",", format = "d")), LB)

strips <- ggplot(A, aes(NES, rk)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55", linetype = "22") +
  geom_point(data = A[A$class == "Other", ], aes(shape = sig), colour = LIGHT, size = 1.05, stroke = 0.3) +
  geom_point(data = A[A$class == "Proliferation", ], aes(shape = sig), colour = DARK, size = 1.2, stroke = 0.4) +
  geom_point(data = A[A$class == "Immune", ], aes(shape = sig), colour = ACC, size = 1.5, stroke = 0.45) +
  geom_segment(data = LL, aes(x = NES, y = rk, xend = xl - 0.15, yend = ly), linewidth = 0.2, colour = "grey60") +
  geom_text(data = LL, aes(x = xl, y = ly, label = lab), family = FONT, size = LAB, colour = ACC, hjust = 0, lineheight = 0.85) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
  facet_wrap(~ layer, nrow = 1, labeller = as_labeller(nl)) +
  scale_y_continuous(breaks = c(1, 50), labels = c("down", "up"), expand = expansion(add = 2)) +
  scale_x_continuous(breaks = c(-2, 0, 2), expand = expansion(mult = c(0.08, 1.25))) +
  coord_cartesian(clip = "off") +
  labs(x = "NES, WGD+ vs WGD−", y = "Rank of 50 Hallmark programmes") +
  theme_classic(base_family = FONT, base_size = TXT) +
  theme(strip.background = element_blank(), strip.text = element_text(size = TXT, colour = "grey10", lineheight = 0.9, margin = margin(0, 0, 1, 0)),
        axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        axis.text = element_text(size = TICK, colour = "grey10"), axis.title = element_text(size = TXT),
        panel.spacing.x = unit(1.5, "mm"), plot.margin = margin(2, 2, 0, 2))

## one-line key, drawn as its own strip (glyph + text, x positions in mm so the spacing is exact)
key_items <- data.frame(txt = c("FDR < 0.05", "n.s.", "immune / IFN", "proliferation", "other"),
                        shape = c(16, 1, 16, 16, 16), col = c("grey45", "grey45", ACC, DARK, LIGHT))
wmm <- function(s) systemfonts::string_width(s, family = FONT, size = TICK * 10, res = 72) / 10 / 72 * 25.4
key_items$w <- vapply(key_items$txt, wmm, 0)
GAP_G <- 1.3; GAP_I <- 3.0                                  # glyph->text, item->item (mm)
x <- 0; gx <- tx <- numeric(nrow(key_items))
for (i in seq_len(nrow(key_items))) {
  if (i == 3) x <- x + 1.5                                  # small break between the FDR pair and the colour triple
  gx[i] <- x + 0.6; tx[i] <- x + 0.6 + GAP_G; x <- tx[i] + key_items$w[i] + GAP_I }
PW <- 96 - 2 * 2 / 72.27 * 25.4                              # key panel width in mm (96 mm minus 2 pt side margins)
off <- (PW - (x - GAP_I)) / 2                               # centre the key under the panel
key_items$gx <- gx + off; key_items$tx <- tx + off
key <- ggplot(key_items) +
  geom_point(aes(gx, 0, shape = I(shape), colour = I(col)), size = 1.2, stroke = 0.4) +
  geom_text(aes(tx, 0, label = txt), hjust = 0, family = FONT, size = TICK / .pt, colour = "grey15") +
  scale_x_continuous(limits = c(0, PW), expand = expansion(0)) + scale_y_continuous(limits = c(-1, 1)) +
  theme_void(base_family = FONT) + theme(plot.margin = margin(0, 2, 1, 2))

p <- strips / wrap_elements(full = key) + plot_layout(heights = c(1, 0.05))   # full width: not aligned to the facet panels

W <- 96; Hh <- 62
ggsave(file.path(OUT, "fig5a.png"), p, width = W, height = Hh, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig5a.pdf"), p, width = W, height = Hh, units = "mm", device = cairo_pdf, bg = "white")
