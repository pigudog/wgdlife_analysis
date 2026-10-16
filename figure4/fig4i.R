# Fig. 4i: where the genome doubled in matched primary–metastasis pairs of MSK-IMPACT with WGD in at least one sample:
# WGD in both samples, in the metastasis only or in the primary tumour only, per cancer type (rows with more than ten such
# pairs) and overall; right columns = metastasis-only : primary-only and the number of pairs.
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
DAT <- "figure4/data"; OUT <- "results/figure4"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
PRI <- "#8FA9C4"; MET <- "#B0403F"
theme_f4 <- function(bs = 7) theme_classic(base_family = FONT, base_size = bs) +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35), axis.text = element_text(size = bs - 0.8, colour = "grey10"),
        axis.title = element_text(size = bs - 0.3), legend.text = element_text(size = bs - 1.3), legend.title = element_text(size = bs - 1),
        legend.key.size = unit(2.6, "mm"), legend.background = element_blank(), plot.margin = margin(3, 4, 2, 2))
SHORT <- c("Non-Small Cell Lung Cancer" = "NSCLC", "Breast Cancer" = "Breast", "Colorectal Cancer" = "Colorectal", "Prostate Cancer" = "Prostate",
           "Pancreatic Cancer" = "Pancreatic", "Bladder Cancer" = "Bladder", "Endometrial Cancer" = "Endometrial", "Ovarian Cancer" = "Ovarian",
           "Esophagogastric Cancer" = "Oesophagogastric", "Hepatobiliary Cancer" = "Hepatobiliary", "Renal Cell Carcinoma" = "Renal", "Head and Neck Cancer" = "Head & neck")
LAB <- c("WGD in both", "Metastasis only", "Primary only"); COL <- c("#454545", MET, PRI)
TB <- fread(file.path(DAT, "fig4i_pair_counts.tsv"))[, .(cancer, n, n1 = WGD_in_both, n2 = metastasis_only, n3 = primary_only)]
CT <- TB[cancer != "All"]; CT[, label := SHORT[cancer]]
CT[, s2 := n2 / n]; setorder(CT, s2); CT[, y := as.numeric(seq_len(.N))]
AL <- TB[cancer == "All"]; AL[, `:=`(label = "All", y = nrow(CT) + 1.45)]
TB <- rbind(CT, AL, fill = TRUE)
seg <- rbindlist(lapply(1:3, function(k) TB[, .(label, y, k, frac = get(paste0("n", k)) / n)]))
seg[, kk := k]; setorder(seg, y, kk)                                          # left -> right: both | met only | pri only
seg[, xmax := 100 * cumsum(frac), by = y][, xmin := xmax - 100 * frac]
pct <- seg[label == "All" & frac >= 0.12, .(x = (xmin + xmax) / 2, y, txt = sprintf("%.0f%%", 100 * frac), col = ifelse(kk == 3, "grey10", "white"))]
TS <- 6; TK <- 5.5; BH <- 0.36
XR <- 105; XN <- 150                                                          # right-hand text columns: ratio | n
p <- ggplot() +
  geom_rect(data = seg, aes(xmin = xmin, xmax = xmax, ymin = y - BH, ymax = y + BH, fill = factor(kk)), colour = "white", linewidth = 0.25) +
  geom_text(data = pct, aes(x, y, label = txt, colour = col), family = FONT, size = TK / .pt) +
  geom_text(data = TB, aes(-2, y, label = label, fontface = ifelse(label == "All", "bold", "plain")), hjust = 1, family = FONT, size = TS / .pt, colour = "grey10") +
  geom_text(data = TB, aes(XR, y, label = sprintf("%d : %d", n2, n3)), hjust = 0, family = FONT, size = TK / .pt, colour = "grey15") +
  geom_text(data = TB, aes(XN, y, label = n), hjust = 1, family = FONT, size = TK / .pt, colour = "grey45") +
  annotate("text", x = XR, y = max(TB$y) + 0.75, label = "Met : Pri", hjust = 0, family = FONT, size = TK / .pt, colour = "grey45") +
  annotate("text", x = XN, y = max(TB$y) + 0.75, label = "n", hjust = 1, family = FONT, size = TK / .pt, colour = "grey45") +
  annotate("segment", x = 0, xend = 100, y = 0.4, yend = 0.4, linewidth = 0.3, colour = "black") +
  scale_fill_manual(values = setNames(COL, 1:3), breaks = 1:3, labels = LAB, name = NULL) +
  scale_colour_identity() +
  scale_x_continuous(breaks = c(0, 50, 100), limits = c(-62, XN + 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, max(TB$y) + 1.1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Pairs with WGD in at least one sample (%)", y = NULL) +
  theme_f4(TS) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_line(linewidth = 0.3), axis.text.x = element_text(size = TK),
        axis.title.x = element_text(size = TS, hjust = 0.58),
        legend.position = "top", legend.justification = c(0.6, 0), legend.text = element_text(size = TK), legend.key.size = unit(2.2, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -4, 0), plot.margin = margin(2, 2, 1, 1))
ggsave(file.path(OUT, "fig4i.png"), p, width = 70, height = 52, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig4i.pdf"), p, width = 70, height = 52, units = "mm", device = cairo_pdf, bg = "white")
