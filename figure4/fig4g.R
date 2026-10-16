# Fig. 4g: genomic events that add to WGD in predicting metastasis in MSK-MET primary tumours. Rows: events retained by
# stability selection in the main-effect model and associated with WGD in Fig. 4f. Columns: mean WGD-association effect size
# of the event in primary-site and metastasis-site samples (green CO / purple ME), the primary-versus-metastasis category of
# Fig. 4f, and the odds ratio for metastasis beyond WGD (purity-adjusted logistic model; dark orange > 1, teal < 1).
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
DAT <- "figure4/data"; OUT <- "results/figure4"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
MET <- "#B0403F"
W <- fread(file.path(DAT, "fig4g_rows.tsv"))
W[, `:=`(OR = OR_beyond_WGD_purity_adj, cat = cat_4f)]
CO <- "#3C8D5A"; ME <- "#7B5EA7"; PB <- "#3F6A9E"                            # CO / ME and the primary blue of Fig. 4f; MET = red
POSC <- "#C8641E"; NEGC <- "#2B8A8A"                                        # OR pair (dark orange > 1 / teal < 1)
NAFILL <- "grey90"; SHARED <- "grey75"; TS <- 6; TK <- 5.5                   # text sizes (pt): body / key ticks + cell numbers
stopifnot(!anyNA(W$OR), all(W$OR_lo > 1 | W$OR_hi < 1))
W[, blk := factor(ifelse(OR > 1, "up", "down"), c("up", "down"))]; setorder(W, blk, -OR)
n <- nrow(W); stopifnot(n == 11); GAP <- 0.45; W[, y := as.numeric(rev(seq_len(n)))]; W[blk == "up", y := y + GAP]        # top row = first
# colour maps (identity fill; ES and log OR each on their own symmetric scale)
L_es <- 100   # ES colour clipped at +-100
OR_LIM <- 2.5; L_or <- log(OR_LIM); stopifnot(max(abs(log(W$OR))) <= L_or)  # symmetric log-OR scale 1/2.5 … 2.5
pal <- function(lo, hi) scales::gradient_n_pal(c(lo, "#F4F4F4", hi), values = c(-1, 0, 1))
div <- function(v, L, lo, hi) ifelse(is.na(v), NAFILL, pal(lo, hi)(pmax(-1, pmin(1, v / L))))
W[, c_P := div(es_P, L_es, ME, CO)]; W[, c_M := div(es_M, L_es, ME, CO)]; W[, c_O := div(log(OR), L_or, NEGC, POSC)]
W[, txt_col := ifelse(abs(log(OR)) / L_or > 0.55, "white", "grey10")]     # printed OR: white on the saturated cells
# geometry (data units): x columns; y rows 1..n (+GAP), headers above, keys below 0.5
X <- c(P = 1, M = 2, cat = 3.15, o = 4.45); TW <- 0.92; TH <- 0.9; TWo <- 1.25
tiles <- rbind(W[, .(x = X["P"], y, fill = c_P, w = TW)], W[, .(x = X["M"], y, fill = c_M, w = TW)], W[, .(x = X["o"], y, fill = c_O, w = TWo)],
               NULL)
cats <- copy(W); cats[, cat_col := fcase(cat %in% c("Mspec", "Menr"), MET, cat %in% c("Pspec", "Penr"), PB, default = SHARED)]
cats[, cat_fill := fcase(cat %in% c("Mspec", "Pspec"), cat_col, cat %in% c("Menr", "Penr"), "white", default = SHARED)]
stopifnot(all(cats$cat %in% c("Mspec", "Menr", "Pspec", "Penr", "none")))
hdr <- data.table(x = X, lab = c("P", "M", "Δ", "OR"))
ytop <- n + GAP + 0.5; HU <- 0.75                                             # one text line
span <- data.table(x0 = c(X["P"] - 0.45, X["o"] - TWo / 2), x1 = c(X["cat"] + 0.3, X["o"] + TWo / 2),
                   lab = c("Mean ES\nwith WGD", "Metastasis\nbeyond WGD"))
ysp <- ytop + 0.2 + HU + 0.25                                                  # span underline
span[, lx := c((X["P"] + X["M"]) / 2, X["o"])]                              # label centres
keybar <- function(x0, x1, yc, lo, hi, k = 40, hh = 0.24) { s <- seq(-1, 1, length.out = k); dx <- (x1 - x0) / k
  data.table(xmin = x0 + (seq_len(k) - 1) * dx, xmax = x0 + seq_len(k) * dx, ymin = yc - hh, ymax = yc + hh, fill = pal(lo, hi)(s)) }
# key block: left column = two colour bars (Mean ES / OR beyond WGD), right column = category items + the grey key
yA <- -1.35; yB <- -3.75; KB <- 2.6; kx <- -4.85; bx <- kx + 0.78               # bar rows; bar length; key column start; bar start
kb <- rbind(keybar(bx, bx + KB, yA, ME, CO), keybar(bx, bx + KB, yB, NEGC, POSC))
ktitle <- data.table(x = kx, y = c(yA, yB) + 0.78, lab = c("Mean ES", "OR beyond WGD"))
ktick <- data.table(x = rep(c(bx, bx + KB / 2, bx + KB), 2), y = rep(c(yA, yB) - 0.7, each = 3),
                    lab = c(sprintf("≤−%d", L_es), "0", sprintf("≥%d", L_es), sprintf("%.1f", 1 / OR_LIM), "1", sprintf("%.1f", OR_LIM)))
kco <- data.table(x = c(kx, bx + KB + 0.1), y = yA, lab = c("ME", "CO"), col = c(ME, CO), hj = c(0, 0))
ckey <- data.table(x = 0.85, y = yA + 0.78 - 0.95 * (1:5), col = c(MET, MET, PB, PB, SHARED), fill = c(MET, "white", PB, "white", SHARED),
                   lab = c("Metastasis-specific", "Metastasis-enriched", "Primary-specific", "Primary-enriched", "Shared"))
ctitle <- data.table(x = 0.85 - 0.21, y = yA + 0.78, lab = "Δ (P vs M)")
ybot <- min(c(ckey$y, yB - 0.7)) - 0.8
SQ <- 0.28                                                                     # half side of a ~1.8 mm square at ~3.2 mm per row
g <- ggplot() +
  geom_rect(data = tiles, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - TH / 2, ymax = y + TH / 2, fill = fill), colour = "grey82", linewidth = 0.12) +
  geom_rect(data = cats, aes(xmin = X["cat"] - 0.21, xmax = X["cat"] + 0.21, ymin = y - SQ, ymax = y + SQ, fill = cat_fill, colour = cat_col), linewidth = 0.45) +
  geom_text(data = W, aes(X["o"], y, label = sprintf("%.2f", OR), colour = txt_col), family = FONT, size = TK / .pt) +
  geom_text(data = W, aes(0.4, y, label = lab), hjust = 1, family = FONT, size = TS / .pt, colour = "grey10") +
  geom_text(data = hdr, aes(x, ytop + 0.2, label = lab), vjust = 0, family = FONT, size = TS / .pt, colour = "grey10") +
  geom_segment(data = span, aes(x = x0, xend = x1, y = ysp, yend = ysp), colour = "grey40", linewidth = 0.3) +
  geom_text(data = span, aes(lx, ysp + 0.35, label = lab), vjust = 0, family = FONT, size = TS / .pt, lineheight = 0.85, colour = "grey10") +
  geom_rect(data = kb, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), colour = NA) +
  annotate("rect", xmin = bx, xmax = bx + KB, ymin = c(yA, yB) - 0.24, ymax = c(yA, yB) + 0.24, fill = NA, colour = "grey60", linewidth = 0.2) +
  geom_text(data = ktitle, aes(x, y, label = lab), hjust = 0, family = FONT, size = TS / .pt, colour = "grey10") +
  geom_text(data = ktick, aes(x, y, label = lab), family = FONT, size = TK / .pt, colour = "grey25") +
  geom_text(data = kco, aes(x, y, label = lab, colour = col, hjust = hj), family = FONT, size = TS / .pt, fontface = "bold") +
  geom_rect(data = ckey, aes(xmin = x - 0.21, xmax = x + 0.21, ymin = y - SQ, ymax = y + SQ, fill = fill, colour = col), linewidth = 0.45) +
  geom_text(data = ckey, aes(x + 0.45, y, label = lab), hjust = 0, family = FONT, size = TS / .pt, colour = "grey10") +
  geom_text(data = ctitle, aes(x, y, label = lab), hjust = 0, family = FONT, size = TS / .pt, colour = "grey10") +
  scale_fill_identity() + scale_colour_identity() +
  coord_cartesian(xlim = c(-4.95, 6.3), ylim = c(ybot, ysp + 2.2), expand = FALSE, clip = "off") +
  theme_void(base_family = FONT) + theme(plot.margin = margin(1, 1, 1, 1), plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(OUT, "fig4g.png"), g, width = 49, height = 70, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig4g.pdf"), g, width = 49, height = 70, units = "mm", device = cairo_pdf, bg = "white")
