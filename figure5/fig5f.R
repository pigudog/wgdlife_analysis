# Fig. 5f: hazard ratios for death after immune checkpoint blockade (Cox models, seven carcinomas): the WGD estimate alone
# and adjusted for TMB, for clinical covariates, for FGA and for MSI-H, then the two competing predictors on their own scale
# (TMB per doubling, FGA per 10 percentage points). Log-scale axis with the dashed null, a "Better | Worse" header and a
# right-hand HR (95% CI) column.
suppressPackageStartupMessages({ library(ggplot2) })
DAT <- "figure5/data"; OUT <- "results/figure5"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
BASE <- 7; TICK <- 6; LW <- 0.45
W <- 53; H <- 60
RED <- "#C0504D"; DARK <- "grey25"

cx <- read.delim(file.path(DAT, "fig5f_cox_models.tsv"))
pick <- function(model, term) { r <- cx[cx$model == model & cx$term == term, ][1, ]; stopifnot(!is.na(r$hr)); r }
row <- function(lab, model, term, col, grp) data.frame(lab = lab, grp = grp, col = col,
                                                       pick(model, term)[, c("hr", "lo", "hi")])
d <- rbind(
  row("WGD alone",        "M0_wgd_alone",        "wgd",    RED,  "WGD"),
  row("WGD + TMB",        "M1_wgd_tmb",          "wgd",    RED,  "WGD"),
  row("WGD + clinical",   "M2_wgd_tmb_clin",     "wgd",    RED,  "WGD"),
  row("WGD + FGA",        "M3_wgd_tmb_fga_clin", "wgd",    RED,  "WGD"),
  row("WGD + MSI-H",      "M2_wgd_tmb_clin_msi", "wgd",    RED,  "WGD"),
  row("TMB per doubling", "M2_wgd_tmb_clin",     "logTMB", DARK, "TMB"),
  row("FGA per 10 pp",    "M3_wgd_tmb_fga_clin", "fga10",  DARK, "FGA"))
d$y  <- c(7.6, 6.6, 5.6, 4.6, 3.6, 2.1, 1.1)                        # 5 WGD rows, gap, the two predictors
d$txt <- sprintf("%.2f (%.2f–%.2f)", d$hr, d$lo, d$hi)

XL <- c(0.70, 3.15); XT <- 1.66   # axis drawn only to XA; the HR column lives inside the panel
XA <- 1.56                                      # HR axis range; x of the HR text column
p <- ggplot(d, aes(hr, y, colour = col)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = LW, colour = "grey50") +
  annotate("segment", x = XL[1], xend = XA, y = 0.45, yend = 0.45, linewidth = LW, colour = "black") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.5, orientation = "y") +
  geom_point(size = 1.3) +
  geom_text(aes(x = XT, label = txt), hjust = 0, family = FONT, size = 1.9, show.legend = FALSE) +
  annotate("text", x = 0.985, y = 8.5, label = "Better", hjust = 1, family = FONT, size = 1.85, colour = "grey40") +
  annotate("text", x = 1.015, y = 8.5, label = "Worse", hjust = 0, family = FONT, size = 1.85, colour = "grey40") +
  annotate("text", x = XT, y = 8.5, label = "HR (95% CI)", hjust = 0, family = FONT, size = 1.85, colour = "grey40") +
  scale_colour_identity() +
  scale_x_continuous(transform = "log", breaks = c(0.75, 1, 1.25, 1.5), labels = c("0.75", "1", "1.25", "1.5")) +
  scale_y_continuous(breaks = d$y, labels = d$lab, expand = c(0, 0)) +
  coord_cartesian(xlim = XL, ylim = c(0.45, 8.9), clip = "off") +
  labs(x = "Hazard ratio for death (log scale)", y = NULL) +
  theme_classic(base_family = FONT, base_size = BASE) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_line(linewidth = LW),
        axis.ticks.length = unit(0.9, "mm"),
        axis.text.y = element_text(size = TICK, colour = "black", hjust = 1, margin = margin(r = 1)),
        axis.text.x = element_text(size = TICK, colour = "black"),
        axis.title.x = element_text(size = TICK, margin = margin(t = 1), hjust = 0.30),
        plot.margin = margin(2, 1.5, 1, 1))
ggsave(file.path(OUT, "fig5f.png"), p, width = W, height = H, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig5f.pdf"), p, width = W, height = H, units = "mm", device = cairo_pdf, bg = "white")
