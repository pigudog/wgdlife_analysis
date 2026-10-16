# Fig. 4a: WGD prevalence along the path to metastasis in the 12 carcinoma types. Left, MSK-IMPACT: primary tumours of
# patients with localized disease, primary tumours of patients with metastatic disease, and metastases (* = trend across the
# three groups, BH-adjusted; All = logistic trend adjusted for cancer type). Right, whole genomes: PCAWG primaries versus
# Hartwig metastases (* = two-sided Fisher test, BH-adjusted; All = Cochran–Mantel–Haenszel).
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(patchwork) })
DAT <- "figure4/data"; OUT <- "results/figure4"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
PRI <- "#8FA9C4"; MET <- "#B0403F"
theme_f4 <- function(bs = 7) theme_classic(base_family = FONT, base_size = bs) +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35), axis.text = element_text(size = bs - 0.8, colour = "grey10"),
        axis.title = element_text(size = bs - 0.3), legend.text = element_text(size = bs - 1.3), legend.title = element_text(size = bs - 1),
        legend.key.size = unit(2.6, "mm"), legend.background = element_blank(), plot.margin = margin(3, 4, 2, 2))
dM <- fread(file.path(DAT, "fig4a_msk_by_disease_state.tsv")); dW <- fread(file.path(DAT, "fig4a_wgs_primary_vs_metastasis.tsv"))
GL <- c("Primary, localized disease", "Primary, metastatic disease", "Metastasis")
lv <- c("All", dM[lab != "All"][order(p3)]$lab); dM[, lab := factor(lab, lv)]; dW[, lab := factor(lab, lv)]
LM <- melt(dM, id.vars = "lab", measure.vars = c("p1", "p2", "p3"), variable.name = "grp", value.name = "pct"); LM[, grp := factor(GL[as.integer(sub("p", "", grp))], GL)]
sty <- theme(legend.position = "bottom", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.5, "mm"),
             legend.key.height = unit(2.2, "mm"), strip.background = element_blank(), strip.text = element_text(size = 6.2))
g1 <- ggplot(dM, aes(y = lab)) + geom_segment(aes(x = pmin(p1, p2, p3), xend = pmax(p1, p2, p3), yend = lab), colour = "grey82", linewidth = 0.9, lineend = "round") +
  geom_point(data = LM, aes(x = pct, colour = grp), size = 1.45) +
  geom_text(aes(x = pmax(p1, p2, p3) + 3.5, label = star), family = FONT, size = 2.3, vjust = 0.75, colour = "grey15") +
  geom_hline(yintercept = 1.5, colour = "grey88", linewidth = 0.3) + facet_wrap(~ "MSK-IMPACT") +
  scale_colour_manual(values = setNames(c(PRI, "#9A7AA0", MET), GL), name = NULL) + scale_x_continuous(limits = c(0, 70), breaks = seq(0, 60, 20), expand = c(0, 0)) +
  labs(x = "WGD+ samples (%)", y = NULL) + theme_f4() + sty + theme(axis.text.y = element_text(size = 5.8, face = ifelse(lv == "All", "bold", "plain")))
g2 <- ggplot(dW, aes(y = lab)) + geom_segment(aes(x = P, xend = M, yend = lab), colour = "grey82", linewidth = 0.9, lineend = "round") +
  geom_point(aes(x = P, colour = "PCAWG primary"), size = 1.45) + geom_point(aes(x = M, colour = "Hartwig metastasis"), size = 1.45) +
  geom_text(aes(x = pmax(P, M) + 3.5, label = star), family = FONT, size = 2.3, vjust = 0.75, colour = "grey15") +
  geom_hline(yintercept = 1.5, colour = "grey88", linewidth = 0.3) + facet_wrap(~ "PCAWG / Hartwig") +
  scale_colour_manual(values = c(`PCAWG primary` = PRI, `Hartwig metastasis` = MET), breaks = c("PCAWG primary", "Hartwig metastasis"), name = NULL) +
  scale_x_continuous(limits = c(0, 85), breaks = seq(0, 75, 25), expand = c(0, 0)) + labs(x = "WGD+ samples (%)", y = NULL) + theme_f4() + sty +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
g <- g1 + g2 + plot_layout(widths = c(1.15, 1))
ggsave(file.path(OUT, "fig4a.png"), g, width = 98, height = 64, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig4a.pdf"), g, width = 98, height = 64, units = "mm", device = cairo_pdf, bg = "white")
