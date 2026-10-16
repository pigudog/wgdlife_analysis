# Fig. 5b: composition of the immune compartment in WGD+ versus WGD− tumours (TCGA, 12 carcinomas): 162 published immune
# features, x = WGD effect in the ploidy-free (compositional) currency, proliferation-adjusted; y = −log10 P on a square-root
# axis. Point fill = significance tier of the features below the empirical FDR 0.05 (1,000 permutations); dashed lines = the
# P at which that FDR is reached and the smallest |effect| that clears it. Labelled features are coloured by cell-type family.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
DAT <- "figure5/data"; OUT <- "results/figure5"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
BASE <- 7; TICK <- 6; LW <- 0.5 / 2.13
SYMFAM <- cairoSymbolFont(FONT, usePUA = FALSE)

th <- theme_classic(base_size = BASE, base_family = FONT) +
  theme(axis.line = element_line(linewidth = LW), axis.ticks = element_line(linewidth = LW),
        axis.ticks.length = unit(1, "mm"), axis.text = element_text(size = TICK, colour = "black"),
        axis.title = element_text(size = BASE), legend.text = element_text(size = TICK),
        legend.title = element_blank(), legend.key.size = unit(2.2, "mm"),
        legend.key.spacing.x = unit(1, "mm"), legend.key.spacing.y = unit(0, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.8, "mm"),
        plot.title = element_text(size = BASE, hjust = 0, margin = margin(0, 0, 1, 0)),
        plot.subtitle = element_text(size = TICK, hjust = 0, colour = "grey30",
                                     margin = margin(0, 0, 1.5, 0)),
        plot.margin = margin(1.5, 1.5, 1, 1.5))

# muted, one hue per family; red is the accent, kept for IFN / cytotoxic / APM
FAMCOL <- c("T / CD8" = "#4C72A0", "NK" = "#5F9EA0", "B / plasma" = "#8C77A8",
            "Myeloid / macrophage" = "#C4803E", "Dendritic" = "#C08BA8", "Stromal" = "#7F9A5E",
            "IFN / cytotoxic / APM" = "#B0403F", "Other" = "#9E9E9E")
SHORT <- c("T / CD8" = "T / CD8", "NK" = "NK", "B / plasma" = "B / plasma",
           "Myeloid / macrophage" = "myeloid", "Dendritic" = "dendritic", "Stromal" = "stromal",
           "IFN / cytotoxic / APM" = "IFN / APM", "Other" = "other")

t <- fread(file.path(DAT, "fig5b_immune_features.tsv"))
st <- read.delim(file.path(DAT, "fig5b_panel_stats.tsv"))   # the dashed thresholds (largest P and smallest |effect| among the features below the permutation FDR 0.05) and the subtitle counts
stopifnot(nrow(t) == st$n_features)
t[, fam := factor(SHORT[family], levels = SHORT)]
t[, sig := fdr_perm_compositional < 0.05]
PCUT <- st$p_cut_fdr05                  # the dashed line: P at which the empirical (1,000-permutation) FDR reaches 0.05
t[, y := -log10(p_compositional)]

# 12 labelled features: most extreme AND interpretable, spread across the source methods
LAB <- c("Wolf|T_cell_PCA_16704732" = "T cell",
         "Wolf|LCK_19272155" = "LCK",
         "Attractors|G_CD3E" = "CD3E",
         "Wolf|STAT1_score" = "STAT1",
         "Cibersort_rel|T.cells.CD4.memory.resting" = "CD4 mem. rest.",
         "Cibersort_rel|T.cells.CD8" = "CD8 T",
         "Bindea|Macrophages" = "macrophages",
         "Cibersort_rel|Macrophages.M0" = "M0 macrophage",
         "Wolf|TREM1_data" = "TREM1",
         "Wolf|CHANG_CORE_SERUM_RESPONSE_UP" = "serum response",
         "Wolf|TAMsurr_TcClassII_ratio" = "TAM : T ratio")
t[, lab := LAB[feature]]
# two down-side labels want the same crowded spot; nudge them apart explicitly
t[, `:=`(nx = 0, ny = 0)]
t[feature == "Wolf|STAT1_score", `:=`(nx = -0.06, ny = -0.6)]
t[feature == "Attractors|G_CD3E", `:=`(nx = 0.10, ny = 1.2)]
t[feature == "Bindea|Macrophages", `:=`(nx = 0.02, ny = 1.4)]
t[feature == "Cibersort_rel|T.cells.CD8", `:=`(nx = -0.10, ny = -0.9)]

# the significance axis is square-root scaled: a handful of features run to P ~ 1e-27 while the
# informative cloud sits between P 1e-1 and 1e-5, which a linear -log10 axis squeezes into a strip.
# square canvas; POINT fill = significance tier; LABEL BOX fill = cell-type family;
# dashed verticals at the smallest |effect| that clears the empirical FDR, dashed horizontal at that FDR.
YMAX <- 29
XL <- max(abs(t$beta_compositional)) * 1.22
SIGLV <- c("P < 0.0001", "P < 0.001", "P < 0.01", "P < 0.05", "n.s.")
SIGCOL <- c("P < 0.0001" = "#F2D64B", "P < 0.001" = "#9CCB6A", "P < 0.01" = "#4FA97F",
            "P < 0.05" = "#2E7C86", "n.s." = "#D9D9D9")
t[, sigtier := factor(fifelse(!sig, "n.s.",
                       fifelse(p_compositional < 1e-4, "P < 0.0001",
                        fifelse(p_compositional < 1e-3, "P < 0.001",
                         fifelse(p_compositional < 1e-2, "P < 0.01", "P < 0.05")))), levels = SIGLV)]
XCUT <- st$effect_cut_fdr05                     # the smallest |effect| that clears the empirical FDR
lb <- t[!is.na(lab)]
p <- ggplot(t, aes(beta_compositional, y)) +
  geom_vline(xintercept = c(-XCUT, XCUT), linetype = "22", linewidth = LW, colour = "grey35") +
  geom_hline(yintercept = -log10(PCUT), linetype = "22", linewidth = LW, colour = "grey35") +
  geom_point(aes(fill = sigtier), shape = 21, size = 1.15, stroke = 0.16, colour = "grey25") +
  geom_label_repel(data = lb[beta_compositional < 0], aes(label = lab, fill = fam),
                   colour = "white", size = 6 / .pt, family = FONT, label.size = 0,
                   label.padding = unit(0.5, "mm"), label.r = unit(0.5, "mm"),
                   segment.size = 0.22, segment.colour = "grey45", min.segment.length = 0,
                   box.padding = 0.35, point.padding = 0.3, max.overlaps = Inf,
                   force = 8, force_pull = 0.3, xlim = c(-XL, -0.02), ylim = c(0.3, 17),
                   nudge_x = lb[beta_compositional < 0, nx], nudge_y = lb[beta_compositional < 0, ny],
                   seed = 3, max.iter = 4e4, show.legend = FALSE) +
  geom_label_repel(data = lb[beta_compositional > 0], aes(label = lab, fill = fam),
                   colour = "white", size = 6 / .pt, family = FONT, label.size = 0,
                   label.padding = unit(0.5, "mm"), label.r = unit(0.5, "mm"),
                   segment.size = 0.22, segment.colour = "grey45", min.segment.length = 0,
                   box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf,
                   force = 14, force_pull = 0.15, xlim = c(0.02, XL), ylim = c(0.3, 21),
                   nudge_x = lb[beta_compositional > 0, nx], nudge_y = lb[beta_compositional > 0, ny],
                   seed = 5, max.iter = 4e4, show.legend = FALSE) +
  scale_fill_manual(values = c(SIGCOL, setNames(FAMCOL, SHORT[names(FAMCOL)])),
                    breaks = SIGLV, name = NULL) +
  scale_x_continuous(limits = c(-XL, XL), breaks = seq(-0.3, 0.3, 0.15), expand = c(0, 0)) +
  scale_y_continuous(transform = "sqrt", limits = c(0, YMAX), breaks = c(0, 2, 5, 10, 20), expand = c(0, 0)) +
  labs(x = "WGD effect, ploidy-free (SD)",
       y = expression(-log[10] * italic(P)),
       title = "Composition of the immune compartment",
       caption = "← increased in WGD−                    increased in WGD+ →",
       subtitle = paste(sprintf("%d published immune features; %d carcinomas, proliferation-adjusted", st$n_features, st$n_cancer_types),
                        sprintf("dashed: empirical FDR 0.05 (%s permutations)", format(st$n_permutations, big.mark = ",")), sep = "\n")) +
  th + theme(plot.caption = element_text(size = TICK, colour = "grey25", hjust = 0.5, margin = margin(0.5, 0, 0, 0)),
             plot.caption.position = "plot",
             legend.position = "right", legend.direction = "vertical",
             legend.key.spacing.y = unit(0.2, "mm")) +
  guides(fill = guide_legend(override.aes = list(size = 1.4, shape = 21, stroke = 0.16, colour = "grey25"),
                             keyheight = unit(2.6, "mm")))

ggsave(file.path(OUT, "fig5b.png"), p, width = 76, height = 72, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig5b.pdf"), p, width = 76, height = 72, units = "mm", device = cairo_pdf, bg = "white", symbolfamily = SYMFAM)
