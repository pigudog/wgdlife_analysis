# Fig. 3d: cancer-type-adjusted consensus of the event timing. Small points: one platform x cancer-type cell (>= 5 patients),
# P(pre-WGD | clonal) (circle MSK-IMPACT, triangle TCGA, square PCAWG). Diamond and bar: DerSimonian-Laird random-effects
# consensus with its 95% CI, one vote per cell; diamond fill = the CI excludes the WGD anchor (blue pre, gold post, white
# unresolved); stars = BH-adjusted test of the consensus against the anchor.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
DAT <- "figure3/data"; OUT <- "results/figure3"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
A <- read.delim(file.path(DAT, "fig3d_consensus.tsv"), stringsAsFactors = FALSE)
B <- read.delim(file.path(DAT, "fig3d_consensus_by_type.tsv"), stringsAsFactors = FALSE)
lab <- function(e) ifelse(grepl("^loss_", e), paste0(sub("^loss_", "", e), " loss"), e)
A$row <- lab(A$event); B$row <- lab(B$event)
S <- A %>% filter(evidence == "snv_driver") %>% arrange(desc(P_MSK_pooled))    # row order: MSK-IMPACT pooled estimate
L <- A %>% filter(evidence == "loss_arm") %>% arrange(desc(P_MSK_pooled))
lev <- c("SNV drivers", S$row, " ", "Arm losses", L$row)
A$row <- factor(A$row, levels = rev(lev)); B$row <- factor(B$row, levels = rev(lev))
A$star <- ifelse(A$q_anchor < 0.001, "***", ifelse(A$q_anchor < 0.01, "**", ifelse(A$q_anchor < 0.05, "*", "ns")))
hdr <- c("SNV drivers", "Arm losses"); is_hdr <- levels(A$row) %in% hdr
B$platform <- factor(B$platform, levels = c("MSK", "TCGA", "PCAWG"))
XI <- 1.08
g <- ggplot() +
  geom_vline(xintercept = 0.5, colour = "grey55", linewidth = 0.4, linetype = "22") +
  geom_point(data = B, aes(x = P_type, y = row, shape = platform), colour = "grey55", fill = "grey80", size = 0.7, stroke = 0.25, alpha = 0.45,
             position = position_jitter(width = 0, height = 0.14, seed = 1)) +
  geom_errorbar(data = A, aes(y = row, xmin = ci_lo, xmax = ci_hi), width = 0, orientation = "y", colour = "grey15", linewidth = 0.55) +
  geom_point(data = A, aes(x = P_consensus, y = row, fill = side), shape = 23, size = 2.0, colour = "grey15", stroke = 0.4) +
  geom_text(data = A, aes(x = XI, y = row, label = star), family = FONT, size = 2.2, colour = "grey20", hjust = 1) +
  annotate("text", x = XI, y = length(lev) + 1.0, label = "vs WGD", family = FONT, size = 2.0, colour = "grey20", hjust = 1, fontface = "bold") +
  scale_fill_manual(values = c(pre = "#6E8FB3", post = "#C9A66B", unresolved = "white"), name = NULL,
                    labels = c(pre = "pre-WGD", post = "post-WGD", unresolved = "unresolved"), breaks = c("pre", "post", "unresolved")) +
  scale_shape_manual(values = c(MSK = 21, TCGA = 24, PCAWG = 22), name = NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), limits = c(-0.02, 1.10), expand = c(0, 0),
                     name = "P(pre-WGD | clonal)") +
  scale_y_discrete(drop = FALSE, expand = expansion(add = c(0.6, 1.6))) +
  coord_cartesian(clip = "off") +
  guides(shape = guide_legend(order = 1, nrow = 1, override.aes = list(size = 1.8)), fill = guide_legend(order = 2, nrow = 1, override.aes = list(shape = 23, size = 1.8))) +
  theme_classic(base_family = FONT, base_size = 9) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 5.3, colour = ifelse(is_hdr, "grey20", "grey10"), face = ifelse(is_hdr, "bold", "plain")),
        axis.line.x = element_line(linewidth = 0.45), axis.ticks.x = element_line(linewidth = 0.45), axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 7.2), legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 5.5), legend.spacing.y = unit(0, "mm"),
        legend.key.size = unit(2.4, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
        plot.margin = margin(3, 3, 1, 1))
ggsave(file.path(OUT, "fig3d.png"), g, width = 73, height = 105, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig3d.pdf"), g, width = 73, height = 105, units = "mm", device = cairo_pdf, bg = "white")
