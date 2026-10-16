# Fig. 2d: cancer-type-specific versus pan-cancer effect size of every retained alteration. x = |log2 ER| in the most
# significant retained cancer type, y = pan-cancer |log2 ER| (each the mean of the two cohorts); fill = where the alteration
# is retained, label colour = direction (green CO / purple ME); cancer-type-only alterations carry their cancer type.
# |log2 ER| is capped at 3 for display; for alterations retained pan-cancer only, x is the most
# significant cancer type tested in both cohorts; when one cohort has no estimate the other cohort's value is used.
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })
DAT <- "figure2/data"; OUT <- "results/figure2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
CO <- "#3A9E6F"; ME <- "#6A4C9C"; CAT <- c(both = "#2C6EA6", `pancancer only` = "#A6A6A6", `type-specific only` = "white")
s <- read.delim(file.path(DAT, "fig2d_effect_sizes.tsv"))
disp <- function(d) { x <- ifelse(is.na(d$display_MSK) | d$display_MSK == "", d$display_TCGA, d$display_MSK)   # focal LOH: cytoband (gene)
  ifelse(d$class == "focal_loh" & !is.na(x) & x != "", x, d$unit) }
lk <- read.delim(file.path(DAT, "fig2_scope_pairs.tsv"), check.names = FALSE); lk <- lk[order(lk$scope != "pancancer"), ]
lk <- lk[!duplicated(lk[c("class", "unit")]), ]; lk$disp <- disp(lk)
s$unit <- lk$disp[match(paste(s$class, s$unit), paste(lk$class, lk$unit))]
cls <- c(broad_gain = "gain", broad_loss = "loss", broad_loh = "LOH", focal_amp = "amp", focal_del = "del", focal_loh = "LOH", snv = "mut")
s$label <- paste(s$unit, cls[s$class])
s$label <- ifelse(s$category_2d == "type-specific only", sprintf("%s (%s)", s$label, s$best_type), s$label)
s$category_2d <- factor(s$category_2d, levels = names(CAT))
stopifnot(all(s$direction %in% c("CO", "ME"))); s$dir <- s$direction
anchors <- c("TP53 mut", "CCNE1 amp", "1q gain", "17p loss")                 # the alterations named in the text
s$score <- s$pan_coord + s$type_coord
cl <- s$type_coord < 0.7 & s$pan_coord < 0.45                                   # the dense lower-left cluster
lab_cl <- s[cl & s$label %in% anchors, ]                                        # only anchors are labelled inside the cluster
lab_out <- s[!cl & rank(-s$score) <= 14, ]                                      # strongest alterations outside the cluster
lab_cl <- lab_cl[order(-lab_cl$type_coord), ]
lab_cl$lx <- 0.03; lab_cl$ly <- seq(0.74, by = -0.04, length.out = nrow(lab_cl))
lab_cl$lx_end <- c(`TP53 mut` = 0.27, `CCNE1 amp` = 0.27, `1q gain` = 0.29, `17p loss` = 0.245)[lab_cl$label]   # leader exits right of the text
p <- ggplot(s, aes(type_coord, pan_coord)) +
  geom_point(aes(fill = category_2d), shape = 21, size = 2.0, stroke = 0.35, colour = "grey20") +
  geom_segment(data = lab_cl, aes(x = lx_end, y = ly, xend = type_coord, yend = pan_coord, colour = dir), linewidth = 0.3, alpha = 0.85) +
  geom_text(data = lab_cl, aes(x = lx, y = ly, label = label, colour = dir), hjust = 0, family = FONT, size = 2.5) +
  geom_text_repel(data = lab_out, aes(label = label, colour = dir, segment.colour = dir),
                  family = FONT, size = 2.5, seed = 4, xlim = c(0.72, NA), ylim = c(0.03, 0.95),
                  box.padding = 0.4, point.padding = 0.25, min.segment.length = 0, segment.size = 0.3, segment.alpha = 0.85,
                  force = 2, max.overlaps = Inf, max.iter = 30000, max.time = Inf) +
  scale_fill_manual(values = CAT, name = "Retained in", labels = c("pan-cancer and cancer type", "pan-cancer only", "cancer type only")) +
  scale_colour_manual(values = c(CO = CO, ME = ME), guide = "none", aesthetics = c("colour", "segment.colour")) +
  labs(x = "Cancer-type-specific effect, |log2 ER|", y = "Pan-cancer effect, |log2 ER|") +
  theme_classic(base_family = FONT, base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.45), axis.ticks = element_line(linewidth = 0.45), axis.text = element_text(size = 8, colour = "grey10"),
        legend.position = "bottom", legend.text = element_text(size = 7.6), legend.title = element_text(size = 7.6), legend.key.size = unit(3.2, "mm"),
        plot.margin = margin(4, 6, 2, 4))
ggsave(file.path(OUT, "fig2d.png"), p, width = 130, height = 118, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig2d.pdf"), p, width = 130, height = 118, units = "mm", device = cairo_pdf, bg = "white")
