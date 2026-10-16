# Fig. 3c: timing of genomic events relative to WGD in WGD+ carcinomas, one row per event. Bars: clonal mutations and arm
# losses timed before (blue) or after (gold) WGD in MSK-IMPACT, the blue length being P(pre-WGD | clonal); arm gains: PCAWG
# three-state posteriors (before / spanning / after WGD). Triangles: the same quantity in TCGA; squares: PCAWG. The coloured
# square at the name is the direction of the event's WGD association (CO / ME). Right: prevalence of the event in the cohort
# it is timed in. Solid lines: cohort-wide baselines over all non-synonymous mutations and over all timed arm losses. Gains
# timed in few patients (low_n in the table) are named in grey.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
DAT <- "figure3/data"; OUT <- "results/figure3"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); FONT <- if (any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "Liberation Sans"
ev <- read.delim(file.path(DAT, "fig3c_event_timing.tsv"), stringsAsFactors = FALSE)
bl0 <- read.delim(file.path(DAT, "fig3c_baselines.tsv"), stringsAsFactors = FALSE)
S <- ev[ev$block == "SNV drivers", ]; L <- ev[ev$block == "Arm losses", ]; G <- ev[ev$block == "Arm gains", ]
hdr <- data.frame(event = c("SNV drivers", " ", "Arm losses", "  ", "Arm gains"), stringsAsFactors = FALSE)
lev <- c("SNV drivers", S$event, " ", "Arm losses", L$event, "  ", "Arm gains", G$event)
A <- bind_rows(ev, hdr); A$event <- factor(A$event, levels = rev(lev)); is_hdr <- levels(A$event) %in% hdr$event[c(1, 3, 5)]
# bar segments as given in the table: clonal-only fractions for mutations and losses; the three states for gains
seg <- A %>% filter(!is.na(block)) %>% select(event, `pre-WGD` = frac_pre, straddling = frac_spanning, `post-WGD` = frac_post) %>%
  pivot_longer(-event, names_to = "cat", values_to = "frac")
seg$cat <- factor(seg$cat, levels = c("pre-WGD", "straddling", "post-WGD"))
cols <- c(`pre-WGD` = "#6E8FB3", straddling = "#D9D9D9", `post-WGD` = "#C9A66B")
tp <- ev %>% filter(!is.na(P_pre_tcga)) %>% transmute(event, x = P_pre_tcga)
pq <- ev %>% filter(!is.na(P_pre_pcawg)) %>% transmute(event, x = P_pre_pcawg)
tag <- ev %>% filter(!is.na(direction)) %>% transmute(event, direction)
XP0 <- 1.07; XPW <- 0.36
yof <- function(e) match(e, levels(A$event))
bl <- data.frame(x = bl0$baseline[match(c("SNV drivers", "Arm losses"), bl0$block)],
                 y0 = c(min(yof(S$event)), min(yof(L$event))) - 0.45,
                 y1 = c(max(yof(S$event)), max(yof(L$event))) + 0.45)
pv <- ev %>% transmute(event, block, prev = prevalence, low_n) %>%
  mutate(x0 = XP0, x1 = XP0 + XPW * prev, ncol = ifelse(low_n, "grey62", "grey20"))
pv$yi <- as.numeric(factor(pv$event, levels = levels(A$event))); stopifnot(!anyNA(pv$yi), !anyNA(pv$prev))
lab_col <- ifelse(is_hdr, "grey20", ifelse(levels(A$event) %in% pv$event[pv$ncol == "grey62"], "grey55", "grey10"))
g <- ggplot(seg, aes(y = event, x = frac, fill = cat)) +
  geom_col(width = 0.74, colour = "white", linewidth = 0.2, position = position_stack(reverse = TRUE)) +
  geom_point(data = tp, aes(y = event, x = x, shape = "TCGA"), inherit.aes = FALSE, size = 1.6, fill = "white", colour = "grey15", stroke = 0.45) +
  geom_point(data = pq, aes(y = event, x = x, shape = "PCAWG"), inherit.aes = FALSE, size = 1.5, fill = "white", colour = "grey15", stroke = 0.45) +
  geom_point(data = tag, aes(y = event, x = -0.04, colour = direction), inherit.aes = FALSE, shape = 15, size = 1.4) +
  geom_rect(data = pv, aes(xmin = x0, xmax = x1, ymin = yi - 0.37, ymax = yi + 0.37), inherit.aes = FALSE, fill = "grey72", colour = NA) +
  geom_segment(data = bl, aes(x = x, xend = x, y = y0, yend = y1, linetype = "baseline"), inherit.aes = FALSE, colour = "grey25", linewidth = 0.45) +
  annotate("segment", x = XP0, xend = XP0 + XPW, y = length(lev) + 0.75, yend = length(lev) + 0.75, colour = "grey30", linewidth = 0.35) +
  annotate("segment", x = c(XP0, XP0 + XPW / 2, XP0 + XPW), xend = c(XP0, XP0 + XPW / 2, XP0 + XPW), y = length(lev) + 0.75, yend = length(lev) + 1.0, colour = "grey30", linewidth = 0.35) +
  annotate("text", x = c(XP0, XP0 + XPW / 2, XP0 + XPW), y = length(lev) + 1.45, label = c("0", "50", "100"), family = FONT, size = 1.9, colour = "grey20") +
  annotate("text", x = XP0 + XPW / 2, y = length(lev) + 2.7, label = "prevalence (%)", family = FONT, size = 1.9, colour = "grey20", hjust = 0.5, vjust = 0.5, fontface = "bold") +
  geom_vline(xintercept = 0.5, colour = "grey55", linewidth = 0.4, linetype = "22") +
  scale_colour_manual(values = c(CO = "#3A9E6F", ME = "#6A4C9C"), name = NULL, labels = c(CO = "CO", ME = "ME")) +
  scale_shape_manual(values = c(TCGA = 24, PCAWG = 22), name = NULL) +
  scale_fill_manual(values = cols, name = NULL, labels = c("pre-WGD", "spans WGD", "post-WGD")) +
  scale_linetype_manual(values = c("baseline" = "solid"), name = NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), limits = c(-0.07, 1.45), expand = c(0, 0),
                     name = "fraction of patients") +
  scale_y_discrete(drop = FALSE, expand = expansion(add = c(0.6, 3.5))) +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 2, nrow = 1), colour = guide_legend(order = 3, nrow = 1, override.aes = list(size = 2)), linetype = guide_legend(order = 4, override.aes = list(colour = "grey25", linewidth = 0.6))) +
  theme_classic(base_family = FONT, base_size = 9) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 5.3, colour = lab_col, face = ifelse(is_hdr, "bold", "plain")),
        axis.line.x = element_line(linewidth = 0.45), axis.ticks.x = element_line(linewidth = 0.45), axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 7.2), legend.position = "bottom", legend.box = "vertical", legend.justification = c(0.5, 0), legend.box.just = "center", legend.text = element_text(size = 5.3),
        legend.key.size = unit(2.4, "mm"), legend.title = element_text(size = 6.2), legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0, "mm"), legend.box.spacing = unit(1, "mm"),
        plot.margin = margin(3, 3, 1, 1))
ggsave(file.path(OUT, "fig3c.png"), g, width = 85, height = 99, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(OUT, "fig3c.pdf"), g, width = 85, height = 99, units = "mm", device = cairo_pdf, bg = "white")
