# TOOLS
require(here)
source(here::here("R/tools.R"))
require(ggpubr)

# DATA
source(here::here("R/DAT_prepare.R"))
dr <- d[bird_ID %in% d[duplicated(bird_ID), bird_ID]]
drw <- reshape(dr[, .(month, bird_ID, VAP, VSL, VCL, motileCount, Morph, age)], idvar = c("bird_ID", "Morph", "age"), timevar = "month", direction = "wide")

# PLOT
g1 <-
  ggplot(drw, aes(x = VCL.May, y = VCL.June)) +
  stat_smooth(method = MASS::rlm, col = 'grey30') +
  geom_point(pch = 21, alpha = .75, aes(fill = Morph, col = Morph)) +
  stat_cor(method = "pearson", size = 2, cor.coef.name = "r", aes(label = ..r.label..)) +
  geom_abline(slope = 1, col = "red", lty = 3) +
  ggtitle("Curvilinear") +
  ylab("Velocity in June [μm/s]") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = fills) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 7, color = "grey20", hjust = 0.5),
    axis.title.x = element_text(size = 8, color = "white"),
    axis.title.y = element_text(size = 8, colour = "grey10"),
    axis.text = element_text(size = 7),
    axis.text.y = element_text(margin = margin(r = -1)),
    axis.text.x = element_text(margin = margin(b = -1)),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.spacing = unit(-0.2, "mm"),
    panel.border = element_rect(color = "grey70")
  )

g2 <-
  ggplot(drw, aes(x = VSL.May, y = VSL.June)) +
  stat_smooth(method = MASS::rlm, col = 'grey30') +
  geom_point(pch = 21, alpha = .75, aes(fill = Morph, col = Morph)) +
  stat_cor(method = "pearson", size = 2, cor.coef.name = "r", aes(label = ..r.label..)) +
  geom_abline(slope = 1, col = "red", lty = 3) +
  ggtitle("Straight line") +
  xlab("Velocity in May [μm/s]") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = fills) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 7, color = "grey20", hjust = 0.5),
    axis.title.x = element_text(size = 8, colour = "grey10", hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 7, margin = margin(r = -1)),
    axis.text.x = element_text(size = 7, margin = margin(b = -1)),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.spacing = unit(-0.2, "mm"),
    panel.border = element_rect(color = "grey70")
  )

g3 <-
  ggplot(drw, aes(x = VAP.May, y = VAP.June)) +
  stat_smooth(method = MASS::rlm, col = 'grey30') +
  geom_point(pch = 21, alpha = .75, aes(fill = Morph, col = Morph)) +
  stat_cor(method = "pearson", size = 2, cor.coef.name = "r", aes(label = ..r.label..)) +
  geom_abline(slope = 1, col = "red", lty = 3) +
  ggtitle("Average path") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = fills) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 7, color = "grey20", , hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 7, margin = margin(b = -1)),
    axis.text.y = element_text(size = 7, margin = margin(r = -1)),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.spacing = unit(-0.2, "mm"),
    panel.border = element_rect(color = "grey70")
  )

gg1 <- ggplotGrob(g1)
gg2 <- ggplotGrob(g2)
gg3 <- ggplotGrob(g3)

gp_ind <- ggscatter(data.frame(x = 1, y = 1), x = "x", y = "y", shape = 21, color = cols[1], fill = ind) +
  theme_transparent() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))
gp_sat <- ggscatter(data.frame(x = 1, y = 1), x = "x", y = "y", shape = 21, color = cols[2], fill = sat) +
  theme_transparent() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))
gp_fae <- ggscatter(data.frame(x = 1, y = 1), x = "x", y = "y", shape = 21, color = cols[3], fill = fae) +
  theme_transparent() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))
gp_ind_grob <- ggplotGrob(gp_ind)
gp_sat_grob <- ggplotGrob(gp_sat)
gp_fae_grob <- ggplotGrob(gp_fae)

g_all <-
  ggarrange(
    g1,
    g2,
    g3,
    ncol = 3, align = "h", widths=c(1.11, 1, 1)
  ) + theme(plot.margin = margin(0, 2, 0, 0, "cm"))

#g_all = cbind(gg1, gg2, gg3, size = "first") + theme(plot.margin = margin(0, 0, 0.8, 0, "cm"))
ymin_ <- 0.3 # 0.55, 0.6
xmin_ <- .47 + 0.54
xmax_ <- .52 + 0.54

g_anot <-
  g_all +
  annotation_custom(gp_ind_grob, xmin = xmin_, xmax = xmax_, ymin = ymin_) +
  annotation_custom(gi, xmin = xmin_, xmax = xmax_, ymin = ymin_ + 0.22) +
  annotation_custom(gp_sat_grob, xmin = xmin_ + .06, xmax = xmax_ + .06, ymin = ymin_) +
  annotation_custom(gs, xmin = xmin_ + .06, xmax = xmax_ + .06, ymin = ymin_ + 0.22) +
  annotation_custom(gp_fae_grob, xmin = xmin_ + .12, xmax = xmax_ + .12, ymin = ymin_) +
  annotation_custom(gf, xmin = xmin_ + .12, xmax = xmax_ + .12, ymin = ymin_ + 0.22)

ggsave(here::here("Outputs/Fig_S4_width-120mm_v2.png"), g_anot, width = 9 / (5 / 7), height = 8 / 2, units = "cm") # #cbind(gg1, gg2, gg3, size = "first")

# END