# tools
   require(here)
   source(here::here("R/tools.R"))
   require(ggpubr)


# data
  source(here::here("R/DAT_prepare.R"))   
  dl <- data.table(melt(d[, .(bird_ID, month, Morph, VAP, VSL, VCL)], id.vars = c("bird_ID", "month", "Morph"), variable.name = "mot"))
  dl[mot == "VAP", mot := "Average path"]
  dl[mot == "VSL", mot := "Straight line"]
  dl[mot == "VCL", mot := "Curvilinear"]
  dl[, mot := factor(mot, levels = c("Curvilinear", "Straight line", "Average path"))]

# plot
  g <- ggplot(dl, aes(x = Morph, y = value)) +
      facet_wrap(~mot, scales = "free_y", nrow = 3, strip.position = "right") +
      geom_dotplot(
          binaxis = "y", stackdir = "center",
          position = position_dodge(), aes(fill = Morph, col = Morph)
      ) +
      geom_boxplot(col = "grey50", fill = NA, outlier.shape = NA) +
      scale_color_manual(values = c(cols, "#a53708")) +
      scale_fill_manual(values = c(ind, sat, fae, "#f89f79")) +
      scale_y_continuous("Velocity [μm/s]", expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      theme_bw() +
      theme(
          axis.ticks = element_blank(),
          axis.title = element_text(size = 10, colour = "grey10"),
          strip.text.y.right = element_text(color = "grey20", margin = margin(1, 1, 1, 1, "mm"), angle = 90),
          strip.background = element_rect(fill = NA, colour = NA, size = 0.25),
          panel.border = element_rect(color = "grey70"),
          plot.margin = margin(14, 3, 1, 1, "mm"),
          legend.position = "none"
          # axis.title.y = element_text(size = 8)
      )
  gg <- ggarrange(g)
  ggExp <- gg +
      annotation_custom(gi, xmin = 0.125+0.025, xmax = 0.275+0.025, ymin = 0.91) +
      annotation_custom(gs, xmin = 0.125 + 0.215, xmax = 0.275 + 0.215, ymin = 0.91) +
      annotation_custom(gf, xmin = 0.125 + 0.4, xmax = 0.275 + 0.4, ymin = 0.906) +
      annotation_custom(gz, xmin = 0.125 + 0.6, xmax = 0.275 + 0.6, ymin = 0.88)
  
  ggsave("Outputs/Fig_S1_width-60mm.png", ggExp, width = 6 / (5 / 7), height = 13, units = "cm", bg = "white", dpi = 600)

  # END