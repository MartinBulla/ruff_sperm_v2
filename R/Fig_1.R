# =============================================================
# ❗ Runs relative to the project's root directory, requires
# testes.xlsx, and exports Fig 1 into ./Outputs/
# =============================================================

# Tools
   require(here)
   source(here::here("R/tools.R"))
   require(ggpubr)

# Data   
   t = data.table(read_excel(here::here("Data/testes.xlsx"), sheet = 1)) # , range = "A1:G161"))
   t[, Morph := factor(Morph, levels = c("Res", "Sat", "Faed"))]
   t[Morph == "Res", Morph := "Independent"]
   t[Morph == "Sat", Morph := "Satellite"]
   t[Morph == "Faed", Morph := "Faeder"]
   t[, soma := Bodymass - Gonadmass]

# plot  
   g0 = ggplot(t, aes(x = Morph, y = Gonadmass)) + 
    geom_boxplot() + 
    annotation_custom(gi, xmin=0.725, xmax=1.3, ymin = 4.665) + 
    annotation_custom(gs, xmin=1.725, xmax=2.3, ymin = 4.665)+#, ymin=4.5, ymax=5.1) + 
    annotation_custom(gfc, xmin=2.745, xmax=3.255,ymin=4.655)+#, ymin=4.68, ymax=4.82) +
    scale_y_continuous('Testes mass [g]', limits = c(4.7,4.8), breaks = c(4.7,4.8))+#scale_y_continuous(limits = c(4.9,5.1))+
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_text(size = 10, color = 'transparent'), 
          axis.text.y = element_text(color = 'transparent'),
          axis.ticks = element_blank(),
          plot.margin = unit(c(1,1.75,0,2), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background  = element_rect(fill='white', color = 'white'),
          legend.position = "none"
          )
   
   g1 = ggplot(t, aes(x = Morph, y = Gonadmass)) + 
    #annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=4.5, ymax=5.1) + 
    #annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=4.5, ymax=5.1) + 
    #annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=4.47, ymax=4.9) +
    geom_boxplot(col = 'grey50', outlier.shape = NA) + 
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(),  aes(col = Morph, fill =Morph), dotsize = 1.1)+ #col = 'darkgrey',
    scale_color_manual(values=c('black','darkgrey','#bf925a'))+ #scale_fill_viridis(discrete=TRUE)+
    scale_fill_manual(values=c(ind,sat,fae))+
    #scale_color_manual(values=c(ind,sat,fae))+
    #scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Testes mass [g]', limits = c(1.5,4.5),expand = c(0, 0))+
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
      axis.title = element_text(size = 10, , colour="grey10"),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = 'grey70'),
      #axis.title.y = element_text(size = 8),
      legend.position = "none")

   g2 = ggplot(t, aes(x = Morph, y = Bodymass)) + 
    #annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=139, ymax=160) + 
    #annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=139, ymax=160) + 
    #annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=140, ymax=151) + 
    geom_boxplot(col = 'grey50', outlier.shape = NA) +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(),  aes(fill =Morph, col = Morph))+
    
    scale_color_manual(values=c('black','darkgrey','#bf925a'))+ #scale_fill_viridis(discrete=TRUE)+
    scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Body mass [g]', limits = c(100,200), expand = c(0, 0))+
    theme_bw()+
    theme(axis.ticks = element_blank(),
      axis.title = element_text(size = 10, colour="grey10"),
      panel.border = element_rect(color = 'grey70'),
      legend.position = "none"
      #axis.title.y = element_text(size = 8)
      )     

   ggA = ggarrange(
    g0,#+theme(plot.margin = unit(c(0,1.75,0,2), "mm")),
    g1+theme(plot.margin = unit(c(0,1.75,0.3,2), "mm")),
    g2+theme(plot.margin = unit(c(1.5,1.75,0.3,2), "mm")),  
    nrow=3, heights=c(1.5, 4, 4.8),  align = 'v'
    )  
   ggA         
   ggsave('Outputs/Fig_1_width-50mm.png',ggA, width = 7, height =13, units = 'cm', bg="white", dpi = 600)

# END