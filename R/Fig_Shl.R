# TOOLS 
  require(here)
  source(here::here('R/tools.R'))
  
  require(ggpubr)
  require(ggsci) 
  require(grid)
  require(gridExtra)
  require(magrittr)
  require(MASS)
  require(PerformanceAnalytics)
  require(stringi)
  require(viridis)

  # constants
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)

    fae = '#d4b691' # 'ffd6af'
    sat = 'white'
    ind = '#303030'
  
  # functions
    getime = function (x) {ifelse(is.na(x), as.numeric(NA), as.numeric(difftime(x, trunc(x,"day"), units = "hours")))}
    
    getDay = function (x) {as.Date(trunc(x, "day"))}
   
   # custom ggplot theme
       theme_MB = theme(  
                title = element_text(size=8, colour="grey30"),
                axis.line = element_blank(),
                #axis.line = element_line(colour="grey70", size=0.25),
                axis.title = element_text(size=7, colour="grey30"),
                axis.title.y = element_text(vjust=3.5),
                axis.title.x = element_text(vjust=1),
                axis.text = element_text(size=6),#, vjust = 0.5, hjust=1),# margin=units(0.5,"mm")),
                axis.ticks.length=unit(0.5,"mm"),
                axis.ticks = element_line(colour = "grey70", size = 0.1),
                #axis.ticks.margin,
                
                strip.text.x = element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
                strip.text.y = element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
                strip.background = element_rect(fill="grey99",colour="grey70", size=0.25),
                  #strip.background = element_blank(), 
                  #strip.text = element_blank(),
                panel.spacing = unit(0, "mm"),
                panel.background=element_blank(),
                panel.border = element_rect(colour="grey70", size=0.1, fill = NA), #panel.border=element_blank(),
                panel.grid = element_blank(),

                legend.text=element_text(size=6),
                legend.title=element_text(size=6),
                legend.key = element_rect(colour = NA, fill = NA),
                legend.key.height= unit(0.5,"line"),
                legend.key.width = unit(0.25, "cm"),
                legend.margin = margin(0,0,0,0, unit="cm"),
                legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
                legend.background = element_blank()
                )  
# DATA 
  source(here::here('R/DAT_prepare.R'))       
 # prepare motility for plotting
    am = a[part=='Acrosome']
    aml = data.table(melt(am[,.(bird_ID,month,Morph,age,HL,motileCount,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph","age","HL","motileCount"), variable.name = "Motility"))
    aml[Motility == 'VAP' ,mot:='Average path']
    aml[Motility == 'VCL' ,mot:='Curvilinear']
    aml[Motility == 'VSL' ,mot:='Straight line']

# prepare panels
  # dummie dataset to create empty panels
    dum = aml[1:4] 
    dum[,mot2 := c('w','x','y','z')]
    aml[, mot2 := factor(mot,levels=c('Curvilinear','Straight line','Average path','w','x','y','z'))]
    amld = rbind(aml,dum)

  gv =
  ggplot(aml, aes(x = HL, y = value)) +
    stat_smooth(method = MASS::rlm, col ='grey30') +
    geom_point(pch = 21, alpha = 0.75, aes(fill = Morph, col = Morph))+
    #stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    facet_wrap(~mot2, scales = 'free_y', nrow = 1,drop=FALSE)+
  
    scale_color_manual(values = cols)+
    scale_fill_manual(values = fills)+
    
    xlab('Homozygousity by locus') +
    ylab('Velocity [μm/s]') +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size=8),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 8,colour="grey10"),
      axis.title.x = element_blank(),
      axis.text = element_text(size=7), 
      axis.text.x = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_text(size = 7, color="grey20",  margin=margin(1,1,1,1,"mm")),
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
      
      panel.border = element_rect(color = 'grey70')
      )  

  ggv = ggplotGrob(gv)
  rm_grobs <- ggv$layout$name %in% c("panel-4-1", "panel-5-1","panel-6-1","panel-7-1", "strip-t-4-1", "strip-t-5-1","strip-t-6-1", "strip-t-7-1")

  ggv$grobs[rm_grobs] <- NULL
  ggv$layout <- ggv$layout[!rm_grobs, ]

  gm =
  ggplot(a, aes(x = HL, y = Length_avg)) +
    stat_smooth(method = MASS::rlm, col ='grey30') +
    geom_point(pch = 21, alpha = 0.75, aes(fill = Morph, col = Morph))+
    #stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    facet_wrap(~part, scales = 'free_y', nrow = 1,drop=FALSE)+
  
    scale_color_manual(values = cols)+
    scale_fill_manual(values = fills)+
    
    xlab('Homozygousity by locus') +
    ylab('Length [µm]') +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size=8),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 8,colour="grey10"),
      axis.title.x = element_blank(),
      axis.text = element_text(size=7), 
      axis.text.x = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      
      strip.text = element_text(size = 7, color="grey20",  margin=margin(1,1,1,1,"mm")),
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
      
      panel.border = element_rect(color = 'grey70')
      )  
  
  gcv = 
  ggplot(cv_, aes(x = HL, y = CV)) +
    stat_smooth(method = MASS::rlm, col ='grey30') +
    geom_point(pch = 21, alpha = 0.75, aes(fill = Morph, col = Morph))+
    #stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    facet_wrap(~part, scales = 'free_y', nrow = 1,drop=FALSE)+
  
    scale_color_manual(values = cols)+
    scale_fill_manual(values = fills)+
    
    xlab('Homozygousity by locus') +
    ylab('Coefficient of variation') +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size=8),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 8,colour="grey10"),
      axis.text = element_text(size=7), 
      axis.text.x = element_text(margin = margin(t = -1)),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_blank(),
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
      
      panel.border = element_rect(color = 'grey70')
      )  

  ggHL = ggarrange(
    ggv,
    gm,
    gcv,  
    nrow=3, align = 'v' #heights=c(1.5, 4, 4.8),  
    )   
  
# prepare legend
    gp_ind = ggscatter(data.frame(x =1, y =1), x = 'x', y = 'y', shape = 21, color =cols[1], fill =ind) +
      theme_transparent()+
      theme(plot.margin = unit(c(0,0,0,0), "mm"))
    gp_sat = ggscatter(data.frame(x =1, y =1), x = 'x', y = 'y', shape = 21, color =cols[2], fill =sat) +
      theme_transparent()+
      theme(plot.margin = unit(c(0,0,0,0), "mm"))
    gp_fae = ggscatter(data.frame(x =1, y =1), x = 'x', y = 'y', shape = 21, color =cols[3], fill =fae) +
      theme_transparent()+
      theme(plot.margin = unit(c(0,0,0,0), "mm")) 
    gp_ind_grob = ggplotGrob(gp_ind)
    gp_sat_grob = ggplotGrob(gp_sat)
    gp_fae_grob = ggplotGrob(gp_fae)

# combine & export
  g_exp =   
    ggHL + 
    annotation_custom(gp_ind_grob, xmin=.47, xmax=.52, ymin = 0.5)+
    annotation_custom(gi, xmin=.47, xmax=.52, ymin = 0.7) +
    annotation_custom(gp_sat_grob, xmin=.47+.06, xmax=.52+.06, ymin = 0.5) +
    annotation_custom(gs, xmin=.47+.06, xmax=.52+.06, ymin = 0.7) +
    annotation_custom(gp_fae_grob, xmin=.47+.12, xmax=.52+.12, ymin = 0.5) +
    annotation_custom(gf, xmin=.47+.12, xmax=.52+.12, ymin = 0.7) 

  ggsave('Outputs/Fig_Shl_140mm.png',g_exp, width = 14/(5/7), height =9, units = 'cm', bg="white", dpi = 600)


# END      