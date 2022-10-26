# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  require(ggnewscale)
  require(ggpubr) 
  require(ggsci)
  require(MuMIn)

  width_ = 0.6 # spacing between error bars
  width__ = 0.4 
  col_ = c(viridis(1, alpha = 1, begin = 0.2, end = 0.2, direction = 1, option = "D"),
         viridis(1, alpha = 1, begin = 0.5, end = 0.5, direction = 1, option = "D"),
         viridis(1, alpha = 1, begin = 0.9, end = 0.9, direction = 1, option = "D")
         )

# DATA
  source(here::here('R/DAT_prepare.R'))       

  dd = d[!Morph%in%'Zebra finch']
    
  # long format
    ddl = data.table(melt(dd[,.(bird_ID,month,Morph,age,motileCount,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph","age","motileCount"), variable.name = "Motility"))
    ddl[Motility == 'VAP' ,mot:='Average path']
    ddl[Motility == 'VCL' ,mot:='Curvilinear']
    ddl[Motility == 'VSL' ,mot:='Straight line']
    ddl[, animal := bird_ID]
    

   # use June values and for 4 males without June, May
      dd1 = dd[month == 'June']
      dd2 = dd[month == 'May']
      ddx = rbind(dd1,dd2[!bird_ID%in%dd1$bird_ID])
      ddxl = data.table(melt(ddx[,.(bird_ID,month,Morph,age,motileCount,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph","age","motileCount"), variable.name = "Motility"))
      ddxl[Motility == 'VAP' ,mot:='Average path']
      ddxl[Motility == 'VCL' ,mot:='Curvilinear']
      ddxl[Motility == 'VSL' ,mot:='Straight line']
      ddxl[, animal := bird_ID]

      aw_ = aw[, c('bird_ID','Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')]

      adw = merge(ddx, aw_, by = 'bird_ID')

      adwl = merge(ddxl, aw_, by = 'bird_ID')
      adwl[Motility == 'VAP' ,Motility:='Average path']
      adwl[Motility == 'VCL' ,Motility:='Curvilinear']
      adwl[Motility == 'VSL' ,Motility:='Straight line']
      adwl[, Motility := factor(Motility, levels=(c("Curvilinear", "Straight line", "Average path")))] 

  adw[,motileCount_ln_z := scale(log(motileCount))]
  adwl[,motileCount_ln_z := scale(log(motileCount))]
  ar[,motileCount_ln_z := scale(log(motileCount))]
 

# Prepare estimates
  effects_xx = c('Head', 'Midpiece', 'Tail')
  l = list()
  for(i in unique(adwl$mot)){
    #i ='Curvilinear'
    ai = adwl[mot == i]
    ms = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Head)+scale(Midpiece)+scale(Tail), ai)
    bsim = sim(ms, n.sim=nsim) 
    v = c(apply(bsim@coef, 2, quantile, prob=c(0.5)))[5:7]
    lwr = c(apply(bsim@coef, 2, quantile, prob=c(0.025)))[5:7]
    upr = c(apply(bsim@coef, 2, quantile, prob=c(0.975)))[5:7]
    ii = data.table(motility=i, effect=effects_xx,estimate=v, lwr=lwr, upr=upr)

    mh = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Head), ai)
    bsim = sim(mh, n.sim=nsim) 
    ii[effect == 'Head', estimate_s := c(apply(bsim@coef, 2, quantile, prob=c(0.5)))[5]]
    ii[effect == 'Head', lwr_s := c(apply(bsim@coef, 2, quantile, prob=c(0.025)))[5]]
    ii[effect == 'Head', upr_s := c(apply(bsim@coef, 2, quantile, prob=c(0.975)))[5]]
    
    mm = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Midpiece), ai)
    bsim = sim(mm, n.sim=nsim) 
    ii[effect == 'Midpiece', estimate_s := c(apply(bsim@coef, 2, quantile, prob=c(0.5)))[5]]
    ii[effect == 'Midpiece', lwr_s := c(apply(bsim@coef, 2, quantile, prob=c(0.025)))[5]]
    ii[effect == 'Midpiece', upr_s := c(apply(bsim@coef, 2, quantile, prob=c(0.975)))[5]]
    
    mt = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Tail), ai)
    bsim = sim(mt, n.sim=nsim) 
    ii[effect == 'Tail', estimate_s := c(apply(bsim@coef, 2, quantile, prob=c(0.5)))[5]]
    ii[effect == 'Tail', lwr_s := c(apply(bsim@coef, 2, quantile, prob=c(0.025)))[5]]
    ii[effect == 'Tail', upr_s := c(apply(bsim@coef, 2, quantile, prob=c(0.975)))[5]]
    
    l[[i]]=ii
  }
  x = do.call(rbind, l)
  x[, motility := factor(motility, levels=rev(c("Curvilinear", "Straight line", "Average path")))] 
  x = x[order(motility)]
  y = x  
  y1 = y[,c('motility', 'effect', 'estimate', 'lwr', 'upr')]  
  y2 = y[,c('motility', 'effect', 'estimate_s', 'lwr_s', 'upr_s')] 
  setnames(y2, old = c('motility', 'effect', 'estimate_s', 'lwr_s', 'upr_s'), new = c('motility', 'effect', 'estimate', 'lwr', 'upr')) 
  y1[, Model:='Head + Midpiece + Tail']
  y2[, Model:='simple, single term']
  y = rbind(y1,y2)

# plot and export
  g = 
  ggplot(y, aes(y = effect, x = estimate, fill = motility, col = motility, shape = Model)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width = width_)) +
    
    scale_x_continuous(limits = c(-0.25, 0.35), expand = c(0, 0))+
    #scale_color_jco(name = 'Motility', guide = guide_legend(reverse = TRUE))+
    #scale_fill_jco(name = 'Motility', guide = guide_legend(reverse = TRUE))+
    scale_color_manual(name = 'Velocity', guide = guide_legend(reverse = TRUE, order = 1,nrow=3,byrow=TRUE), values = col_)+
    scale_fill_manual(name = 'Velocity', guide = guide_legend(reverse = TRUE, order = 1,nrow=3,byrow=TRUE), values = col_)+
    scale_shape_manual(name = 'Model', values =c(23,21), guide = guide_legend(reverse = TRUE, override.aes = list(fill = c('grey30'), col = 'grey30'),order = 0, nrow=2,byrow=TRUE))+
    labs(y = NULL, x = "Standardized effect size")+

    theme_bw() +
    theme(plot.subtitle = element_text(size=9, color = 'grey30'),
        legend.title=element_text(size=8, color = 'grey30'),
        legend.text=element_text(size=7.5, color = 'grey30'),
        #legend.spacing.y = unit(1, 'mm'),
        legend.key.height= unit(0.5,"line"),
        legend.margin=margin(0,0,0,0),
        #legend.position=c(0.45,1.6),

        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 10, color ='grey10'),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank(),

        plot.margin = margin(3,3,1,1, "mm")
        )  
  ggsave('Outputs/Fig_Sxx_viridis.png',g, width = 7/(5/7), height =4/(5/7), units = 'cm', bg="transparent", dpi = 600)

# END