# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  require(ggnewscale)
  require(ggpubr) 
  require(ggsci)
  require(MuMIn)

  width_ = 0.6 # spacing between error bars

  col_ = c(viridis(1, alpha = 1, begin = 0.2, end = 0.2, direction = 1, option = "D"),
         viridis(1, alpha = 1, begin = 0.5, end = 0.5, direction = 1, option = "D"),
         viridis(1, alpha = 1, begin = 0.9, end = 0.9, direction = 1, option = "D")
         )
  #col_ = c(viridis(1, alpha = 1, begin = 0.3, end = 0.3, direction = 1, option = "D"),
  #      viridis(1, alpha = 1, begin = 0.6, end = 0.6, direction = 1, option = "D"),
  #       viridis(1, alpha = 1, begin = 0.9, end = 0.9, direction = 1, option = "D")
  #       )
  col_2 = rev(col_)
  #show_col(col_) 
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
 
# prepare estimates and pred for univariate models
  effects_ = c('intercept','motileCount_ln', 'morphSat', 'morphFae', 'pred')
  lvx = list()
  lvpx =list()
  for(i in c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')){
    #i = 'Acrosome'
    xi = adwl[, i, with = FALSE] 
    adwl[, Length_avg := xi] 
    for(k in unique(adwl$Motility)){
      #k = 'Average path'
      adwlk = adwl[Motility == k]
      # get estimates
      m = lm(scale(value) ~ scale(log(motileCount)) + Morph + scale(Length_avg), adwlk)

      bsim = sim(m, n.sim=nsim) 
      v = c(apply(bsim@coef, 2, quantile, prob=c(0.5)))
      lwr = c(apply(bsim@coef, 2, quantile, prob=c(0.025)))
      upr = c(apply(bsim@coef, 2, quantile, prob=c(0.975)))
      
      lvx[[paste(k,i)]]=data.frame(response=k, trait = i, effect=effects_,estimate=v, lwr=lwr, upr=upr)
      
      # get predictions
      m = lm(value ~ motileCount_ln_z + Morph + Length_avg, adwlk)
      bsim = sim(m, n.sim=nsim) 
      v = apply(bsim@coef, 2, quantile, prob=c(0.5))
      newD=data.frame(motileCount_ln_z = mean(adwlk$motileCount_ln_z), Morph = unique(b$Morph)[2], Length_avg = seq(min(adwlk$Length_avg), max(adwlk$Length_avg), length.out = 100)) # values to predict for
      X <- model.matrix(~ motileCount_ln_z + Morph + Length_avg,data=newD) # exactly the model which was used has to be specified here 
      newD$pred <-(X%*%v) 
      predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
      for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                  predmatrix[predmatrix < 0] <- 0
                  newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                  newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                  #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
      newD$response = k
      newD$trait = i
      lvpx[[paste(k,i)]] = data.table(newD)
      print(paste(k,i))   
      }
    }
         
  llvx = data.table(do.call(rbind,lvx) ) 
  llvx[, Motility:=response]
  llvx[, Motility := factor(Motility, levels=rev(c("Curvilinear", "Straight line", "Average path")))] 
  llvx[, trait := factor(trait, levels=rev(c(c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel'))))] 

  llvpx = data.table(do.call(rbind,lvpx) ) 
  llvpx[ , Motility:=response]
  llvpx[, Motility := factor(Motility, levels=(c("Curvilinear", "Straight line", "Average path")))] 
  llvpx[, trait := factor(trait, levels=rev(c(c('Acrosome', 'Nucleus', 'Midpiece', 'Tail', 'Total', 'Head','Flagellum', 'Midpiece_rel', 'Flagellum_rel'))))] 

# (a) & (b) labels top
# (a) - effect sizes
  llvx_ = llvx[effect == 'pred']
  llvx_[trait == 'Midpiece_rel', trait:='Midpiece\n(relative)']
  llvx_[trait == 'Flagellum_rel', trait:='Flagellum\n(relative)']
  llvx_[, trait := factor(trait, levels=rev(c("Acrosome", "Nucleus", "Midpiece","Tail","Total","Head", "Flagellum","Midpiece\n(relative)","Flagellum\n(relative)")))] 
  llvx_[, Velocity:=Motility]

  gEvir =
    ggplot(llvx_, aes(y = trait, x = estimate, shape = Velocity, col = Velocity)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width =width_)) +
    scale_x_continuous(limits = c(-.5, .5), expand = c(0, 0))+
    scale_color_manual(values = col_)  +
    scale_fill_manual(values = col_) + 
 
    labs(y = NULL, x = "Standardized effect size",  tag = '(a)')+
    guides(col=guide_legend(nrow=3,byrow=TRUE,reverse = TRUE),shape = guide_legend(nrow=3,byrow=TRUE,reverse = TRUE))+
    theme_bw() +
    theme(
        plot.margin = margin(18,3,0,1, "mm"),
        plot.subtitle = element_text(size=9, color = 'grey30'),
        plot.tag.position = c(0.036, 1.19),
        plot.tag = element_text(face='bold',size =10),
        #legend.position = "none",
        legend.title = element_text(size=9, color = 'grey30'),
        legend.text=element_text(size=7.5, color = 'grey30'),
        legend.key.height= unit(0.2,"line"),
        legend.margin=margin(0,0,0,0),
        legend.position=c(0.5,1.125),
        #legend.position=c(0.5,1.6),

        axis.title.x = element_text(size = 10, color ='grey10'),
        axis.ticks = element_blank(),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank()
        )    
  #ggsave('Outputs/Fig_4a_width-50mnm_viridis_v2.png',gEvir, width = 4/(5/7), height =10, units = 'cm', bg="white", dpi = 600)     
# (b) - predictions with raw data Fig ER v2 - x-axis labels - illustrations - MOTIL-COL
  size_ =1.2
  line_col = 'red' #grey30
  llvpx_a = llvpx[trait == 'Acrosome']
  llvpx_a[,value:=pred]
  llvpx_n = llvpx[trait == 'Nucleus']
  llvpx_n[,value:=pred]
  llvpx_m = llvpx[trait == 'Midpiece']
  llvpx_m[,value:=pred]
  llvpx_t = llvpx[trait == 'Tail']
  llvpx_t[,value:=pred]
  llvpx_h = llvpx[trait == 'Head']
  llvpx_h[,value:=pred]
  llvpx_f = llvpx[trait == 'Flagellum']
  llvpx_f[,value:=pred]
  llvpx_o = llvpx[trait == 'Total']
  llvpx_o[,value:=pred]
  llvpx_mr = llvpx[trait == 'Midpiece_rel']
  llvpx_mr[,value:=pred]
  llvpx_fr = llvpx[trait == 'Flagellum_rel']
  llvpx_fr[,value:=pred]

  gA =
  ggplot(adwl, aes(x = Acrosome, y = value)) +
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    
    geom_ribbon(data = llvpx_a, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    new_scale_color() +
    geom_line(data = llvpx_a, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+
    
    scale_y_continuous('Velocity [μm/s]', expand = c(0, 0))+
    scale_x_continuous(limits = c(3,5), breaks = c(3,4,5))+
    #coord_cartesian(clip = 'off')+ 
    labs(subtitle = 'Length [μm]')+#, tag = "(b)")+
    #xlab('Acrosome') +
    labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.065, 1.225),
      plot.tag = element_text(face='bold',size =10),

      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      #axis.title.x = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      )  

  gN =
  ggplot(adwl, aes(x = Nucleus, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_n, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_n, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+

    new_scale_color() +
    geom_line(data = llvpx_n, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(26,28,30))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),
      
      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      )    
  
  gM =
  ggplot(adwl, aes(x = Midpiece, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_m, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_m, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(21,24,27))+

    new_scale_color() +
    geom_line(data = llvpx_m, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gT =
  ggplot(adwl, aes(x = Tail, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_t, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_t, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(72, 82, 92))+

    new_scale_color() +
    geom_line(data = llvpx_t, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gH =
  ggplot(adwl, aes(x = Head, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_h, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_h, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(29, 32, 35))+

    new_scale_color() +
    geom_line(data = llvpx_h, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gF =
  ggplot(adwl, aes(x = Flagellum, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_f, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_f, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(95,105, 115))+

    new_scale_color() +
    geom_line(data = llvpx_f, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gO =
  ggplot(adwl, aes(x = Total, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_o, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_o, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(138, 148))+

    new_scale_color() +
    geom_line(data = llvpx_o, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gMR =
  ggplot(adwl, aes(x = Midpiece_rel, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
   stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_mr, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_mr, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(0.15, 0.17, 0.19), labels = c(15,17,19))+

    new_scale_color() +
    geom_line(data = llvpx_mr, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    xlab('Midpiece') +
    labs(subtitle = 'Total length %')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gFR =
  ggplot(adwl, aes(x = Flagellum_rel, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_fr, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_fr, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1, strip.position="right") +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(0.73, 0.76, 0.79), labels = c(73,73,79))+

    new_scale_color() +
    geom_line(data = llvpx_fr, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    xlab('Flagellum') +
    labs(subtitle = '')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),

      strip.text.y.right = element_text(color="grey20",  margin=margin(1,1,1,1,"mm"), angle=90), #size = 7.5
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
  
      panel.border = element_rect(color = 'grey70')
      ) 


  ggR = ggarrange(
    gA, gN, gM, gT, gO, gH, gF, gMR, gFR,
    ncol=9,  widths=c(1.34,1,1,1,1,1,1,1,1.225)
    ) 
  #ggsave('Outputs/Fig_4b.png',ggR, width = 15/(5/7), height = 8, units = 'cm', bg="white", dpi = 600)  
# mix virids the two & export
  blank = ggplot() + theme_void() 
  gB = ggarrange(blank, ggR, nrow=2, heights=c(7.16-5.7,5.7))
  gAll = ggarrange(gEvir, gB, ncol=2, widths=c(4,15))  
  
  # add legend
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

  ymin_ = 0.525 # 0.55, 0.6
  g_anot =   
    gAll + 
    annotation_custom(gp_ind_grob, xmin=.47, xmax=.52, ymin =ymin_)+
    annotation_custom(gi, xmin=.47, xmax=.52, ymin = ymin_+0.2) +
    annotation_custom(gp_sat_grob, xmin=.47+.06, xmax=.52+.06, ymin = ymin_) +
    annotation_custom(gs, xmin=.47+.06, xmax=.52+.06, ymin = ymin_+0.2) +
    annotation_custom(gp_fae_grob, xmin=.47+.12, xmax=.52+.12, ymin = ymin_) +
    annotation_custom(gf, xmin=.47+.12, xmax=.52+.12, ymin = ymin_+0.2) 
  
  ggsave('Outputs/Fig_4_width-190mm_v10.png',g_anot, width = 19/(5/7), height =10, units = 'cm', bg="white", dpi = 600)


#####

# not used (a) & (b) top aligned graph lines
# (a) - effect sizes
  llvx_ = llvx[effect == 'pred']
  llvx_[trait == 'Midpiece_rel', trait:='Midpiece\n(relative)']
  llvx_[trait == 'Flagellum_rel', trait:='Flagellum\n(relative)']
  llvx_[, trait := factor(trait, levels=rev(c("Acrosome", "Nucleus", "Midpiece","Tail","Total","Head", "Flagellum","Midpiece\n(relative)","Flagellum\n(relative)")))] 
  llvx_[, Velocity:=Motility]
  
  gEvir =
    ggplot(llvx_, aes(y = trait, x = estimate, shape = Velocity, col = Velocity)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width =width_)) +
    scale_x_continuous(limits = c(-.5, .5), expand = c(0, 0))+
    scale_color_manual(values = col_)  +
    scale_fill_manual(values = col_) + 
 
    labs(y = NULL, x = "Standardized effect size",  tag = '(a)')+
    guides(col=guide_legend(nrow=3,byrow=TRUE,reverse = TRUE),shape = guide_legend(nrow=3,byrow=TRUE,reverse = TRUE))+
    theme_bw() +
    theme(
        plot.margin = margin(18,3,0,1, "mm"),
        plot.subtitle = element_text(size=9, color = 'grey30'),
        plot.tag.position = c(0.036, 1.03),
        plot.tag = element_text(face='bold',size =10),
        #legend.position = "none",
        legend.title = element_text(size=9, color = 'grey30'),
        legend.text=element_text(size=7.5, color = 'grey30'),
        legend.key.height= unit(0.2,"line"),
        legend.margin=margin(0,0,0,0),
        legend.position=c(0.5,1.125),
        #legend.position=c(0.5,1.6),

        axis.title.x = element_text(size = 10, color ='grey10'),
        axis.ticks = element_blank(),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank()
        )    
  #ggsave('Outputs/Fig_4a_width-50mnm_viridis_v2.png',gEvir, width = 4/(5/7), height =10, units = 'cm', bg="white", dpi = 600)     
# (b) - predictions with raw data Fig ER v2 - x-axis labels - illustrations - MOTIL-COL
  size_ =1.2
  line_col = 'red' #grey30
  llvpx_a = llvpx[trait == 'Acrosome']
  llvpx_a[,value:=pred]
  llvpx_n = llvpx[trait == 'Nucleus']
  llvpx_n[,value:=pred]
  llvpx_m = llvpx[trait == 'Midpiece']
  llvpx_m[,value:=pred]
  llvpx_t = llvpx[trait == 'Tail']
  llvpx_t[,value:=pred]
  llvpx_h = llvpx[trait == 'Head']
  llvpx_h[,value:=pred]
  llvpx_f = llvpx[trait == 'Flagellum']
  llvpx_f[,value:=pred]
  llvpx_o = llvpx[trait == 'Total']
  llvpx_o[,value:=pred]
  llvpx_mr = llvpx[trait == 'Midpiece_rel']
  llvpx_mr[,value:=pred]
  llvpx_fr = llvpx[trait == 'Flagellum_rel']
  llvpx_fr[,value:=pred]

  gA =
  ggplot(adwl, aes(x = Acrosome, y = value)) +
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    
    geom_ribbon(data = llvpx_a, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    new_scale_color() +
    geom_line(data = llvpx_a, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+
    
    scale_y_continuous('Velocity [μm/s]', expand = c(0, 0))+
    scale_x_continuous(limits = c(3,5), breaks = c(3,4,5))+
    #coord_cartesian(clip = 'off')+ 
    labs(subtitle = 'Length [μm]')+#, tag = "(b)")+
    #xlab('Acrosome') +
    labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.065, 0.97),
      plot.tag = element_text(face='bold',size =10),

      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      #axis.title.x = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      )  

  gN =
  ggplot(adwl, aes(x = Nucleus, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_n, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_n, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+

    new_scale_color() +
    geom_line(data = llvpx_n, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(26,28,30))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),
      
      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      )    
  
  gM =
  ggplot(adwl, aes(x = Midpiece, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_m, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_m, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(21,24,27))+

    new_scale_color() +
    geom_line(data = llvpx_m, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gT =
  ggplot(adwl, aes(x = Tail, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_t, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_t, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(72, 82, 92))+

    new_scale_color() +
    geom_line(data = llvpx_t, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gH =
  ggplot(adwl, aes(x = Head, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_h, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_h, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(29, 32, 35))+

    new_scale_color() +
    geom_line(data = llvpx_h, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gF =
  ggplot(adwl, aes(x = Flagellum, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_f, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_f, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(95,105, 115))+

    new_scale_color() +
    geom_line(data = llvpx_f, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gO =
  ggplot(adwl, aes(x = Total, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_o, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_o, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(138, 148))+

    new_scale_color() +
    geom_line(data = llvpx_o, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gMR =
  ggplot(adwl, aes(x = Midpiece_rel, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
   stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_mr, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_mr, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(0.15, 0.17, 0.19), labels = c(15,17,19))+

    new_scale_color() +
    geom_line(data = llvpx_mr, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    xlab('Midpiece') +
    labs(subtitle = 'Total length %')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gFR =
  ggplot(adwl, aes(x = Flagellum_rel, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_fr, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_fr, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1, strip.position="right") +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(0.73, 0.76, 0.79), labels = c(73,73,79))+

    new_scale_color() +
    geom_line(data = llvpx_fr, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    xlab('Flagellum') +
    labs(subtitle = '')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),

      strip.text.y.right = element_text(color="grey20",  margin=margin(1,1,1,1,"mm"), angle=90), #size = 7.5
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
  
      panel.border = element_rect(color = 'grey70')
      ) 


  ggR = ggarrange(
    gA, gN, gM, gT, gO, gH, gF, gMR, gFR,
    ncol=9,  widths=c(1.34,1,1,1,1,1,1,1,1.225)
    ) 
  #ggsave('Outputs/Fig_4b.png',ggR, width = 15/(5/7), height = 8, units = 'cm', bg="white", dpi = 600)  
# mix virids the two & export
  blank = ggplot() + theme_void() 
  gB = ggarrange(blank, ggR, blank, nrow=3, heights=c(0.93,5.7,0.53))
  gAll = ggarrange(gEvir, gB, ncol=2, widths=c(4,15))  

  # add legend
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

  ymin_ = 0.65 # 0.55, 0.6
  g_anot =   
    gAll + 
    annotation_custom(gp_ind_grob, xmin=.47, xmax=.52, ymin =ymin_)+
    annotation_custom(gi, xmin=.47, xmax=.52, ymin = ymin_+0.2) +
    annotation_custom(gp_sat_grob, xmin=.47+.06, xmax=.52+.06, ymin = ymin_) +
    annotation_custom(gs, xmin=.47+.06, xmax=.52+.06, ymin = ymin_+0.2) +
    annotation_custom(gp_fae_grob, xmin=.47+.12, xmax=.52+.12, ymin = ymin_) +
    annotation_custom(gf, xmin=.47+.12, xmax=.52+.12, ymin = ymin_+0.2) 
  
  ggsave('Outputs/Fig_4_width-190mm_v11.png',g_anot, width = 19/(5/7), height =10, units = 'cm', bg="white", dpi = 600)


# not used alternative (a) & (b) positions
# (a) - effect sizes
  llvx_ = llvx[effect == 'pred']
  llvx_[trait == 'Midpiece_rel', trait:='Midpiece\n(relative)']
  llvx_[trait == 'Flagellum_rel', trait:='Flagellum\n(relative)']
  llvx_[, trait := factor(trait, levels=rev(c("Acrosome", "Nucleus", "Midpiece","Tail","Total","Head", "Flagellum","Midpiece\n(relative)","Flagellum\n(relative)")))] 

  gEvir =
    ggplot(llvx_, aes(y = trait, x = estimate, shape = Motility, col = Motility)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width =width_)) +
    scale_x_continuous(limits = c(-.5, .5), expand = c(0, 0))+
    scale_color_manual(values = col_)  +
    scale_fill_manual(values = col_) + 
 
    labs(y = NULL, x = "Standardized effect size",  tag = '(a)')+
    guides(col=guide_legend(nrow=3,byrow=TRUE,reverse = TRUE),shape = guide_legend(nrow=3,byrow=TRUE,reverse = TRUE))+
    theme_bw() +
    theme(
        plot.margin = margin(18,3,0,1, "mm"),
        plot.subtitle = element_text(size=9, color = 'grey30'),
        plot.tag.position = c(0.035, 1.19),
        plot.tag = element_text(face='bold',size =10),
        #legend.position = "none",
        legend.title = element_text(size=9, color = 'grey30'),
        legend.text=element_text(size=7.5, color = 'grey30'),
        legend.key.height= unit(0.2,"line"),
        legend.margin=margin(0,0,0,0),
        legend.position=c(0.5,1.125),
        #legend.position=c(0.5,1.6),

        axis.title.x = element_text(size = 10, color ='grey10'),
        axis.ticks = element_blank(),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank()
        )    
  #ggsave('Outputs/Fig_4a_width-50mnm_viridis_v2.png',gEvir, width = 4/(5/7), height =10, units = 'cm', bg="white", dpi = 600)     
# (b) - predictions with raw data Fig ER v2 - x-axis labels - illustrations - MOTIL-COL
  size_ =1.2
  line_col = 'red' #grey30
  llvpx_a = llvpx[trait == 'Acrosome']
  llvpx_a[,value:=pred]
  llvpx_n = llvpx[trait == 'Nucleus']
  llvpx_n[,value:=pred]
  llvpx_m = llvpx[trait == 'Midpiece']
  llvpx_m[,value:=pred]
  llvpx_t = llvpx[trait == 'Tail']
  llvpx_t[,value:=pred]
  llvpx_h = llvpx[trait == 'Head']
  llvpx_h[,value:=pred]
  llvpx_f = llvpx[trait == 'Flagellum']
  llvpx_f[,value:=pred]
  llvpx_o = llvpx[trait == 'Total']
  llvpx_o[,value:=pred]
  llvpx_mr = llvpx[trait == 'Midpiece_rel']
  llvpx_mr[,value:=pred]
  llvpx_fr = llvpx[trait == 'Flagellum_rel']
  llvpx_fr[,value:=pred]

  gA =
  ggplot(adwl, aes(x = Acrosome, y = value)) +
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    
    geom_ribbon(data = llvpx_a, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    new_scale_color() +
    geom_line(data = llvpx_a, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+
    
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(limits = c(3,5), breaks = c(3,4,5))+
    #coord_cartesian(clip = 'off')+ 
    labs(subtitle = 'Length [μm]')+#, tag = "(b)")+
    #xlab('Acrosome') +
    labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.07, 0.97),
      plot.tag = element_text(face='bold',size =10),

      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      #axis.title.x = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      )  

  gN =
  ggplot(adwl, aes(x = Nucleus, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_n, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_n, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+

    new_scale_color() +
    geom_line(data = llvpx_n, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(26,28,30))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),
      
      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      )    
  
  gM =
  ggplot(adwl, aes(x = Midpiece, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_m, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_m, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(21,24,27))+

    new_scale_color() +
    geom_line(data = llvpx_m, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gT =
  ggplot(adwl, aes(x = Tail, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_t, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_t, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(72, 82, 92))+

    new_scale_color() +
    geom_line(data = llvpx_t, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gH =
  ggplot(adwl, aes(x = Head, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_h, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_h, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(29, 32, 35))+

    new_scale_color() +
    geom_line(data = llvpx_h, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gF =
  ggplot(adwl, aes(x = Flagellum, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_f, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_f, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(95,105, 115))+

    new_scale_color() +
    geom_line(data = llvpx_f, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gO =
  ggplot(adwl, aes(x = Total, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_o, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_o, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(138, 148))+

    new_scale_color() +
    geom_line(data = llvpx_o, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gMR =
  ggplot(adwl, aes(x = Midpiece_rel, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
   stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_mr, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_mr, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(0.15, 0.17, 0.19), labels = c(15,17,19))+

    new_scale_color() +
    geom_line(data = llvpx_mr, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    xlab('Midpiece') +
    labs(subtitle = 'Total length %')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gFR =
  ggplot(adwl, aes(x = Flagellum_rel, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_fr, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    #geom_line(data = llvpx_fr, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1, strip.position="right") +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(0.73, 0.76, 0.79), labels = c(73,73,79))+

    new_scale_color() +
    geom_line(data = llvpx_fr, aes(x = Length_avg, y =pred, col=Motility))+
    scale_color_manual(NULL, values = col_2)+

    xlab('Flagellum') +
    labs(subtitle = '')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),

      strip.text.y.right = element_text(color="grey20",  margin=margin(1,1,1,1,"mm"), angle=90), #size = 7.5
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
  
      panel.border = element_rect(color = 'grey70')
      ) 


  ggR = ggarrange(
    gA, gN, gM, gT, gO, gH, gF, gMR, gFR,
    ncol=9,  widths=c(1.34,1,1,1,1,1,1,1,1.225)
    ) 
  ggsave('Outputs/Fig_M_v4.png',ggR, width = 15/(5/7), height = 8, units = 'cm', bg="white", dpi = 600)  
# mix virids the two & export
  blank = ggplot() + theme_void() 
  gB = ggarrange(blank, ggR, nrow=2, heights=c(7.16-5.7,5.7))
  gAll = ggarrange(gEvir, gB, ncol=2, widths=c(4,15))  
  
  # add legend
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

  ymin_ = 0.525 # 0.55, 0.6
  g_anot =   
    gAll + 
    annotation_custom(gp_ind_grob, xmin=.47, xmax=.52, ymin =ymin_)+
    annotation_custom(gi, xmin=.47, xmax=.52, ymin = ymin_+0.2) +
    annotation_custom(gp_sat_grob, xmin=.47+.06, xmax=.52+.06, ymin = ymin_) +
    annotation_custom(gs, xmin=.47+.06, xmax=.52+.06, ymin = ymin_+0.2) +
    annotation_custom(gp_fae_grob, xmin=.47+.12, xmax=.52+.12, ymin = ymin_) +
    annotation_custom(gf, xmin=.47+.12, xmax=.52+.12, ymin = ymin_+0.2) 
  ggsave('Outputs/Fig_4_width-190mm_illust_viridis_virids.png',g_anot, width = 19/(5/7), height =10, units = 'cm', bg="white", dpi = 600)

# not used (b) - predictions with raw data Fig ER v2 - x-axis labels - illustrations - RED
  size_ =1.2
  line_col = 'red' #grey30
  llvpx_a = llvpx[trait == 'Acrosome']
  llvpx_a[,value:=pred]
  llvpx_n = llvpx[trait == 'Nucleus']
  llvpx_n[,value:=pred]
  llvpx_m = llvpx[trait == 'Midpiece']
  llvpx_m[,value:=pred]
  llvpx_t = llvpx[trait == 'Tail']
  llvpx_t[,value:=pred]
  llvpx_h = llvpx[trait == 'Head']
  llvpx_h[,value:=pred]
  llvpx_f = llvpx[trait == 'Flagellum']
  llvpx_f[,value:=pred]
  llvpx_o = llvpx[trait == 'Total']
  llvpx_o[,value:=pred]
  llvpx_mr = llvpx[trait == 'Midpiece_rel']
  llvpx_mr[,value:=pred]
  llvpx_fr = llvpx[trait == 'Flagellum_rel']
  llvpx_fr[,value:=pred]

  gA =
  ggplot(adwl, aes(x = Acrosome, y = value)) +
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    
    geom_ribbon(data = llvpx_a, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_a, aes(x = Length_avg, y =pred), col = 'red')+
    
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(limits = c(3,5), breaks = c(3,4,5))+
    #coord_cartesian(clip = 'off')+ 
    labs(subtitle = 'Length [μm]')+#, tag = "(b)")+
    #xlab('Acrosome') +
    labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.07, 0.97),
      plot.tag = element_text(face='bold',size =10),

      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      #axis.title.x = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      )  

  gN =
  ggplot(adwl, aes(x = Nucleus, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_n, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_n, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(26,28,30))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),
      
      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      )    
  
  gM =
  ggplot(adwl, aes(x = Midpiece, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_m, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_m, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(21,24,27))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gT =
  ggplot(adwl, aes(x = Tail, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_t, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_t, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(72, 82, 92))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gH =
  ggplot(adwl, aes(x = Head, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_h, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_h, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(29, 32, 35))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gF =
  ggplot(adwl, aes(x = Flagellum, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_f, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_f, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(95,105, 115))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gO =
  ggplot(adwl, aes(x = Total, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_o, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_o, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(138, 148))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gMR =
  ggplot(adwl, aes(x = Midpiece_rel, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
   stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_mr, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_mr, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(0.15, 0.17, 0.19), labels = c(15,17,19))+
    xlab('Midpiece') +
    labs(subtitle = 'Total length %')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  
  gFR =
  ggplot(adwl, aes(x = Flagellum_rel, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_fr, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_fr, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1, strip.position="right") +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    scale_x_continuous(breaks = c(0.73, 0.76, 0.79), labels = c(73,73,79))+
    xlab('Flagellum') +
    labs(subtitle = '')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30', vjust = -1),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),

      strip.text.y.right = element_text(color="grey20",  margin=margin(1,1,1,1,"mm"), angle=90), #size = 7.5
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
  
      panel.border = element_rect(color = 'grey70')
      ) 


  ggR = ggarrange(
    gA, gN, gM, gT, gO, gH, gF, gMR, gFR,
    ncol=9,  widths=c(1.34,1,1,1,1,1,1,1,1.225)
    ) 
  ggsave('Outputs/Fig_M_v3.png',ggR, width = 15/(5/7), height = 8, units = 'cm', bg="white", dpi = 600)
# TEMP test like in Knief et al
    i = 'Nucleus'
    xi = adwl[, i, with = FALSE] 
    adwl[, Length_avg := xi] 
    k = 'Average path'
    adwlk = adwl[Motility == k]
    m = lm(scale(value) ~ scale(log(motileCount)) + Morph + scale(Length_avg), adwlk)
    adwlk[,predictions:=predict(m)]

    ggplot(adwlk, aes(y=value, x = predictions, col = Morph)) + geom_point(alpha=0.5)
    ggplot(adwlk, aes(y=value, x = predictions, col = Morph)) + geom_point(pch = 1)

# END