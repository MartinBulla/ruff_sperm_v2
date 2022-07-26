# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  require(ggpubr) 
  require(ggsci)
  require(MuMIn)

  width_ = 0.8 # spacing between error bars

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
      adwl[, Motility := factor(Motility, levels=rev(c("Curvilinear", "Straight line", "Average path")))] 

  adw[,motileCount_ln_z := scale(log(motileCount))]
  adwl[,motileCount_ln_z := scale(log(motileCount))]
  ar[,motileCount_ln_z := scale(log(motileCount))]
 
# compare AIC
   l = list()
   for(i in unique(adwl$mot)){
        #i ='Straight line'
        ai = adwl[mot == i]
        
        ma = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Acrosome), ai)
        map = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Acrosome,2), ai)
        aic = data.table(AICc(ma, map))
        aic[, part:='Acrosome']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mn = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Nucleus), ai)
        mnp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Nucleus,2), ai)
        aic = data.table(AICc(mn, mnp))
        aic[, part:='Nucleus']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mm = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Midpiece), ai)
        mmp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Midpiece,2), ai)
        aic = data.table(AICc(mm, mmp))
        aic[, part:='Midpiece']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mt = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Tail), ai)
        mtp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Tail,2), ai)
        aic = data.table(AICc(mt, mtp))
        aic[, part:='Tail']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mo = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Total), ai)
        mop = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Total,2), ai)
        aic = data.table(AICc(mo, mop))
        aic[, part:='Total']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic
        
        mh = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Head), ai)
        mhp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Head,2), ai)
        aic = data.table(AICc(mh, mhp))
        aic[, part:='Head']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic
        

        mf = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Flagellum), ai)
        mfp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Flagellum,2), ai)
        aic = data.table(AICc(mf, mfp))
        aic[, part:='Flagellum']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mmr = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Midpiece_rel), ai)
        mmrp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Midpiece_rel,2), ai)
        aic = data.table(AICc(mmr, mmrp))
        aic[, part:='Midpiece_rel']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mfr = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Flagellum_rel), ai)
        mfrp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Flagellum_rel,2), ai)
        aic = data.table(AICc(mfr, mfrp))
        aic[, part:='Flagellum_rel']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic
       
        }    
   aic = do.call(rbind, l)
   #aic[df == 6, poly:='no']
   #aic[df == 7, poly:='yes']
   aic[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Midpiece","Tail","Total","Head","Flagellum","Midpiece_rel","Flagellum_rel")))] 
   aic[, mot := factor(mot, levels=rev(c("Curvilinear", "Straight line", "Average path")))] 
   aic = merge(aic[df==6], aic[df==7], by = c('part', 'mot'))
   setnames(aic, old = c('AICc.x', 'AICc.y'), new = c('simple', 'quadratic'))
   aic$df.x = aic$df.y = NULL
   aic = aic[,.(mot, part, simple, quadratic)]
   aic[,simple:=round(simple,2)]
   aic[,quadratic:=round(quadratic,2)]
   aic[, deltaAICc:=round(quadratic-simple,2)]
   aic[, prob := round(exp(-0.5*deltaAICc)/sum(exp(-0.5*deltaAICc)),2)]
   aic[, ER := round((1-prob)/prob, 2)]    
   aic = aic[order(mot,part, decreasing = TRUE)]   
   aic_w = reshape(aic, idvar = c('part'), timevar = 'mot', direction = "wide")  
   fwrite(aic, file = 'Outputs/Table_SpolyAIC_l.csv') 
   fwrite(aic_w, file = 'Outputs/Table_SpolyAIC_w.csv')
# prepare estimates and pred for univariate models
  effects_ = c('intercept','motileCount_ln', 'morphSat', 'morphFae', 'pred')
  lvx = list()
  lvpx =list()
  for(i in c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')){
    #i = 'Acrosome'
    xi = adwl[, i, with = FALSE] 
    adwl[, Length_avg := xi] 
    for(k in unique(adwl$Motility)){
      #k = 'VCL'
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
  llvx[response == 'VAP' ,Motility:='Average path']
  llvx[response == 'VCL' ,Motility:='Curvilinear']
  llvx[response == 'VSL' ,Motility:='Straight line']
  llvx[, Motility := factor(Motility, levels=rev(c("Curvilinear", "Straight line", "Average path")))] 
  llvx[, trait := factor(trait, levels=rev(c(c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel'))))] 

  llvpx = data.table(do.call(rbind,lvpx) ) 
  llvpx[response == 'VAP' ,Motility:='Average path']
  llvpx[response == 'VCL' ,Motility:='Curvilinear']
  llvpx[response == 'VSL' ,Motility:='Straight line']
  llvpx[, Motility := factor(Motility, levels=(c("Curvilinear", "Straight line", "Average path")))] 
  llvpx[, trait := factor(trait, levels=rev(c(c('Acrosome', 'Nucleus', 'Midpiece', 'Tail', 'Total', 'Head','Flagellum', 'Midpiece_rel', 'Flagellum_rel'))))] 

# plot effect sizes
  llvx_ = llvx[effect == 'pred']
  gE = 
    ggplot(llvx_, aes(y = trait, x = estimate, shape = Motility, col = Motility)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width =width_)) +
    scale_x_continuous(limits = c(-.5, .5), expand = c(0, 0))+
    scale_color_jco()+
    labs(y = NULL, x = "Standardized effect size")+
    guides(col=guide_legend(nrow=3,byrow=TRUE,reverse = TRUE),shape = guide_legend(nrow=3,byrow=TRUE,reverse = TRUE))+
    theme_bw() +
    theme(
        plot.subtitle = element_text(size=9, color = 'grey30'),
        #legend.position = "none",
        legend.title = element_text(size=9, color = 'grey30'),
        legend.text=element_text(size=7.5, color = 'grey30'),
        legend.key.height= unit(0.2,"line"),
        legend.margin=margin(0,0,0,0),
        legend.position=c(0.5,1.1),
        #legend.position=c(0.5,1.6),

        axis.title.x = element_text(size = 10, color ='grey10'),
        axis.ticks = element_blank(),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank(),

        plot.margin = margin(18,3,1,1, "mm")
        )     
# START HERE AND FIX THE PLOT MARGINS TO MAKE LITTLE GAB BETWEEN PLOTSplot predictions with raw data Fig ER v2 - x-axis labels - illustrations
  size_ =1.2
  
  llvpx_a = llvpx[trait == 'Acrosome']
  llvpx_a[,value:=pred]
  gA =
  ggplot(adwl, aes(x = Acrosome, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    geom_ribbon(data = llvpx_a, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_a, aes(x = Length_avg, y =pred), col ='grey30')+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    labs(subtitle = 'Length [μm]')+#, tag = "(b)")+
    #xlab('Acrosome') +
    labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),

      plot.subtitle = element_text(size=9, color = 'grey30'),

      axis.title = element_text(size = 10, , colour="grey10"),
      #axis.title.x = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      )  
  llvpx_n = llvpx[trait == 'Nucleus']
  llvpx_n[,value:=pred]
  gN =
  ggplot(adwl, aes(x = Nucleus, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    geom_ribbon(data = llvpx_n, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_n, aes(x = Length_avg, y =pred), col ='grey30')+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30'),
      
      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      )    
  llvpx_m = llvpx[trait == 'Midpiece']
  llvpx_m[,value:=pred]
  gM =
  ggplot(adwl, aes(x = Midpiece, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    geom_ribbon(data = llvpx_m, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_m, aes(x = Length_avg, y =pred), col ='grey30')+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30'),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  llvpx_t = llvpx[trait == 'Tail']
  llvpx_t[,value:=pred]
  gT =
  ggplot(adwl, aes(x = Tail, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    geom_ribbon(data = llvpx_t, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_t, aes(x = Length_avg, y =pred), col ='grey30')+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30'),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  llvpx_h = llvpx[trait == 'Head']
  llvpx_h[,value:=pred]
  gH =
  ggplot(adwl, aes(x = Head, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    geom_ribbon(data = llvpx_h, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_h, aes(x = Length_avg, y =pred), col ='grey30')+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30'),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  llvpx_f = llvpx[trait == 'Flagellum']
  llvpx_f[,value:=pred]
  gF =
  ggplot(adwl, aes(x = Flagellum, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    geom_ribbon(data = llvpx_f, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_f, aes(x = Length_avg, y =pred), col ='grey30')+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30'),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  llvpx_o = llvpx[trait == 'Total']
  llvpx_o[,value:=pred]
  gO =
  ggplot(adwl, aes(x = Total, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    geom_ribbon(data = llvpx_o, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_o, aes(x = Length_avg, y =pred), col ='grey30')+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    labs(subtitle = '')+
    #labs(tag = '(b)')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30'),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  llvpx_mr = llvpx[trait == 'Midpiece_rel']
  llvpx_mr[,value:=pred]
  gMR =
  ggplot(adwl, aes(x = Midpiece_rel, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    geom_ribbon(data = llvpx_mr, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_mr, aes(x = Length_avg, y =pred), col ='grey30')+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    xlab('Midpiece') +
    labs(subtitle = 'Relative to total length')+
    theme_bw() +
    theme(
      legend.position = "none",
      #plot.margin = margin(c(l = 0, r = 0), unit = "mm"),
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30'),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
  
      panel.border = element_rect(color = 'grey70')
      ) 
  llvpx_fr = llvpx[trait == 'Flagellum_rel']
  llvpx_fr[,value:=pred]
  gFR =
  ggplot(adwl, aes(x = Flagellum_rel, y = value)) +
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    geom_ribbon(data = llvpx_fr, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_fr, aes(x = Length_avg, y =pred), col ='grey30')+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1, strip.position="right") +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Motility [μm/s]', expand = c(0, 0))+
    xlab('Flagellum') +
    labs(subtitle = '')+
    theme_bw() +
    theme(
      legend.position = "none",
      plot.tag.position = c(0.005, 1),
      plot.tag = element_text(face='bold',size =10),
      plot.subtitle = element_text(size=9, color = 'grey30'),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.y = element_blank(), 
      #axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),

      strip.text.y.right = element_text(color="grey20",  margin=margin(1,1,1,1,"mm"), angle=90), #size = 7.5
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
  
      panel.border = element_rect(color = 'grey70')
      ) 

    grid.draw(cbind(
    ggplotGrob(gA), ggplotGrob(gN), ggplotGrob(gM),ggplotGrob(gT), ggplotGrob(gO), 
    ggplotGrob(gH),ggplotGrob(gF), ggplotGrob(gMR), ggplotGrob(gFR),
    size = "first")
    )

   ggR = ggarrange(
    gA, gN, gM, gT, gO, gH, gF, gMR, gFR,
    ncol=9,  widths=c(1.195,1,1,1,1,1,1,1,1.13)
    ) 

   ggsave('Outputs/Fig_M_v1.png',ggR, width = 15/(5/7), height = 8, units = 'cm', bg="white", dpi = 600)
   
  right = 0.38
  right2 = 0.01
  ymin_ = -0.92
    
  ggExp =  
  ggR + 
    annotation_custom(gi, xmin=0.064, xmax=0.116, ymin=ymin_) + 
    annotation_custom(gs, xmin=0.064+0.05, xmax=0.116+0.05, ymin=ymin_) + 
    annotation_custom(gf, xmin=0.064+0.1, xmax=0.116+0.1, ymin=ymin_)+

    annotation_custom(gi, xmin=0.064+right, xmax=0.116+right, ymin=ymin_) + 
    annotation_custom(gs, xmin=0.064+0.05+right, xmax=0.116+0.05+right, ymin=ymin_) + 
    annotation_custom(gf, xmin=0.064+0.1+right, xmax=0.116+0.1+right, ymin=ymin_) +

    annotation_custom(gi, xmin=0.064+2*right+right2, xmax=0.116+2*right+right2, ymin=ymin_) + 
    annotation_custom(gs, xmin=0.064+0.05+2*right+right2, xmax=0.116+0.05+2*right+right2, ymin=ymin_) + 
    annotation_custom(gf, xmin=0.064+0.1+2*right+right2, xmax=0.116+0.1+2*right+right2, ymin=ymin_) 
  #ggExp   
  ggsave('Outputs/Fig_ER_130mm_v2.png',ggExp, width = 13/(5/7), height =13, units = 'cm', bg="white", dpi = 600)




# START HERE

  setnames(bw,old = c('Length_µm.Acrosome', 'Length_µm.Nucleus','Length_µm.Head','Length_µm.Midpiece','Length_µm.Tail','Length_µm.Flagellum', 'Length_µm.Total'), new = c('Acrosome', 'Nucleus', 'Head','Midpiece', 'Tail','Flagellum','Total'))


   cbind(aic[mot == 'Curvilinear'], aic[mot == 'Straight line',.('part', 'simple', 'quadratic', 'deltaAICc', 'prob', 'ER')],  aic[mot == 'Average path',.('part','simple', 'quadratic', 'deltaAICc', 'prob', 'ER')], by = c('part'))
    
    aic[mot == 'Average path'],

    m = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Head) + scale(Midpiece)+scale(Tail), adw)
    mp = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Head,2) + poly(Midpiece,2)+poly(Tail,2), adw)
    mi = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Head)*scale(Midpiece)+
            scale(Head)*scale(Tail) +
            scale(Midpiece)*scale(Tail), adw)
    mpi = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Head,2)*poly(Midpiece,2)+
            poly(Head,2)*poly(Tail,2) +
            poly(Midpiece,2)*poly(Tail,2), adw)
       


  }

   lz = list()
  lpz =list()
  lprz =list()

  # VAP
     for(i in unique(a$part)){
        #i ='Tail'
        ai = a[part == i]
        m = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Length_avg), ai)
        mp = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(scale(Length_avg),2), ai)
        mp = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Length_avg,2), ai)
        

        ma = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Acrosome), adw)
        map = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Acrosome,2), adw)

        mn = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Nucleus), adw)
        mnp = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Nucleus,2), adw)
        
        mo = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Total), adw)
        mop = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Total,2), adw)

        mf = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Flagellum_rel), adw)
        mfp = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Flagellum_rel,2), adw)

        mh = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Head), adw)
        mhp = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Head,2), adw)
        mm = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Midpiece), adw)
        mmp = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Midpiece,2), adw)
        mt = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Tail), adw)
        mtp = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Tail,2), adw)
        m = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Head) + scale(Midpiece)+scale(Tail), adw)
        mp = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Head,2) + poly(Midpiece,2)+poly(Tail,2), adw)
        mi = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Head)*scale(Midpiece)+
                scale(Head)*scale(Tail) +
                scale(Midpiece)*scale(Tail), adw)
        mpi = lm(scale(VAP) ~ motileCount_ln_z+ Morph+poly(Head,2)*poly(Midpiece,2)+
                poly(Head,2)*poly(Tail,2) +
                poly(Midpiece,2)*poly(Tail,2), adw)
       
        aic = data.table(AICc(mo, mop, mf, mfp, mh, mhp, mm, mmp, mt, mtp, m, mp, mi,mpi))
        aic[,model := c('Total', 'Total_poly', 'Flagellum_rel', 'Flagellum_rel_poly','Head', 'Head_poly', 'Midpiece','Midpiece_poly','Tail', 'Tail_poly', 'Head_Mid_Tail', 'Head_Mid_Tail_poly','Head_Mid_Tail_2way', 'Head_Mid_Tail_poly_2way' )]
        aic = aic[order(AICc)]
        aic[, deltaAICc:=AICc-min(AICc)]
        aic[, prob := round(exp(-0.5*deltaAICc)/sum(exp(-0.5*deltaAICc)),2)]
        aic[, ER := round(max(prob)/prob, 2)]
        aic

        aic = data.table(AICc(ma, mn, mo, mf, mh, mm, mt, m, mi, m2))
        aic[,model := c('Acrosome','Nucleus','Total',  'Flagellum_rel', 'Head', 'Midpiece','Tail',  'Head_Mid_Tail','Head_Mid_Tail_2way','Head_Mid_Tail_Tot' )]
        aic = aic[order(AICc)]
        aic[, deltaAICc:=AICc-min(AICc)]
        aic[, prob := round(exp(-0.5*deltaAICc)/sum(exp(-0.5*deltaAICc)),2)]
        aic[, ER := round(max(prob)/prob, 2)]
        aic
        require(car)
        vif(mi)
        vif(mi, terms = 'marginal')
        vif(mi, terms = 'high-order')



        plot(allEffects(mi))
        plot(allEffects(mo))
        #summary(mp)
        #plot(allEffects(mp))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lz[[paste('VAP',i)]]=data.frame(mot = 'VAP', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VAP ~ motileCount_ln_z + Morph+Length_avg, ai)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        nd = list()
        for(j in unique(ai$Morph)){
          #j=1
          nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
          }
        newD=do.call(rbind,nd)

        # values to predict for
        X <- model.matrix(~ motileCount_ln_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
        newD$VAP <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                        predmatrix[predmatrix < 0] <- 0
                        newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                        newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                        #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$part=i
        newD$mot = 'VAP'
        setnames(newD, old = 'VAP', new = 'motility')
        lpz[[paste(i,newD$mot[1])]] = data.table(newD)

        print(paste(i,newD$mot[1]))     
        }          
     for(i in unique(ar$part)){
        #i ='Nucleus'
        if(i == 'Midpiece'){ii = 'Midpiece_rel'}
        if(i == 'Flagellum'){ii = 'Flagellum_rel'}
        ai = ar[part == i]
        m = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Length_rel), ai)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lz[[paste('VAP',ii)]]=data.frame(mot = 'VAP', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VAP ~ motileCount_ln_z + Morph+Length_rel, ai)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        nd = list()
        for(j in unique(ai$Morph)){
          #j=1
          nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
          }
        newD=do.call(rbind,nd)

        # values to predict for
        X <- model.matrix(~ motileCount_ln_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
        newD$VAP <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                        predmatrix[predmatrix < 0] <- 0
                        newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                        newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                        #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$part=ii
        newD$mot = 'VAP'
        setnames(newD, old = 'VAP', new = 'motility')
        lpz[[paste(ii,newD$mot[1])]] = data.table(newD)

        print(paste(ii,newD$mot[1]))     
        }            
  # VSL
     for(i in unique(a$part)){
        #i ='Nucleus'
        ai = a[part == i]
        m = lm(scale(VSL) ~ motileCount_ln_z + Morph+scale(Length_avg), ai)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lz[[paste('VSL',i)]]=data.frame(mot = 'VSL', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VSL ~ motileCount_ln_z + Morph+Length_avg, ai)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        nd = list()
        for(j in unique(ai$Morph)){
          #j=1
          nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
          }
        newD=do.call(rbind,nd)
        X <- model.matrix(~ motileCount_ln_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
        newD$VSL <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                        predmatrix[predmatrix < 0] <- 0
                        newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                        newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                        #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$part=i
        newD$mot = 'VSL'
        setnames(newD, old = 'VSL', new = 'motility')
        lpz[[paste(i,newD$mot[1])]] = data.table(newD)

        print(paste(i,newD$mot[1]))     
        }          
     for(i in unique(ar$part)){
        #i ='Nucleus'
        if(i == 'Midpiece'){ii = 'Midpiece_rel'}
        if(i == 'Flagellum'){ii = 'Flagellum_rel'}
        ai = ar[part == i]
        m = lm(scale(VSL) ~ motileCount_ln_z+ Morph+scale(Length_rel), ai)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lz[[paste('VSL',ii)]]=data.frame(mot = 'VSL', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VSL ~ motileCount_ln_z + Morph+Length_rel, ai)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        nd = list()
        for(j in unique(ai$Morph)){
          #j=1
          nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
          }
        newD=do.call(rbind,nd)

        # values to predict for
        X <- model.matrix(~ motileCount_ln_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
        newD$VSL <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                        predmatrix[predmatrix < 0] <- 0
                        newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                        newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                        #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$part=ii
        newD$mot = 'VSL'
        setnames(newD, old = 'VSL', new = 'motility')
        lpz[[paste(ii,newD$mot[1])]] = data.table(newD)

        print(paste(ii,newD$mot[1]))     
        }            
  # VCL
     for(i in unique(a$part)){
        #i ='Nucleus'
        ai = a[part == i]
        m = lm(scale(VCL) ~ motileCount_ln_z + Morph+scale(Length_avg), ai)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lz[[paste('VCL',i)]]=data.frame(mot = 'VCL', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VCL ~ motileCount_ln_z + Morph+Length_avg, ai)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        nd = list()
        for(j in unique(ai$Morph)){
          #j=1
          nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
          }
        newD=do.call(rbind,nd)
        X <- model.matrix(~ motileCount_ln_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
        newD$VCL <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                        predmatrix[predmatrix < 0] <- 0
                        newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                        newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                        #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$part=i
        newD$mot = 'VCL'
        setnames(newD, old = 'VCL', new = 'motility')
        lpz[[paste(i,newD$mot[1])]] = data.table(newD)

        print(paste(i,newD$mot[1]))     
        }          
     for(i in unique(ar$part)){
        #i ='Nucleus'
        if(i == 'Midpiece'){ii = 'Midpiece_rel'}
        if(i == 'Flagellum'){ii = 'Flagellum_rel'}
        ai = ar[part == i]
        m = lm(scale(VCL) ~ motileCount_ln_z+ Morph+scale(Length_rel), ai)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lz[[paste('VCL',ii)]]=data.frame(mot = 'VCL', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VCL ~ motileCount_ln_z + Morph+Length_rel, ai)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        nd = list()
        for(j in unique(ai$Morph)){
          #j=1
          nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
          }
        newD=do.call(rbind,nd)

        # values to predict for
        X <- model.matrix(~ motileCount_ln_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
        newD$VCL <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                        predmatrix[predmatrix < 0] <- 0
                        newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                        newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                        #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$part=ii
        newD$mot = 'VCL'
        setnames(newD, old = 'VCL', new = 'motility')
        lpz[[paste(ii,newD$mot[1])]] = data.table(newD)

        print(paste(ii,newD$mot[1]))     
        }          
           
  llz = data.table(do.call(rbind,lz) ) 
  llz[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
  llz[effect == '(Intercept)', effect:='Independent\n(Intercept)']
  llz[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
  llz[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
  llz[effect == 'scale(Length_avg)', effect := 'Length_µm']
  llz[effect == 'scale(Length_rel)', effect := 'Length_µm'] # dummy variable
  
  llz[, effect := factor(effect, levels=c('Length_µm','motileCount_ln_z',"Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

  llz[mot=='VCL', mot:='Curvilinear (VCL)']
  llz[mot=='VSL', mot:='Straight-line (VSL)']
  llz[mot=='VAP', mot:='Average-path (VAP)']
  llz[, mot := factor(mot, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
