# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  require(ggpubr) 
  require(ggsci)
  require(MuMIn)

  width_ = 0.6 # spacing between error bars

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
 
# compare AIC for simple and quadratic
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

# compare AIC for simple and interactions
  l = list()
  for(i in unique(adwl$mot)){
    #i ='Curvilinear'
    ai = adwl[mot == i]
    ms = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Head)+scale(Midpiece)+scale(Tail), ai)
    mi = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Head)*scale(Midpiece)+scale(Head)*scale(Tail) + scale(Midpiece)*scale(Tail), ai)
    aic = data.table(AICc(ms,mi))
    aic[, mot:=i]
    l[[i]] = aic
  }
  aic = do.call(rbind, l)
  aic[, mot := factor(mot, levels=rev(c("Curvilinear", "Straight line", "Average path")))] 
  aic = merge(aic[df==8], aic[df==11], by = c('mot'))
  setnames(aic, old = c('AICc.x', 'AICc.y'), new = c('simple', 'interactions'))
  aic$df.x = aic$df.y = NULL
  aic[,simple:=round(simple,2)]
  aic[,interactions:=round(interactions,2)]
  aic[, deltaAICc:=round(interactions-simple,2)]
  aic[, prob := round(exp(-0.5*deltaAICc)/sum(exp(-0.5*deltaAICc)),2)]
  aic[, ER := round((1-prob)/prob, 2)]    
  aic = aic[order(mot, decreasing = TRUE)]   
  fwrite(aic, file = 'Outputs/Table_SintAIC.csv') 
    
# Table Sxx
  effects_xx = c('Head', 'Midpiece', 'Tail')
  l = list()
  for(i in unique(adwl$mot)){
    #i ='Curvilinear'
    ai = adwl[mot == i]
    ms = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Head)+scale(Midpiece)+scale(Tail), ai)
    bsim = sim(ms, n.sim=nsim) 
    v = format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.5)))[5:7],2), nsmall = 2)
    lwr = format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.025)))[5:7],2), nsmall = 2)
    upr = format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.975)))[5:7],2), nsmall = 2)
    ii = data.table(motility=i, effect=effects_xx,estimate=v, lwr=lwr, upr=upr)

    mh = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Head), ai)
    bsim = sim(mh, n.sim=nsim) 
    ii[effect == 'Head', estimate_s := format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.5)))[5],2), nsmall=2)]
    ii[effect == 'Head', lwr_s := format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.025)))[5],2), nsmall=2)]
    ii[effect == 'Head', upr_s := format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.975)))[5],2), nsmall=2)]
    
    mm = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Midpiece), ai)
    bsim = sim(mm, n.sim=nsim) 
    ii[effect == 'Midpiece', estimate_s := format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.5)))[5],2), nsmall=2)]
    ii[effect == 'Midpiece', lwr_s := format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.025)))[5],2), nsmall=2)]
    ii[effect == 'Midpiece', upr_s := format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.975)))[5],2), nsmall=2)]
    
    mt = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Tail), ai)
    bsim = sim(mt, n.sim=nsim) 
    ii[effect == 'Tail', estimate_s := format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.5)))[5],2), nsmall=2)]
    ii[effect == 'Tail', lwr_s := format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.025)))[5],2), nsmall=2)]
    ii[effect == 'Tail', upr_s := format(round(c(apply(bsim@coef, 2, quantile, prob=c(0.975)))[5],2), nsmall=2)]
    
    l[[i]]=ii
  }
  x = do.call(rbind, l)
  x[, motility := factor(motility, levels=(c("Curvilinear", "Straight line", "Average path")))] 
  x = x[order(motility)]
  x[, es_ci_i := paste0(estimate, ' (', lwr,' -',upr,')')]
  x[, es_ci_s := paste0(estimate_s, ' (', lwr_s,' -',upr_s,')')]
  fwrite(x[,.(motility, es_ci_i, es_ci_s, file = 'Outputs/Table_Sxx.csv')

# Figure Sxx
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

  width__ = 0.4 
  cols_=pal_jco()(3)
  g = 
  ggplot(y, aes(y = effect, x = estimate, fill = motility, col = motility, shape = Model)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width = width_)) +
    
    scale_x_continuous(limits = c(-0.25, 0.35), expand = c(0, 0))+
    #scale_color_jco(name = 'Motility', guide = guide_legend(reverse = TRUE))+
    #scale_fill_jco(name = 'Motility', guide = guide_legend(reverse = TRUE))+
    scale_color_jco(name = 'Motility', guide = guide_legend(reverse = TRUE, order = 1,nrow=3,byrow=TRUE))+
    scale_fill_jco(name = 'Motility', guide = guide_legend(reverse = TRUE, order = 1,nrow=3,byrow=TRUE))+
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
    ggsave('Outputs/Fig_Sxx.png',g, width = 7/(5/7), height =4/(5/7), units = 'cm', bg="transparent", dpi = 600)
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
  llvx[, Motility:=response]
  llvx[, Motility := factor(Motility, levels=rev(c("Curvilinear", "Straight line", "Average path")))] 
  llvx[, trait := factor(trait, levels=rev(c(c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel'))))] 

  llvpx = data.table(do.call(rbind,lvpx) ) 
  llvpx[ , Motility:=response]
  llvpx[, Motility := factor(Motility, levels=(c("Curvilinear", "Straight line", "Average path")))] 
  llvpx[, trait := factor(trait, levels=rev(c(c('Acrosome', 'Nucleus', 'Midpiece', 'Tail', 'Total', 'Head','Flagellum', 'Midpiece_rel', 'Flagellum_rel'))))] 

# Table Smm
  t = copy(llvx)
  t[effect=='pred', effect := trait]
  t[, Motility := factor(Motility, levels=(c("Curvilinear", "Straight line", "Average path")))] 
  t[, trait := factor(trait, levels=(c('Acrosome', 'Nucleus','Midpiece', 'Tail','Total', 'Head', 'Flagellum', 'Midpiece_rel', 'Flagellum_rel')))] 
  t = t[order(trait,Motility)]


  t[,estimate := format(round(estimate,2), nsmall = 2)]
  t[,lwr := format(round(lwr,2), nsmall = 2)]
  t[,upr := format(round(upr,2), nsmall = 2)]
  #t[,CI := paste0(lwr,' -',upr)]
  t[, es_ci := paste0(estimate, ' (', lwr,' -',upr,')')]
  

  t_w = reshape(t[,.(trait, Motility, effect, es_ci)], idvar = c('trait','effect'), timevar = 'Motility', direction = "wide") 
  fwrite(t_w, file = 'Outputs/Table_Smm.csv')

# plot effect sizes
  llvx_ = llvx[effect == 'pred']
  llvx_[trait == 'Midpiece_rel', trait:='Midpiece\n(relative)']
  llvx_[trait == 'Flagellum_rel', trait:='Flagellum\n(relative)']
  llvx_[, trait := factor(trait, levels=rev(c("Acrosome", "Nucleus", "Midpiece","Tail","Total","Head", "Flagellum","Midpiece\n(relative)","Flagellum\n(relative)")))] 
  
  gE = 
    ggplot(llvx_, aes(y = trait, x = estimate, shape = Motility, col = Motility)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width =width_)) +
    scale_x_continuous(limits = c(-.5, .5), expand = c(0, 0))+
    scale_color_jco()+
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
  ggsave('Outputs/Fig_Ma_width-50mnm.png',gE, width = 4/(5/7), height =10, units = 'cm', bg="white", dpi = 600)     

  # viridis col
    col_ = c(viridis(1, alpha = 1, begin = 0.3, end = 0.3, direction = 1, option = "D"),
         viridis(1, alpha = 1, begin = 0.6, end = 0.6, direction = 1, option = "D"),
         viridis(1, alpha = 1, begin = 0.9, end = 0.9, direction = 1, option = "D")
         )

    col_ = c(viridis(1, alpha = 1, begin = 0.2, end = 0.2, direction = 1, option = "D"),
         viridis(1, alpha = 1, begin = 0.5, end = 0.5, direction = 1, option = "D"),
         viridis(1, alpha = 1, begin = 0.9, end = 0.9, direction = 1, option = "D")
         )

    #show_col(col_) 
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
  ggsave('Outputs/Fig_Ma_width-50mnm_viridis_v2.png',gEvir, width = 4/(5/7), height =10, units = 'cm', bg="white", dpi = 600)     

# plot predictions with raw data Fig ER v2 - x-axis labels - illustrations
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
    geom_point(aes(col = Morph, fill =Morph), pch =21, alpha = 0.8)+
    stat_cor(method="pearson",size = 2.75, cor.coef.name = 'r',aes(label = ..r.label..), label.x.npc = 'left', label.y.npc = 'bottom') +
    geom_ribbon(data = llvpx_a, aes(x=Length_avg, ymin=lwr, ymax=upr), fill = 'grey30', alpha = 0.2, show.legend = NA)+
    geom_line(data = llvpx_a, aes(x = Length_avg, y =pred), col =line_col)+
    facet_wrap(~Motility, scales = 'free_y', ncol = 1) +
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
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
   
# mix the two & export
  blank = ggplot() + theme_void() 
  gB = ggarrange(blank, ggR, nrow=2, heights=c(7.16-5.7,5.7))
  gAll = ggarrange(gE, gB, ncol=2, widths=c(4,15))  

  #ggsave('Outputs/Fig_4_width-190mm.png',gAll, width = 19/(5/7), height =10, units = 'cm', bg="white", dpi = 600)

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
  ggsave('Outputs/Fig_4_width-190mm_illust_v3_red.png',g_anot, width = 19/(5/7), height =10, units = 'cm', bg="white", dpi = 600)
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
  ggsave('Outputs/Fig_4_width-190mm_illust_viridis_red.png',g_anot, width = 19/(5/7), height =10, units = 'cm', bg="white", dpi = 600)


# END