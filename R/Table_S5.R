# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  require(ggpubr) 
  require(ggsci)
# DATA
  source(here::here('R/DAT_prepare.R'))       

# Prepare estimates
  # motility - only one measure per bird (June)
    effects_ = c('Intercept (independent)','Motile count','Satellite relative to independent', 'Faeder relative to independent')
    lvx = list()
    lvpx =list()
    d[, motileCount_ln:=scale(log(motileCount))]
    dd = d[!Morph%in%'Zebra finch']
   
    # use June values and for 4 males without June, May
      dd1 = dd[month == 'June']
      dd2 = dd[month == 'May']
      ddx = rbind(dd1,dd2[!bird_ID%in%dd1$bird_ID])
    
    # VAP
      m = lm(scale(VAP) ~ scale(log(motileCount)) + Morph, ddx)
      #summary(m)
      #plot(allEffects(m))
      bsim = sim(m, n.sim=nsim) 
    
      v = apply(bsim@coef, 2, quantile, prob=c(0.5))
      lwr = apply(bsim@coef, 2, quantile, prob=c(0.025))
      upr = apply(bsim@coef, 2, quantile, prob=c(0.975))
      
      lvx[['VAP']]=data.frame(response='Average path',effect=effects_,estimate=v, lwr=lwr, upr=upr)      
    # VSL
      m = lm(scale(VSL) ~ scale(log(motileCount))  + Morph, ddx)
      #summary(m)
      #plot(allEffects(m)) 
      bsim = sim(m, n.sim=nsim) 
    
      v = apply(bsim@coef, 2, quantile, prob=c(0.5))
      lwr = apply(bsim@coef, 2, quantile, prob=c(0.025))
      upr = apply(bsim@coef, 2, quantile, prob=c(0.975))
      
      lvx[['VSL']]=data.frame(response='Straight line',effect=effects_,estimate=v, lwr=lwr, upr=upr)
    # VCL
      m = lm(scale(VCL) ~scale(log(motileCount)) + Morph, ddx)
      #summary(m)
      #plot(allEffects(m))
      bsim = sim(m, n.sim=nsim) 
    
      v = apply(bsim@coef, 2, quantile, prob=c(0.5))
      lwr = apply(bsim@coef, 2, quantile, prob=c(0.025))
      upr = apply(bsim@coef, 2, quantile, prob=c(0.975))

      lvx[['VCL']]=data.frame(response='Curvilinear',effect=effects_,estimate=v, lwr=lwr, upr=upr)            
    
    llvx = data.table(do.call(rbind,lvx) ) 
    llvx[, response := factor(response, levels=(c("Curvilinear", "Straight line", "Average path")))] 
    llvx = llvx[order(response)]
  # morpho - averages
    effects_ = c('Intercept (independent)','Satellite relative to independent', 'Faeder relative to independent')
    l = list()
    lcv = list()
    
    for(i in unique(a$part)){
      #i ='Nucleus'
      m = lm(scale(Length_avg) ~ Morph, a[part == i])
      #summary(m)
      #plot(allEffects(m))
      bsim = sim(m, n.sim=nsim) 
      mb = data.table(bsim@coef)
      v = apply(mb, 2, quantile, prob=c(0.5))
      ci = apply(mb, 2, quantile, prob=c(0.025,0.975)) 
      l[[i]]=data.frame(response=i,effect=effects_,estimate=v, lwr=ci[1,], upr=ci[2,])
      print(i)     
      }          
    for(i in unique(ar$part)){
      if(i == 'Midpiece'){ii = 'Midpiece_rel'}
      if(i == 'Flagellum'){ii = 'Flagellum_rel'}
      m = lm(scale(Length_rel) ~ Morph, ar[part == i])
      #summary(m)
      #plot(allEffects(m))
      bsim = sim(m, n.sim=nsim) 
      mb = data.table(bsim@coef)
      v = apply(mb, 2, quantile, prob=c(0.5))
      ci = apply(mb, 2, quantile, prob=c(0.025,0.975)) 
      l[[ii]]=data.frame(response=ii,effect=effects_,estimate=v, lwr=ci[1,], upr=ci[2,])
      print(i)     
      }          
    for(i in unique(cv_$part)){
      #i ='Nucleus'
      m = lm(scale(CV) ~ Morph, cv_[part == i])
      #summary(m)
      #plot(allEffects(m))
      bsim = sim(m, n.sim=nsim) 
      mb = data.table(bsim@coef)
      v = apply(mb, 2, quantile, prob=c(0.5))
      ci = apply(mb, 2, quantile, prob=c(0.025,0.975)) 
      lcv[[i]]=data.frame(response=i,effect=effects_,estimate=v, lwr=ci[1,], upr=ci[2,])
      print(i)     
      }          
    
    ll = data.table(do.call(rbind,l) ) 
    ll[, response := factor(response, levels=(c("Acrosome", "Nucleus", "Midpiece","Tail","Total","Head","Flagellum","Midpiece_rel","Flagellum_rel")))] 
    ll = ll[order(response)]

    llcv = data.table(do.call(rbind,lcv) ) 
    llcv[, response := factor(response, levels=(c("Acrosome", "Nucleus",  "Midpiece","Tail","Total","Head","Flagellum")))] 
    llcv = llcv[order(response)]
    
# Export
  t = copy(llvx)
  t[effect=='pred', effect := trait]
  t[, response := factor(response, levels = (c("Curvilinear", "Straight line", "Average path")))]

  t[,estimate := format(round(estimate,2), nsmall = 2)]
  t[,lwr := format(round(lwr,2), nsmall = 2)]
  t[,upr := format(round(upr,2), nsmall = 2)]
  t[, es_ci := paste0(estimate, ' (', lwr,' -',upr,')')]
  
  t_w = reshape(t[,.(response, effect, es_ci)], idvar = c('effect'), timevar = 'response', direction = "wide") 
  fwrite(t_w, file = 'Outputs/Table_S5a.csv')

  tl = copy(ll)
  tc = copy(llcv)
  tl[,unit := 'length']
  tc[,unit := 'CV']
  tt = rbind(tl,tc)
  tt[,estimate := format(round(estimate,2), nsmall = 2)]
  tt[,lwr := format(round(lwr,2), nsmall = 2)]
  tt[,upr := format(round(upr,2), nsmall = 2)]
  tt[, es_ci := paste0(estimate, ' (', lwr,' -',upr,')')]
  tt_w = reshape(tt[,.(response, effect, unit, es_ci)], idvar = c('response','effect'), timevar = 'unit', direction = "wide") 
  fwrite(tt_w, file = 'Outputs/Table_S5b.csv')


  fwrite(rbind(llvx,ll,llcv), file = 'Outputs/Table_S1_long.csv')

# END