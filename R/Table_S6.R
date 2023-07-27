  # =============================================================
  # ❗ Runs relative to the project's root directory and
  # exports Table S6 into ./Outputs/
  # =============================================================

# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  require(ggnewscale)
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

  fwrite(t_w, file = 'Outputs/Table_S6.csv')

# END