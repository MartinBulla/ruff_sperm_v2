# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(ggpubr) 

  # constants
    set.seed = 5
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    fae = '#d4b691' # 'ffd6af'
    sat = 'white'
    ind = '#303030'
    colors = c(ind,sat,fae)
  
  # functions

# DATA
  source(here::here('R/DAT_prepare.R'))       

# Prepare estimates
  # motility - only one measure per bird
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
      ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
      lvx[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

      # get predictions
      m = lm(VAP ~ motileCount_ln + Morph, ddx)
      bsim = sim(m, n.sim=nsim) 
      v = apply(bsim@coef, 2, quantile, prob=c(0.5))
      newD=data.frame(motileCount_ln = mean(ddx$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
      X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
      newD$pred <-(X%*%v) 
      predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
      for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                  predmatrix[predmatrix < 0] <- 0
                  newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                  newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                  #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
      newD$motility = 'Average-path (VAP)'
      lvpx[['vap']] = data.table(newD)   
    # VSL
      m = lm(scale(VSL) ~ scale(log(motileCount))  + Morph, ddx)
      #summary(m)
      #plot(allEffects(m))
      bsim = sim(m, n.sim=nsim) 
      v = apply(bsim@coef, 2, quantile, prob=c(0.5))
      ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
      lvx[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

      # get predictions
      m = lm(VSL ~ motileCount_ln + Morph, ddx)
      bsim = sim(m, n.sim=nsim) 
      v = apply(bsim@coef, 2, quantile, prob=c(0.5))
      newD=data.frame(motileCount_ln = mean(ddx$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
      X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
      newD$pred <-(X%*%v) 
      predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
      for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                  predmatrix[predmatrix < 0] <- 0
                  newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                  newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                  #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
      newD$motility = 'Straight-line (VSL)'
      lvpx[['VSL']] = data.table(newD)
    # VCL
      m = lm(scale(VCL) ~scale(log(motileCount)) + Morph, ddx)
      #summary(m)
      #plot(allEffects(m))
      bsim = sim(m, n.sim=nsim) 
      v = apply(bsim@coef, 2, quantile, prob=c(0.5))
      ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
      lvx[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

      # get predictions
      m = lm(VCL ~ motileCount_ln + Morph, ddx)
      bsim = sim(m, n.sim=nsim) 
      v = apply(bsim@coef, 2, quantile, prob=c(0.5))
      newD=data.frame(motileCount_ln = mean(ddx$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
      X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
      newD$pred <-(X%*%v) 
      predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
      for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                  predmatrix[predmatrix < 0] <- 0
                  newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                  newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                  #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
      newD$motility = 'Curvilinear (VCL)'
      lvpx[['VCL']] = data.table(newD) 
               
    llvx = data.table(do.call(rbind,lvx) ) 
    llvx = llvx[effect != 'scale(log(motileCount))' ]
    llvx[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
    llvx[effect == '(Intercept)', effect:='Independent\n(Intercept)']
    llvx[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
    llvx[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
    llvx[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

    llvpx = data.table(do.call(rbind,lvpx) ) 
    llvpx[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
    llvpx[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
  



    # motility
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lv = list()
      lvp =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      #dd = a[part =='Acrosome']
      #dd[, motileCount_ln:=scale(log(motileCount))]
      # VAP
        m = lm(scale(VAP) ~ scale(log(motileCount)) + Morph, dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lv[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VAP ~ motileCount_ln + Morph, dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvp[['vap']] = data.table(newD)   
      # VSL
        m = lm(scale(VSL) ~ scale(log(motileCount))  + Morph, dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lv[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VSL ~ motileCount_ln + Morph, dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvp[['VSL']] = data.table(newD)
      # VCL
        m = lm(scale(VCL) ~scale(log(motileCount)) + Morph, dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lv[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VCL ~ motileCount_ln + Morph, dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvp[['VCL']] = data.table(newD) 
                 
      llv = data.table(do.call(rbind,lv) ) 
      llv = llv[effect != 'scale(log(motileCount))' ]
      llv[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llv[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llv[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llv[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      
      llv[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvp = data.table(do.call(rbind,lvp) ) 
      llvp[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvp[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
# motility 3 - mixed model
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvq = list()
      lvpq =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvq[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln + Morph+ (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpq[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))  + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvq[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpq[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount)) + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvq[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln + Morph + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpq[['VCL']] = data.table(newD) 
                 
      llvq = data.table(do.call(rbind,lvq) ) 
      llvq = llvq[effect != 'scale(log(motileCount))' ]
      llvq[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvq[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvq[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvq[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvq[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpq = data.table(do.call(rbind,lvpq) ) 
      llvpq[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpq[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 4 - mixed model with month
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvm = list()
      lvpm =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + month + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvm[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln+ month + Morph+ (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpm[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))+ month  + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvm[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph+ month + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5,Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpm[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount))+ month + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvm[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln+ month + Morph + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpm[['VCL']] = data.table(newD) 
                 
      llvm = data.table(do.call(rbind,lvm) ) 
      llvm = llvm[!effect %in% c('scale(log(motileCount))','monthMay') ]
      llvm[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvm[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvm[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvm[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvm[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpm = data.table(do.call(rbind,lvpm) ) 
      llvpm[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpm[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 5 - mixed model, only data without issues
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvi = list()
      lvpi =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      ddi = dd[issues == 'zero']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + Morph + (1|bird_ID), ddi)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvi[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln + Morph+ (1|bird_ID), ddi)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddi$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpi[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))  + Morph + (1|bird_ID), ddi)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvi[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph + (1|bird_ID), ddi)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddi$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpi[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount)) + Morph + (1|bird_ID), ddi)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvi[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln + Morph + (1|bird_ID), ddi)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddi$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpi[['VCL']] = data.table(newD)              
      llvi = data.table(do.call(rbind,lvi) ) 
      llvi = llvi[effect != 'scale(log(motileCount))' ]
      llvi[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvi[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvi[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvi[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvi[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpi = data.table(do.call(rbind,lvpi) ) 
      llvpi[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpi[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 6- mixed model with month and issuee
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvs = list()
      lvps =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      dd[issues == 'zero', issue:='no']
      dd[issues != 'zero', issue:='yes']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + month + issue +  Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvs[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln+ month + issue + Morph+ (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, issue = 0.5,  Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + issue + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvps[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))+ month + issue  + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvs[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph+ month + issue + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, issue = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + issue + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvps[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount))+ month + issue + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvs[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln+ month + Morph + issue + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, issue = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + issue + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvps[['VCL']] = data.table(newD) 
                 
      llvs = data.table(do.call(rbind,lvs) ) 
      llvs = llvs[!effect %in% c('scale(log(motileCount))','monthMay','issueyes') ]
      llvs[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvs[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvs[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvs[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvs[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvps = data.table(do.call(rbind,lvps) ) 
      llvps[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvps[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
      

    # morpho  
      # for averages
        l = list()
        lp =list()
        lpr = list()
        lcv = list()
        lpcv = list()
        
        for(i in unique(a$part)){
          #i ='Nucleus'
          m = lm(scale(Length_avg) ~ Morph, a[part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
          l[[i]]=data.frame(response=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

          # get predictions
          m = lm(Length_avg ~ Morph, a[part == i])
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          newD=data.frame(Morph = unique(a$Morph)) # values to predict for
          X <- model.matrix(~ Morph,data=newD) # exactly the model which was used has to be specified here 
          newD$Length_avg <-(X%*%v) 
          predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
          for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                          predmatrix[predmatrix < 0] <- 0
                          newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                          newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                          #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
          newD$part=i
          lp[[i]] = data.table(newD)

          print(i)     
          }          
        for(i in unique(ar$part)){
          if(i == 'Midpiece'){ii = 'Midpiece_rel'}
          if(i == 'Flagellum'){ii = 'Flagellum_rel'}
          m = lm(scale(Length_rel) ~ Morph, ar[part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
          l[[ii]]=data.frame(response=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

          # get predictions
          m = lm(Length_rel ~ Morph, ar[part == i])
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          newD=data.frame(Morph = unique(a$Morph)) # values to predict for
          X <- model.matrix(~ Morph,data=newD) # exactly the model which was used has to be specified here 
          newD$Length_rel <-(X%*%v) 
          predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
          for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                          predmatrix[predmatrix < 0] <- 0
                          newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                          newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                          #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
          newD$part=i
          lpr[[i]] = data.table(newD)

          print(i)     
          }          
        for(i in unique(cv_$part)){
          #i ='Nucleus'
          m = lm(scale(CV) ~ Morph, cv_[part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
          lcv[[i]]=data.frame(response=paste('CV',i),effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

          # get predictions
          m = lm(CV ~ Morph, cv_[part == i])
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          newD=data.frame(Morph = unique(a$Morph)) # values to predict for
          X <- model.matrix(~ Morph,data=newD) # exactly the model which was used has to be specified here 
          newD$Length_avg <-(X%*%v) 
          predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
          for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                          predmatrix[predmatrix < 0] <- 0
                          newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                          newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                          #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
          newD$part=i
          lpcv[[i]] = data.table(newD)

          print(i)     
          }          
        
        ll = data.table(do.call(rbind,l) ) 
        ll[, response := factor(response, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        ll[, effect := factor(effect, levels=rev(c("(Intercept)", "MorphSatellite", "MorphFaeder")))]
        ll[, unit := 'male average']

        llcv = data.table(do.call(rbind,lcv) ) 
        llcv[, response := factor(response, levels=rev(c("CV Acrosome", "CV Nucleus", "CV Head", "CV Midpiece","CV Tail","CV Flagellum","CV Total")))] 
        llcv[, effect := factor(effect, levels=rev(c("(Intercept)", "MorphSatellite", "MorphFaeder")))]
        
        llp = data.table(do.call(rbind,lp) ) 
        llp[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total")))] 
        llp[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

        llpr = data.table(do.call(rbind,lpr) ) 
        llpr[, part := factor(part, levels=rev(c("Midpiece","Flagellum")))] 
        llpr[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

      # for single values 
       

        ls = list()
        lps = list()
        lpsr = list()
        for(i in unique(a$part)){
          #i ='Nucleus'
          m = lmer(scale(Length_µm) ~ Morph + (1|bird_ID), b[part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=5000) 
          v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
          ls[[i]]=data.frame(response=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])
         
          # get predictions
          m = lmer(Length_µm ~ Morph + (1|bird_ID), b[part == i])
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          newD=data.frame(Morph = unique(a$Morph)) # values to predict for
          X <- model.matrix(~ Morph,data=newD) # exactly the model which was used has to be specified here 
          newD$Length_µm <-(X%*%v) 
          predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
          for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                          predmatrix[predmatrix < 0] <- 0
                          newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                          newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                          #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
          newD$part=i
          lps[[i]] = data.table(newD)

          print(i)         
          }  
        for(i in unique(ar$part)){
          if(i == 'Midpiece'){ii = 'Midpiece_rel'}
          if(i == 'Flagellum'){ii = 'Flagellum_rel'}
          m = lmer(scale(Length_rel) ~ Morph + (1|bird_ID), br[part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=5000) 
          v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
          ls[[ii]]=data.frame(response=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])
         
          # get predictions
          m = lmer(Length_rel ~ Morph + (1|bird_ID), br[part == i])
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          newD=data.frame(Morph = unique(ar$Morph)) # values to predict for
          X <- model.matrix(~ Morph,data=newD) # exactly the model which was used has to be specified here 
          newD$Length_rel <-(X%*%v) 
          predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
          for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                          predmatrix[predmatrix < 0] <- 0
                          newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                          newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                          #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
          newD$part=i
          lpsr[[i]] = data.table(newD)

          print(i)     
          }          
        
        lls = data.table(do.call(rbind,ls) ) 
        lls[, response := factor(response, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        lls[, effect := factor(effect, levels=rev(c("(Intercept)", "MorphSatellite", "MorphFaeder")))]
        lls[, unit := 'single sperm']

        llps = data.table(do.call(rbind,lps) ) 
        llps[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        llps[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

        llpsr = data.table(do.call(rbind,lpsr) ) 
        llpsr[, part := factor(part, levels=rev(c("Midpiece","Flagellum")))] 
        llpsr[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

    # motility
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lv = list()
      lvp =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      #dd = a[part =='Acrosome']
      #dd[, motileCount_ln:=scale(log(motileCount))]
      # VAP
        m = lm(scale(VAP) ~ scale(log(motileCount)) + Morph, dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lv[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VAP ~ motileCount_ln + Morph, dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvp[['vap']] = data.table(newD)   
      # VSL
        m = lm(scale(VSL) ~ scale(log(motileCount))  + Morph, dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lv[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VSL ~ motileCount_ln + Morph, dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvp[['VSL']] = data.table(newD)
      # VCL
        m = lm(scale(VCL) ~scale(log(motileCount)) + Morph, dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lv[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VCL ~ motileCount_ln + Morph, dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvp[['VCL']] = data.table(newD) 
                 
      llv = data.table(do.call(rbind,lv) ) 
      llv = llv[effect != 'scale(log(motileCount))' ]
      llv[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llv[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llv[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llv[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      
      llv[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvp = data.table(do.call(rbind,lvp) ) 
      llvp[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvp[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 2 - only one measure per bird
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvx = list()
      lvpx =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      dd1 = dd[month == 'June']
      dd2 = dd[month == 'May']
      ddx = rbind(dd1,dd2[!bird_ID%in%dd1$bird_ID])
      # VAP
        m = lm(scale(VAP) ~ scale(log(motileCount)) + Morph, ddx)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lvx[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VAP ~ motileCount_ln + Morph, ddx)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddx$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpx[['vap']] = data.table(newD)   
      # VSL
        m = lm(scale(VSL) ~ scale(log(motileCount))  + Morph, ddx)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lvx[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VSL ~ motileCount_ln + Morph, ddx)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddx$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpx[['VSL']] = data.table(newD)
      # VCL
        m = lm(scale(VCL) ~scale(log(motileCount)) + Morph, ddx)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lvx[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VCL ~ motileCount_ln + Morph, ddx)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddx$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpx[['VCL']] = data.table(newD) 
                 
      llvx = data.table(do.call(rbind,lvx) ) 
      llvx = llvx[effect != 'scale(log(motileCount))' ]
      llvx[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvx[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvx[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvx[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvx[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpx = data.table(do.call(rbind,lvpx) ) 
      llvpx[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpx[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 3 - mixed model
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvq = list()
      lvpq =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvq[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln + Morph+ (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpq[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))  + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvq[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpq[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount)) + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvq[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln + Morph + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpq[['VCL']] = data.table(newD) 
                 
      llvq = data.table(do.call(rbind,lvq) ) 
      llvq = llvq[effect != 'scale(log(motileCount))' ]
      llvq[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvq[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvq[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvq[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvq[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpq = data.table(do.call(rbind,lvpq) ) 
      llvpq[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpq[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 4 - mixed model with month
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvm = list()
      lvpm =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + month + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvm[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln+ month + Morph+ (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpm[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))+ month  + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvm[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph+ month + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5,Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpm[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount))+ month + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvm[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln+ month + Morph + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpm[['VCL']] = data.table(newD) 
                 
      llvm = data.table(do.call(rbind,lvm) ) 
      llvm = llvm[!effect %in% c('scale(log(motileCount))','monthMay') ]
      llvm[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvm[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvm[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvm[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvm[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpm = data.table(do.call(rbind,lvpm) ) 
      llvpm[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpm[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 5 - mixed model, only data without issues
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvi = list()
      lvpi =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      ddi = dd[issues == 'zero']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + Morph + (1|bird_ID), ddi)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvi[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln + Morph+ (1|bird_ID), ddi)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddi$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpi[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))  + Morph + (1|bird_ID), ddi)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvi[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph + (1|bird_ID), ddi)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddi$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpi[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount)) + Morph + (1|bird_ID), ddi)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvi[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln + Morph + (1|bird_ID), ddi)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddi$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpi[['VCL']] = data.table(newD)              
      llvi = data.table(do.call(rbind,lvi) ) 
      llvi = llvi[effect != 'scale(log(motileCount))' ]
      llvi[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvi[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvi[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvi[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvi[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpi = data.table(do.call(rbind,lvpi) ) 
      llvpi[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpi[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 6- mixed model with month and issuee
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvs = list()
      lvps =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      dd[issues == 'zero', issue:='no']
      dd[issues != 'zero', issue:='yes']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + month + issue +  Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvs[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln+ month + issue + Morph+ (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, issue = 0.5,  Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + issue + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvps[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))+ month + issue  + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvs[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph+ month + issue + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, issue = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + issue + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvps[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount))+ month + issue + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvs[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln+ month + Morph + issue + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, issue = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + issue + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvps[['VCL']] = data.table(newD) 
                 
      llvs = data.table(do.call(rbind,lvs) ) 
      llvs = llvs[!effect %in% c('scale(log(motileCount))','monthMay','issueyes') ]
      llvs[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvs[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvs[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvs[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvs[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvps = data.table(do.call(rbind,lvps) ) 
      llvps[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvps[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    


# perhaps useful
    require(arm) 
    require(effects)
    require(ggpubr)
    require(ggsci) 
    require(grid)
    require(gridExtra)
    require(magrittr)
    require(multcomp)
    require(PerformanceAnalytics)
    require(rptR) 
    require(stringi)
    require(viridis)
    require(readxl)

# END