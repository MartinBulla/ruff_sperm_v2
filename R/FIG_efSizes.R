# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  require(ggpubr) 
  require(ggsci)

  width_ = 0.4 # spacing between error bars

# DATA
  source(here::here('R/DAT_prepare.R'))       

# Prepare estimates
  effects_ = c('Satellite relative to independent','Faeder relative to independent', 'Faeder relative to satellite')
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
      
      mb = data.table(bsim@coef)
      names(mb) = c('int','n','s','f')
      mb[, FrelS:=(int+f)-(int+s)]
      v = c(apply(bsim@coef, 2, quantile, prob=c(0.5))[3:4], median(mb$FrelS)) 
      lwr = c(apply(bsim@coef, 2, quantile, prob=c(0.025))[3:4], quantile(mb$FrelS, prob = 0.025)) 
      upr = c(apply(bsim@coef, 2, quantile, prob=c(0.975))[3:4], quantile(mb$FrelS, prob = 0.975)) 
      
      lvx[['VAP']]=data.frame(response='Average path',effect=effects_,estimate=v, lwr=lwr, upr=upr)

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
      newD$motility = 'Average path'
      lvpx[['vap']] = data.table(newD)   
    # VSL
      m = lm(scale(VSL) ~ scale(log(motileCount))  + Morph, ddx)
      #summary(m)
      #plot(allEffects(m))
      bsim = sim(m, n.sim=nsim) 
      mb = data.table(bsim@coef)
      names(mb) = c('int','n','s','f')
      mb[, FrelS:=(int+f)-(int+s)]
      v = c(apply(bsim@coef, 2, quantile, prob=c(0.5))[3:4], median(mb$FrelS)) 
      lwr = c(apply(bsim@coef, 2, quantile, prob=c(0.025))[3:4], quantile(mb$FrelS, prob = 0.025)) 
      upr = c(apply(bsim@coef, 2, quantile, prob=c(0.975))[3:4], quantile(mb$FrelS, prob = 0.975)) 
      
      lvx[['VSL']]=data.frame(response='Straight line',effect=effects_,estimate=v, lwr=lwr, upr=upr)

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
      newD$motility = 'Straight line'
      lvpx[['VSL']] = data.table(newD)
    # VCL
      m = lm(scale(VCL) ~scale(log(motileCount)) + Morph, ddx)
      #summary(m)
      #plot(allEffects(m))
      bsim = sim(m, n.sim=nsim) 
      mb = data.table(bsim@coef)
      names(mb) = c('int','n','s','f')
      mb[, FrelS:=(int+f)-(int+s)]
      v = c(apply(bsim@coef, 2, quantile, prob=c(0.5))[3:4], median(mb$FrelS)) 
      lwr = c(apply(bsim@coef, 2, quantile, prob=c(0.025))[3:4], quantile(mb$FrelS, prob = 0.025)) 
      upr = c(apply(bsim@coef, 2, quantile, prob=c(0.975))[3:4], quantile(mb$FrelS, prob = 0.975)) 
      
      lvx[['VCL']]=data.frame(response='Curvilinear',effect=effects_,estimate=v, lwr=lwr, upr=upr)

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
      newD$motility = 'Curvilinear'
      lvpx[['VCL']] = data.table(newD) 
               
    llvx = data.table(do.call(rbind,lvx) ) 
    llvx[, effect := factor(effect, levels=c("Faeder relative to satellite","Faeder relative to independent","Satellite relative to independent"))] 
    llvx[, response := factor(response, levels=rev(c("Curvilinear", "Straight line", "Average path")))] 

    llvpx = data.table(do.call(rbind,lvpx) ) 
    llvpx[, motility := factor(motility, levels=(c("Curvilinear", "Straight line", "Average path")))] 
    llvpx[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
  # morpho - averages
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
      mb = data.table(bsim@coef)
      names(mb) = c('int','s','f')
      mb[, FrelS:=(int+f)-(int+s)]
      mb$int = NULL
      v = apply(mb, 2, quantile, prob=c(0.5))
      ci = apply(mb, 2, quantile, prob=c(0.025,0.975)) 
      l[[i]]=data.frame(response=i,effect=effects_,estimate=v, lwr=ci[1,], upr=ci[2,])

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
      mb = data.table(bsim@coef)
      names(mb) = c('int','s','f')
      mb[, FrelS:=(int+f)-(int+s)]
      mb$int = NULL
      v = apply(mb, 2, quantile, prob=c(0.5))
      ci = apply(mb, 2, quantile, prob=c(0.025,0.975)) 
      l[[ii]]=data.frame(response=ii,effect=effects_,estimate=v, lwr=ci[1,], upr=ci[2,])

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
      mb = data.table(bsim@coef)
      names(mb) = c('int','s','f')
      mb[, FrelS:=(int+f)-(int+s)]
      mb$int = NULL
      v = apply(mb, 2, quantile, prob=c(0.5))
      ci = apply(mb, 2, quantile, prob=c(0.025,0.975)) 
      lcv[[i]]=data.frame(response=i,effect=effects_,estimate=v, lwr=ci[1,], upr=ci[2,])

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
    ll[, effect := factor(effect, levels=c("Faeder relative to satellite","Faeder relative to independent","Satellite relative to independent"))] 
    ll[, unit := 'male average']
    ll_ = ll[response %in%c("Acrosome", "Nucleus","Midpiece","Tail","Total")]

    llcv = data.table(do.call(rbind,lcv) ) 
    llcv[, response := factor(response, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total")))] 
    llcv[, effect := factor(effect, c("Faeder relative to satellite","Faeder relative to independent","Satellite relative to independent"))] 
    llcv_ = llcv[response %in%c("Acrosome", "Nucleus","Midpiece","Tail","Total")]
    
    llp = data.table(do.call(rbind,lp) ) 
    llp[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total")))] 
    llp[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

    llpr = data.table(do.call(rbind,lpr) ) 
    llpr[, part := factor(part, levels=rev(c("Midpiece","Flagellum")))] 
    llpr[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

# not used - Fig E v1
  gV = 
    ggplot(llvx, aes(y = response, x = estimate, shape = effect, col = effect)) +  
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width = width_)) +
    scale_x_continuous(limits = c(-1.5, 2), expand = c(0, 0))+
    scale_color_jco(guide = guide_legend(reverse = TRUE))+
    scale_shape(guide = guide_legend(reverse = TRUE))+
    labs(y = NULL, x = "Standardized effect size", subtitle = 'Motility')+
    theme_bw() +
    theme(plot.subtitle = element_text(size=9, color = 'grey30'),
        legend.title = element_blank(),
        legend.text=element_text(size=7),
        legend.key.height= unit(0.75,"line"),
        legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(-10,-10,-10,-10),
        #legend.position = "none"
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank(),

        plot.margin = margin(3,3,1,1, "mm")
        )  

  gM = 
    ggplot(ll_, aes(y = response, x = estimate, shape = effect, col = effect)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width = width_)) +
    scale_x_continuous(limits = c(-1.5, 2), expand = c(0, 0))+
    scale_color_jco(guide = guide_legend(reverse = TRUE))+
    scale_shape(guide = guide_legend(reverse = TRUE))+
    labs(y = NULL, x = "Standardized effect size", subtitle = 'Morphology')+
    theme_bw() +
    theme(plot.subtitle = element_text(size=9, color = 'grey30'),
        legend.title = element_blank(),
        legend.text=element_text(size=7),
        legend.key.height= unit(0.75,"line"),
        legend.margin=margin(0,0,0,0),
        #legend.position = "none"
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank(),

        plot.margin = margin(3,3,1,1, "mm")
        )  

  gCV = 
    ggplot(llcv_, aes(y = response, x = estimate, shape = effect, col = effect)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width =width_)) +
    scale_x_continuous(limits = c(-1.5, 2), expand = c(0, 0))+
    scale_color_jco(guide = guide_legend(reverse = TRUE))+
    scale_shape(guide = guide_legend(reverse = TRUE))+
    labs(y = NULL, x = "Standardized effect size", subtitle = 'Coefficient of variation')+
    theme_bw() +
    theme(
        plot.subtitle = element_text(size=9, color = 'grey30'),
        legend.title = element_blank(),
        legend.text=element_text(size=7),
        legend.key.height= unit(0.75,"line"),
        legend.margin=margin(0,0,0,0),
        #legend.position = "none"
        axis.title.x = element_text(size = 10, color ='grey10'),
        axis.ticks = element_blank(),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank(),

        plot.margin = margin(3,3,1,1, "mm")
        )    
  
  ggA = ggarrange(
    gV,gM,gCV, 
    nrow=3, heights=c(2.5,3.78,4.22), align = 'v'
    )  

  blank = ggplot() + geom_rect(aes(xmin = 1, xmax = 3, ymin = 10, ymax = 15), fill = 'white') +theme_transparent()
  blank_grob = ggplotGrob(blank)
  
  g_exp = ggA + annotation_custom(blank_grob, xmin=.55, xmax=1.1, ymin = 0, ymax = 0.7)

  ggsave('Outputs/Fig_E.png',g_exp, width = 8/(5/7), height =15, units = 'cm', bg="white", dpi = 600)
# Fig E v2 - legend top
  cols_=pal_jco()(3)
  gV = 
    ggplot(llvx, aes(y = response, x = estimate, shape = effect, col = effect)) +  
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width = width_)) +
    scale_x_continuous(limits = c(-1.5, 2), expand = c(0, 0))+
    scale_color_jco()+
    scale_shape(guide = guide_legend(reverse = TRUE))+
    labs(y = NULL, x = "Standardized effect size", subtitle = 'Motility')+
    guides(col=guide_legend(nrow=3,byrow=TRUE,reverse = TRUE),shape = guide_legend(nrow=3,byrow=TRUE,reverse = TRUE))+
    #annotate(geom="text", x=0.65, y=3.13, label="Satellite\nrelative to\nindependent", color=cols_[3],hjust = 0, size = 3.25) +
    #annotate(geom="text", x=0.65, y=2, label="Faeder\nrelative to\nindependent", color=cols_[2],hjust = 0, size = 3.25) +
    #annotate(geom="text", x=0.65, y=.87, label="Faeder\nrelative to\nsatellite", color=cols_[1],hjust = 0, size = 3.25) +
    theme_bw() +
    theme(plot.subtitle = element_text(size=9, color = 'grey30'),
        legend.title = element_blank(),
        legend.text=element_text(size=7),
        legend.key.height= unit(0.25,"line"),
        legend.margin=margin(0,0,0,0),
        legend.position=c(0.45,1.6),

        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank(),

        plot.margin = margin(15.5,3,1,1, "mm")
        )  
  
  gM = 
    ggplot(ll_, aes(y = response, x = estimate, shape = effect, col = effect)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width = width_)) +
    scale_x_continuous(limits = c(-1.5, 2), expand = c(0, 0))+
    scale_color_jco(guide = guide_legend(reverse = TRUE))+
    scale_shape(guide = guide_legend(reverse = TRUE))+
    labs(y = NULL, x = "Standardized effect size", subtitle = 'Morphology')+
    theme_bw() +
    theme(plot.subtitle = element_text(size=9, color = 'grey30'),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank(),

        plot.margin = margin(3,3,1,1, "mm")
        )  

  gCV = 
    ggplot(llcv_, aes(y = response, x = estimate, shape = effect, col = effect)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width =width_)) +
    scale_x_continuous(limits = c(-1.5, 2), expand = c(0, 0))+
    scale_color_jco(guide = guide_legend(reverse = TRUE))+
    scale_shape(guide = guide_legend(reverse = TRUE))+
    labs(y = NULL, x = "Standardized effect size", subtitle = 'Coefficient of variation')+
    theme_bw() +
    theme(
        plot.subtitle = element_text(size=9, color = 'grey30'),
        legend.position = "none",
        axis.title.x = element_text(size = 10, color ='grey10'),
        axis.ticks = element_blank(),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank(),

        plot.margin = margin(3,3,1,1, "mm")
        )    
  
  ggA = ggarrange(
    gV,gM,gCV, 
    nrow=3, heights=c(3.5,3.78,4.22), align = 'v'
    )
  ggsave('Outputs/Fig_E_width-50mnm.png',ggA, width = 5/(5/7), height =16, units = 'cm', bg="white", dpi = 600)


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
  # morphp
  # CV
      
# Further estimates


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