# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  #require(facetscales)
  require(ggpubr) 
  require(ggsci)

  width_ = 0.4 # spacing between error bars

# DATA
  source(here::here('R/DAT_prepare.R'))       
  # prepare motility for plotting
    am = a[part=='Acrosome']
    aml = data.table(melt(am[,.(bird_ID,month,Morph,age,HL,motileCount,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph","age","HL","motileCount"), variable.name = "Motility"))
    aml[Motility == 'VAP' ,mot:='Average path']
    aml[Motility == 'VCL' ,mot:='Curvilinear']
    aml[Motility == 'VSL' ,mot:='Straight line']

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
      bsim = sim(object = m, n.sims = nsim)
      
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
  # morpho 
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

    llpcv = data.table(do.call(rbind,lpcv) ) 
    llpcv[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total")))] 
    llpcv[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

# (a) - legend top
  cols_=pal_jco()(3)
  
  gV = 
    ggplot(llvx, aes(y = response, x = estimate, shape = effect, col = effect)) +  
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width = width_)) +
    scale_x_continuous(limits = c(-1.5, 2), expand = c(0, 0))+
    scale_color_jco()+
    scale_shape(guide = guide_legend(reverse = TRUE))+
    labs(y = NULL, x = "Standardized effect size", subtitle = 'Velocity', tag = "(A)")+
    guides(col=guide_legend(nrow=3,byrow=TRUE,reverse = TRUE),shape = guide_legend(nrow=3,byrow=TRUE,reverse = TRUE))+
    #annotate(geom="text", x=0.65, y=3.13, label="Satellite\nrelative to\nindependent", color=cols_[3],hjust = 0, size = 3.25) +
    #annotate(geom="text", x=0.65, y=2, label="Faeder\nrelative to\nindependent", color=cols_[2],hjust = 0, size = 3.25) +
    #annotate(geom="text", x=0.65, y=.87, label="Faeder\nrelative to\nsatellite", color=cols_[1],hjust = 0, size = 3.25) +
    theme_bw() +
    theme(
        plot.tag.position = c(0.027, 0.98),
        plot.tag = element_text(face='bold', size = 10),
        plot.subtitle = element_text(size=9, color = 'grey30'),
        legend.title = element_blank(),
        legend.text=element_text(size=7.5, color = 'grey30'),
        legend.key.height= unit(0.2,"line"),
        legend.margin=margin(0,0,0,0),
        legend.position=c(0.5,1.6),

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
    labs(y = NULL, x = "Standardized effect size", subtitle = 'Morphology\n(length)')+
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
    labs(y = NULL, x = "Standardized effect size", subtitle = 'Morphology\n(coefficient of variation)')+
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
  #ggsave('Outputs/Fig_3a_width-50mnm.png',ggA, width = 5/(5/7), height =16, units = 'cm', bg="white", dpi = 600)
# (b) - x-axis labels - illustrations
  aml[, mot2 := factor(mot,levels=c('Curvilinear','Straight line','Average path','w','x'))]
  #aml[,summary(value), by = mot2]
  llvpx[, value:=pred]
  llvpx[, mot2:=motility]

  a_m = a[!part%in%c('Head','Flagellum')]
  a_m[, part := factor(part, levels=(c("Acrosome", "Nucleus", "Midpiece","Tail","Total")))] 
  llp_m = llp[!part%in%c('Head','Flagellum')]
  llp_m[, part := factor(part, levels=(c("Acrosome", "Nucleus", "Midpiece","Tail","Total")))] 

  cv_m =cv_[!part%in%c('Head','Flagellum')]
  cv_m[, part := factor(part, levels=(c("Acrosome", "Nucleus", "Midpiece","Tail","Total")))] 
  # for visualization purpose one extreme tail CV (and hence also Total CV) value brought down and labeled
    cv_m[part== 'Tail' & CV>15, CV:=8] 
    cv_m[part == 'Total' & CV>7.5, CV:=6] 
    #cv_m[,summary(CV), by = part]
  llpcv_ = llpcv[!part%in%c('Head','Flagellum')]
  llpcv_[, part := factor(part, levels=(c("Acrosome", "Nucleus", "Midpiece","Tail","Total")))] 
  llpcv_[,CV:=Length_avg]
  size_ =1.2
  gv =
  ggplot(aml, aes(x = Morph, y = value)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(),  aes(col = Morph, fill =Morph), dotsize = 1.1)+
    geom_boxplot(width = 0.25, col = 'grey50', outlier.shape = NA, fill = NA) + 
    geom_errorbar(data = llvpx, aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25), col = 'red' ) +
    geom_point(data = llvpx, aes(x = Morph, y =value), position = position_dodge(width = 0.25), col = 'red', size = size_) +
    facet_wrap(~mot2, scales = 'free_y', nrow = 1,drop=FALSE)+
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    scale_y_continuous('Velocity [μm/s]', expand = c(0, 0))+
    xlab('Morph') +
    labs(tag = '(B)')+
    guides(x =  guide_axis(angle = -15)) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.tag.position = c(0.006, 1),
      plot.tag = element_text(face='bold',size =10),

      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.x = element_blank(), axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_text(size = 7.5, color="grey20",  margin=margin(1,1,1,1,"mm")),
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
      
      panel.border = element_rect(color = 'grey70')
      )  

    # adjust scales for morpho - https://stackoverflow.com/questions/51735481/ggplot2-change-axis-limits-for-each-individual-facet-panel
    facet_bounds_v = data.frame(
      mot2 = c('Curvilinear','Average path', 'Straight line','w','x'),
      ymin = c(25, 10, 5,1,1), 
      ymax = c(80, 50, 35,5,5), 
      breaks = c(10,10,10,1,1), stringsAsFactors=FALSE)

    ff_v <- with(facet_bounds_v,
           data.frame(value=c(ymin,ymax),
                      mot2=c(mot2,mot2)))
    ff_v$mot2 = factor(ff_v$mot2, levels=(c("Curvilinear", "Average path", "Straight line",'w','x'))) 
    gv = gv + geom_point(data=ff_v,x=NA)
   # remove panels  
    ggv = ggplotGrob(gv)
    rm_grobs <- ggv$layout$name %in% c("panel-4-1", "panel-5-1","strip-t-4-1", "strip-t-5-1","axis-l-1-4", "axis-l-1-5")
    ggv$grobs[rm_grobs] <- NULL
    ggv$layout <- ggv$layout[!rm_grobs, ]
  gm=
  ggplot(a_m, aes(x = Morph, y = Length_avg)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(),  aes(col = Morph, fill =Morph), dotsize = 1.1)+
    geom_boxplot(width = 0.25, col = 'grey50', outlier.shape = NA, fill = NA) + 
    geom_errorbar(data = llp_m, aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25), col = 'red' ) +
    geom_point(data = llp_m, aes(x = Morph, y =Length_avg), position = position_dodge(width = 0.25), col = 'red', size = size_) +
    facet_wrap(~part, scales = 'free_y', nrow = 1,drop=FALSE)+
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    xlab('Morph') +
    scale_y_continuous('Length [µm]', expand = c(0, 0))+
    guides(x =  guide_axis(angle = -15)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.x = element_blank(), axis.text.x = element_blank(), 
      axis.ticks = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_text(size = 7.5, color="grey20",  margin=margin(1,1,1,1,"mm")),
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
      
      panel.border = element_rect(color = 'grey70')
      ) 
    # adjust scales for morpho - https://stackoverflow.com/questions/51735481/ggplot2-change-axis-limits-for-each-individual-facet-panel
    facet_bounds <- read.table(header=TRUE,
        text=                           
        "part ymin ymax breaks
        Acrosome     3     5    .5
        Nucleus     24     32    2
        Midpiece     20     28    2
        Tail     70     95    5
        Total     125    150    5",
        stringsAsFactors=FALSE)

    ff <- with(facet_bounds,
           data.frame(Length_avg=c(ymin,ymax),
                      part=c(part,part)))
    ff$part = factor(ff$part, levels=(c("Acrosome", "Nucleus", "Midpiece","Tail","Total"))) 
    gm = gm + geom_point(data=ff,x=NA)
  gcv= 
  ggplot(cv_m, aes(x = Morph, y = CV)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(),  aes(col = Morph, fill =Morph), dotsize = 1.1)+
    geom_boxplot(width = 0.25, col = 'grey50', outlier.shape = NA, fill = NA) + 
    geom_errorbar(data = llpcv_, aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25), col = 'red' ) +
    geom_point(data = llpcv_, aes(x = Morph, y =CV), position = position_dodge(width = 0.25), col = 'red', size = size_) +
    facet_wrap(~part, scales = 'free_y', nrow = 1,drop=FALSE)+
    scale_color_manual(values=cols)+ 
    scale_fill_manual(values=fills)+
    xlab('Morph') +
    scale_y_continuous('Coefficient of variation', expand = c(0, 0))+
    guides(x =  guide_axis(angle = -15)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10, , colour="grey10"),
      axis.title.x = element_blank(), 
      axis.text.x = element_blank(), 
      axis.text.y = element_text(margin = margin(r = -1)),
      axis.ticks = element_blank(),
      strip.text = element_blank(),
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
      
      plot.margin = margin(b=10.4,1,1,1, "mm"),
      panel.border = element_rect(color = 'grey70')
      ) 

    # adjust scales for morpho - https://stackoverflow.com/questions/51735481/ggplot2-change-axis-limits-for-each-individual-facet-panel
    facet_bounds_cv <- read.table(header=TRUE,
        text=                           
        "part ymin ymax breaks
        Acrosome     0     30    10
        Nucleus     2     12    2
        Midpiece     1     7    1
        Tail     0     8    2
        Total     1    6    1",
        stringsAsFactors=FALSE)

    ff_cv <- with(facet_bounds_cv,
           data.frame(CV=c(ymin,ymax),
                      part=c(part,part)))

    gcv = gcv + geom_point(data=ff_cv,x=NA) 

  ggR = ggarrange(
    ggv,
    gm,
    gcv,  
    nrow=3, align = 'v', heights=c(3,3,3.2)  
    ) 
  right = 0.38
  right2 = 0.01
  ymin_ = -0.92
    
  ggExp =  
  ggR + 
    #annotation_custom(gi, xmin=0.064, xmax=0.116, ymin=ymin_) + 
    #annotation_custom(gs, xmin=0.064+0.05, xmax=0.116+0.05, ymin=ymin_) + 
    #annotation_custom(gf, xmin=0.064+0.1, xmax=0.116+0.1, ymin=ymin_)#+

    annotation_custom(gi, xmin=0.064+right, xmax=0.116+right, ymin=ymin_) + 
    annotation_custom(gs, xmin=0.064+0.05+right, xmax=0.116+0.05+right, ymin=ymin_) + 
    annotation_custom(gf, xmin=0.064+0.1+right, xmax=0.116+0.1+right, ymin=ymin_) #+

    #annotation_custom(gi, xmin=0.064+2*right+right2, xmax=0.116+2*right+right2, ymin=ymin_) + 
    #annotation_custom(gs, xmin=0.064+0.05+2*right+right2, xmax=0.116+0.05+2*right+right2, ymin=ymin_) + 
    #annotation_custom(gf, xmin=0.064+0.1+2*right+right2, xmax=0.116+0.1+2*right+right2, ymin=ymin_) 
  #ggExp   
  #ggsave('Outputs/Fig_3b_130mm_v1.png',ggExp, width = 13/(5/7), height =13, units = 'cm', bg="white", dpi = 600)
# mix (a) & (b) and export
  blank = ggplot() + theme_void() 
  ggB = ggarrange(blank, ggExp, blank, nrow=3, heights=c(10,92.9, 7.1+4.3))
  ggAll = ggarrange(ggA, ggB, ncol=2, widths=c(5,13))  
  ggsave('Outputs/Fig_2_width-180mm_v5.png',ggAll, width = 18/(5/7), height =16, units = 'cm', bg="white", dpi = 600)
  ggsave('Outputs/Fig_2_width-180mm_v5.jpg',ggAll, width = 18/(5/7), height =16, units = 'cm', bg="white", dpi = 600)
# END