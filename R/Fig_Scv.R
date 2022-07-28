# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  require(ggpubr) 
  require(ggsci)
  require(gtools)
  
  require(admisc)
  require(gap)
  require(related)
  require(reshape2)

  width_ = 0.8 # spacing between error bars
  effects_ = c('Satellite relative to independent','Faeder relative to independent', 'Faeder relative to satellite')
# load data
  source(here::here('R/DAT_prepare.R'))      

  a[, animal := bird_ID]
  a[, id := bird_ID]
  b[, animal := bird_ID]
  b[, id := bird_ID]
  cv_[, animal := bird_ID]
  cv_[, id := bird_ID]
  # prepare motility
    dd = d[!Morph%in%'Zebra finch']
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

# load CV estimates from models controlled for relatedness (cauchy and default priors give same results)
    load(file = 'Outputs/CV_rel_control_cauchy_5000.Rdata') #load(file = 'Outputs/CV_rel_control_default_5000.Rdata')# #mi_, mi_co_
    s = mi_co_
    setnames(s, old = c('Estimate','CI.Lower','CI.Upper'), new = c('estimate','lwr','upr'))
    s[,model:='yes']
    s[, model := factor(model, levels=rev(c("no", "yes")))] 
    s[,response:=substring(response,4)]
    s[, response := factor(response, levels=rev(c("Acrosome", "Nucleus","Midpiece","Tail", "Head", "Flagellum","Total")))] 
    s[, effect := factor(effect, c("Faeder relative to satellite","Faeder relative to independent","Satellite relative to independent"))] 

# prepare CV estimates from simple models
    lcv = list()
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
      print(i)     
      }          
    llcv = data.table(do.call(rbind,lcv) ) 
    llcv[, response := factor(response, levels=rev(c("Acrosome", "Nucleus","Midpiece","Tail", "Head", "Flagellum","Total")))] 
    llcv[, effect := factor(effect, c("Faeder relative to satellite","Faeder relative to independent","Satellite relative to independent"))] 
    #llcv_ = llcv[response %in%c("Acrosome", "Nucleus","Midpiece","Tail","Total")]
    llcv[,model:='no']
    llcv[, model := factor(model, levels=rev(c("no", "yes")))] 
# combine
    cs = rbind(llcv,s[,.(response, effect, estimate, lwr, upr, model)])
    
# plot
  cols_=pal_jco()(3)  
  gCV = 
    ggplot(cs, aes(y = response, x = estimate, col = effect, fill = effect, shape = model)) +
    geom_vline(xintercept = 0, col = "grey60", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = width_) ) +
    geom_point(position = position_dodge(width =width_)) +
    scale_x_continuous(limits = c(-1.5, 2), expand = c(0, 0))+
    scale_color_jco(name = 'Contrast', guide = guide_legend(reverse = TRUE, order = 1,nrow=3,byrow=TRUE))+
    scale_fill_jco(name = 'Contrast', guide = guide_legend(reverse = TRUE, order = 1,nrow=3,byrow=TRUE))+
    scale_shape_manual(name = 'Model controlled for relatedness', values =c(23,21), guide = guide_legend(reverse = TRUE, override.aes = list(fill = c('grey30'), col = 'grey30'), order = 0,nrow=2,byrow=TRUE))+
    labs(y = NULL, x = "Standardized effect size", subtitle = 'Coefficient of variation')+
    theme_bw() +   
    theme(
        plot.subtitle = element_text(size=9, color = 'grey30'),
        legend.title=element_text(size=8, color = 'grey30'),
        legend.text=element_text(size=7.5, color = 'grey30'),
        #legend.spacing.y = unit(1, 'mm'),
        legend.key.height= unit(0.5,"line"),
        legend.margin=margin(0,0,0,-1),
        #legend.position=c(0.45,1.6),

        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 10, color ='grey10'),

        panel.border = element_rect(color = 'grey70'),
        panel.grid.minor = element_blank(),

        plot.margin = margin(3,3,1,1, "mm")
        )  
  ggsave('Outputs/Fig_Scv_90mm.png',gCV, width = 9/(5/7), height =8/(5/7), units = 'cm', bg="white", dpi = 600)
  #ggsave('Outputs/Fig_Scv_default_90mm.png',gCV, width = 9/(5/7), height =8/(5/7), units = 'cm', bg="white", dpi = 600)

# end    