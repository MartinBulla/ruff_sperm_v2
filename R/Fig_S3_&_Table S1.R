# =============================================================
# ❗ Runs relative to the project's root directory,
# and exports Fig S3 & Table S1 into ./Outputs/
# =============================================================

# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(ggpubr) 
  require(rptR) 

  # constants
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    fae = '#d4b691' # 'ffd6af'
    sat = 'white'
    ind = '#303030'
    colors = c(ind,sat,fae)
  
  # functions
     # for repeatability output based on sim
        R_out = function(name = "define", model = m, nsim = 5000){
         bsim <- sim(model, n.sim=nsim)  
         l=data.frame(summary(model)$varcor)
         l = l[is.na(l$var2),]
         l$var1 = ifelse(is.na(l$var1),"",l$var1)
         l$pred = paste(l$grp,l$var1)

         q50={}
         q025={}
         q975={}
         pred={}
         
         # variance of random effects
         for (ran in names(bsim@ranef)) {
           #ran =names(bsim@ranef)[1]
           ran_type = l$var1[l$grp == ran]
           for(i in ran_type){
              # i = ran_type[2]
            q50=c(q50,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.5)))
            q025=c(q025,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.025)))
            q975=c(q975,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.975)))
            pred= c(pred,paste(ran, i))
            }
           }
         # residual variance
         q50=c(q50,quantile(bsim@sigma^2, prob=c(0.5)))
         q025=c(q025,quantile(bsim@sigma^2, prob=c(0.025)))
         q975=c(q975,quantile(bsim@sigma^2, prob=c(0.975)))
         pred= c(pred,'Residual')

         ci = c(round(100*q025/sum(q025))[1], round(100*q975/sum(q975))[1])
         ci = ci[order(ci)]
         
         ri=data.table(model = name, repeatability=paste0(round(100*q50/sum(q50)),'%')[1], CI = paste0(paste(ci[1], ci[2], sep ="-"), '%'))
         
         
         return(ri)
         }
# DATA
  source(here::here('R/DAT_prepare.R'))       

# estimate Velocity
  dr = d[bird_ID%in%d[duplicated(bird_ID), bird_ID]]
  drw = reshape(dr[,.(month,bird_ID,VAP,VSL,VCL, motileCount, Morph, age)], idvar = c('bird_ID','Morph','age'), timevar = 'month', direction = "wide") 
  R = rpt(VAP ~ (1 | bird_ID), grname = "bird_ID", data = dr, datatype = "Gaussian")
  RR = data.table(merge(data.frame(trait = "Velocity", velocity = "Average path"), paste0(round(R$R * 100), "%"))) %>% setnames(new = c("trait", "specification", "repeatability"))
  RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
  RR[, pred := 100*R$R]
  RR[, lwr := 100*R$CI_emp[1]]
  RR[, upr := 100*R$CI_emp[2]]
  VAP = RR

  R = rpt(VSL ~ (1 | bird_ID), grname = "bird_ID", data = dr, datatype = "Gaussian")
  RR = data.table(merge(data.frame(trait = "Velocity", velocity = "Straight line"), paste0(round(R$R * 100), "%"))) %>% setnames(new = c("trait", "specification", "repeatability"))
  RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
  RR[, pred := 100*R$R]
  RR[, lwr := 100*R$CI_emp[1]]
  RR[, upr := 100*R$CI_emp[2]]
  VSL = RR

  R = rpt(VCL ~ (1 | bird_ID), grname = "bird_ID", data = dr, datatype = "Gaussian")
  RR = data.table(merge(data.frame(trait = "Velocity", velocity = "Curvilinear"), paste0(round(R$R * 100), "%"))) %>% setnames(new = c("trait", "specification", "repeatability"))
  RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
  RR[, pred := 100*R$R]
  RR[, lwr := 100*R$CI_emp[1]]
  RR[, upr := 100*R$CI_emp[2]]
  VCL = RR

  r = rbind(VCL, VSL, VAP)
  
# estimate Morphology
    lfrpt = list()
    for(i in c('Total','Flagellum','Head','Tail','Midpiece','Nucleus','Acrosome')){
      part_ = i
      # part_ = "Acrosome"
      bi = b[part == part_]
      R = rpt(Length_µm ~ (1 | bird_ID), grname = "bird_ID", data = bi, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(trait = "Length", name = part_), paste0(round(R$R * 100), "%"))) %>% setnames(new = c("trait", "specification", "repeatability"))
      lfrpt[[i]] = RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      print(i)
    }
    y = do.call(rbind,lfrpt)
    y[, pred := as.numeric(substr(repeatability, 1, 2))]
    y[, lwr := as.numeric(substr(CI, 1, 2))]
    y[, upr := as.numeric(substr(CI, 4, 5))]

# Table S1
  fwrite(rbind(r, y[order(nrow(y):1)])[, 1:4], file = "Outputs/Table_S1.csv")

# Figure S3
  r[, specification := factor(specification, levels = rev(c("Curvilinear", "Straight line", "Average path")))]
  gv = 
  ggplot(r, aes(x = specification, y = pred)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25) ) +
    #ggtitle ("Sim based")+
    geom_point(position = position_dodge(width = 0.25)) +
    #scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0,100, by = 20), expand = c(0, 0))+
    labs(x = NULL, y = "Repeatability [%]\n(within male)", subtitle = 'Velocity')+
    coord_flip()+
    theme_bw() +
    theme(plot.subtitle = element_text(size=9, color = 'grey30'),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(color = 'grey70'),
          plot.margin = margin(3,3,1,1, "mm"),
          legend.position = "none"
          )  
  
  y[, specification := factor(specification, levels = rev(c("Acrosome", "Nucleus", "Midpiece", "Tail", "Total", "Head", "Flagellum")))]
  gm = 
  ggplot(y, aes(x = specification, y = pred)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25) ) +
    #ggtitle ("Sim based")+
    geom_point(position = position_dodge(width = 0.25)) +
    #scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0,100, by = 20),expand = c(0, 0))+
    labs(x = NULL, y = "Repeatability [%]\n(within male)", subtitle = 'Morphology')+
    #annotate(geom = "text", y = 85, x = 2, label = "Composite\ntraits", color = "grey50",
    #         angle = 90, size  = 3)+
    coord_flip()+
    theme_bw() +
    theme(plot.subtitle = element_text(size=9, color = 'grey30'),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 10, color ='grey10'),
          panel.border = element_rect(color = 'grey70'),
          plot.margin = margin(3,3,1,1, "mm"),
          legend.position = "none"
          )  

  ggA = ggarrange(
    gv,gm, 
    nrow=2, heights=c(3,7), align = 'v'
    )  
  ggA
  ggsave(here::here('Outputs/Fig_S3_width-43mm.png'),ggA, width = 4.3/(5/7), height =9, units = 'cm', dpi = 600)

# END