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
  dr = d[bird_ID%in%d[duplicated(bird_ID), bird_ID]]

# Velocity
  R = rpt(VAP ~ (1 | bird_ID), grname = "bird_ID", data = dr, datatype = "Gaussian")
  RR = data.table(merge(data.frame(velocity ='VAP'), paste0(round(R$R*100),'%'))) %>% setnames(new = c('velocity', 'repeatability'))
  RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
  RR[, pred := 100*R$R]
  RR[, lwr := 100*R$CI_emp[1]]
  RR[, upr := 100*R$CI_emp[2]]
  VAP = RR

  R = rpt(VSL ~ (1 | bird_ID), grname = "bird_ID", data = dr, datatype = "Gaussian")
  RR = data.table(merge(data.frame(velocity ='VSL'), paste0(round(R$R*100),'%'))) %>% setnames(new = c('velocity', 'repeatability'))
  RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
  RR[, pred := 100*R$R]
  RR[, lwr := 100*R$CI_emp[1]]
  RR[, upr := 100*R$CI_emp[2]]
  VSL = RR

  R = rpt(VCL ~ (1 | bird_ID), grname = "bird_ID", data = dr, datatype = "Gaussian")
  RR = data.table(merge(data.frame(velocity ='VCL'), paste0(round(R$R*100),'%'))) %>% setnames(new = c('velocity', 'repeatability'))
  RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
  RR[, pred := 100*R$R]
  RR[, lwr := 100*R$CI_emp[1]]
  RR[, upr := 100*R$CI_emp[2]]
  VCL = RR

  r = rbind(VCL,VAP,VSL) 
  r[velocity == 'VCL', motility := 'Curvilinear']
  r[velocity == 'VAP', motility := 'Average path']
  r[velocity == 'VSL', motility := 'Straight line']

  fwrite(r[,.(motility, repeatability, CI)],file = 'Outputs/Table_SRa.csv')

  gv = 
    ggplot(r, aes(x = motility, y = pred)) +
        geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25) ) +
        #ggtitle ("Sim based")+
        geom_point(position = position_dodge(width = 0.25)) +
        #scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
        #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
        labs(x = NULL, y = "Repeatability [%]\n(within male)", subtitle = 'Motility')+
        ylim(c(0,100))+
        coord_flip()+
        theme_bw() +
        theme(plot.subtitle = element_text(size=9),
          legend.position = "none")

# Morphology
    lfrpt = list()
    for(i in c('Total','Flagellum','Head','Tail','Midpiece','Nucleus','Acrosome')){
      part_ = i
      # part_ = "Acrosome"
      bi = b[part == part_]
      R = rpt(Length_µm ~ (1 | bird_ID), grname = "bird_ID", data = bi, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
      lfrpt[[i]] = RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      print(i)
    }
    y = do.call(rbind,lfrpt)
    
    fwrite(y[order(nrow(y):1)], file = 'Outputs/Table_SRb.csv')

    y[, pred:= as.numeric(substr(Repeatability,1,2))]
    y[nchar(CI) == 5, CI := paste0(0,CI) ]
    y[, lwr:= as.numeric(substr(CI,1,2))]
    y[, upr:= as.numeric(substr(CI,4,5))]
    names(y)[2] = tolower( names(y)[2])
    y[, part := factor(part, levels=c('Total','Flagellum','Head','Tail','Midpiece','Nucleus','Acrosome'))] 
    
    g1 = 
      ggplot(y, aes(x = part, y = pred)) +
        geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25) ) +
        #ggtitle ("Sim based")+
        geom_point(position = position_dodge(width = 0.25)) +
        #scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
        #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
        labs(x = NULL, y = "Repeatability [%]\n(within male)", subtitle = 'Morphology')+
        ylim(c(0,100))+
        coord_flip()+
        theme_bw() +
        theme(plot.title = element_text(size=9),
          legend.position = "none")

    g2 = 
      ggplot(y, aes(x = part, y = pred)) +
        geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25) ) +
        #ggtitle ("Sim based")+
        geom_point(position = position_dodge(width = 0.25)) +
        #scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
        #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
        labs(x = NULL, y = "Repeatability [%]", subtitle = "within male")+
        ylim(c(0,100))+
        coord_flip()+
        theme_MB

    ggsave(file ='Outputs/Fig_R.png', g1, units = 'cm', width = 6, height = 6 )
# Combine 1
  gv = 
  ggplot(r, aes(x = motility, y = pred)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25) ) +
    #ggtitle ("Sim based")+
    geom_point(position = position_dodge(width = 0.25)) +
    #scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
    labs(x = NULL, y = "Repeatability [%]\n(within male)", subtitle = 'Motility')+
    coord_flip()+
    theme_bw() +
    theme(plot.subtitle = element_text(size=9),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(color = 'grey70'),
          plot.margin = margin(3,3,1,1, "mm"),
          legend.position = "none"
          )  

  gm = 
  ggplot(y, aes(x = part, y = pred)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25) ) +
    #ggtitle ("Sim based")+
    geom_point(position = position_dodge(width = 0.25)) +
    #scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
    labs(x = NULL, y = "Repeatability [%]\n(within male)", subtitle = 'Morphology')+
    coord_flip()+
    theme_bw() +
    theme(plot.subtitle = element_text(size=9),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 10),
          panel.border = element_rect(color = 'grey70'),
          plot.margin = margin(3,3,1,1, "mm"),
          legend.position = "none"
          )  

  ggA = ggarrange(
    gv,gm, 
    nrow=2, heights=c(3,7), align = 'v'
    )  
  ggA
  ggsave(here::here'Outputs/Fig_R.png',ggA, width = 6, height =9, units = 'cm', dpi = 600)
# Combine 2
    r[, part:=motility]
    r[, what := 'Motility']
    y[, what := 'Sperm length']
    ry=rbind(r[,.(what,part,pred,lwr,upr)], y[,.(what,part,pred,lwr,upr)])
    g =
    ggplot(ry, aes(y = part, x = pred)) +
          geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = 0.01) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.01)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          #scale_color_manual(values = col_,guide = guide_legend(reverse = TRUE))  +
          #scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          #scale_fill_manual(values = col_,guide = guide_legend(reverse = TRUE)) + 
          #scale_y_discrete(limits=rev, position = "right")+
          #coord_fixed(ratio = 0.05)+ #,xlim = c(-0.23, 0.15)
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          scale_x_continuous(limits = c(0, 100), expand = c(0, 0))+#, breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Repeatability [%]\n(within male) ") + #title = "a)",Effect of Period (before/during shutdown)
          #ylim(c(0,100))+
          #coord_flip()+
          facet_wrap(~what, nrow = 2, scales = 'free_y')+
          theme_MB +
          theme( legend.position ="none",
                plot.title = element_text(size=7),
                plot.tag = element_text(size=7),
                legend.title=element_text(size=7), 
                legend.text=element_text(size=6),
                ##legend.spacing.y = unit(0.1, 'cm'), 
                legend.key.height= unit(0.5,"line"),
                #plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                strip.background = element_blank(),
                strip.text = element_text(hjust = 0),
                axis.line = element_line(colour = ax_lines, size = 0.25),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x= element_line( colour = ax_lines, size = 0.25),
                axis.ticks.length = unit(1, "pt"),
                axis.text.x = element_text(colour="grey30", size = 6),
                axis.text.y=element_text(colour="grey30", size = 6),
                axis.title=element_text(size=7)
                )
    g
    ggsave(here::here('Outputs/Fig_R_v2.png'),g, width = 6, height =9, units = 'cm', dpi = 600)   

# END      