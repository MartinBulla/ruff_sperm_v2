#' ---
#' title: "Lack of major differences in sperm traits of highly divergent male reproductive morphs in the ruff"
#' author: "Martin Bulla"
#' date: "`r Sys.time()`"
#' output: 
#'     html_document:
#'         toc: true
#'         toc_float: true
#'         toc_depth: 5
#'         code_folding: hide
#'         bibliography: ruff_sperm.xml
#'         link-citations: yes
#' ---

#+ r setup, include=FALSE 
knitr::opts_chunk$set(message = FALSE, warning = FALSE, cache = TRUE)

#' # Code to load tools and data
  # figs checked for colorblind friendliness with https://www.color-blindness.com/coblis-color-blindness-simulator/
  # TOOLS 
    require(here)
    source(here::here('R/tools.R'))
    
    require(arm) 
    require(effects)
    require(ggpubr)
    require(ggsci) 
    require(grid)
    require(gridExtra)
    require(magrittr)
    require(MASS)
    require(multcomp)
    require(PerformanceAnalytics)
    require(readxl)
    require(rptR) 
    require(stringi)
    require(viridis)


    # constants
      round_ = 3 # number of decimal places to round model coefficients
      nsim = 5000 # number of simulations to extract estimates and 95%CrI
      ax_lines = "grey60" # defines color of the axis lines
      colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)

      fae = '#d4b691' # 'ffd6af'
      sat = 'white'
      ind = '#303030'
    
    # functions
      getime = function (x) {ifelse(is.na(x), as.numeric(NA), as.numeric(difftime(x, trunc(x,"day"), units = "hours")))}
      
      getDay = function (x) {as.Date(trunc(x, "day"))}
     
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
      # custom ggplot theme
         theme_MB = theme(  
                  title = element_text(size=8, colour="grey30"),
                  axis.line = element_blank(),
                  #axis.line = element_line(colour="grey70", size=0.25),
                  axis.title = element_text(size=7, colour="grey30"),
                  axis.title.y = element_text(vjust=3.5),
                  axis.title.x = element_text(vjust=1),
                  axis.text = element_text(size=6),#, vjust = 0.5, hjust=1),# margin=units(0.5,"mm")),
                  axis.ticks.length=unit(0.5,"mm"),
                  axis.ticks = element_line(colour = "grey70", size = 0.1),
                  #axis.ticks.margin,
                  
                  strip.text.x = element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
                  strip.text.y = element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
                  strip.background = element_rect(fill="grey99",colour="grey70", size=0.25),
                    #strip.background = element_blank(), 
                    #strip.text = element_blank(),
                  panel.spacing = unit(0, "mm"),
                  panel.background=element_blank(),
                  panel.border = element_rect(colour="grey70", size=0.1, fill = NA), #panel.border=element_blank(),
                  panel.grid = element_blank(),

                  legend.text=element_text(size=6),
                  legend.title=element_text(size=6),
                  legend.key = element_rect(colour = NA, fill = NA),
                  legend.key.height= unit(0.5,"line"),
                  legend.key.width = unit(0.25, "cm"),
                  legend.margin = margin(0,0,0,0, unit="cm"),
                  legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
                  legend.background = element_blank()
                  )  
  # DATA 
    source(here::here('R/DAT_prepare.R'))       
   # prepare motility for plotting
      am = a[part=='Acrosome']
      aml = data.table(melt(am[,.(bird_ID,month,Morph,age,HL,motileCount,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph","age","HL","motileCount"), variable.name = "Motility"))
      aml[Motility == 'VAP' ,mot:='Average path']
      aml[Motility == 'VCL' ,mot:='Curvilinear']
      aml[Motility == 'VSL' ,mot:='Straight line']

# Fig 1
  t = data.table(read_excel(here::here('Data/testes.xlsx'), sheet = 1))#, range = "A1:G161"))
  t[, Morph := factor(Morph, levels=c("Res", "Sat", "Faed"))] 
  t[Morph == 'Res', Morph := 'Independent']
  t[Morph == 'Sat', Morph := 'Satellite']
  t[Morph == 'Faed', Morph := 'Faeder']
  t[, soma := Bodymass - Gonadmass]

  # not used - images - within the plot
   g1 = ggplot(t, aes(x = Morph, y = Gonadmass)) + 
    annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=4.5, ymax=5.1) + 
    annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=4.5, ymax=5.1) + 
    annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=4.47, ymax=4.9) +
    geom_boxplot() + 
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'darkgrey', aes(fill =Morph), dotsize = 1.1)+
    scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Testes mass [g]', limits = c(1.96,5.1), breaks = c(2,2.5,3,3.5,4,4.5))+
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
      #axis.title.y = element_text(size = 8),
      legend.position = "none")

   g2 = ggplot(t, aes(x = Morph, y = Bodymass)) + 
    annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=139, ymax=160) + 
    annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=139, ymax=160) + 
    annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=140, ymax=151) + 
    geom_boxplot() +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'darkgrey', aes(fill =Morph))+
    
    scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Bodymass [g]', limits = c(100,200), expand = c(0, 0))+
    theme_bw()+theme(legend.position = "none"
      #axis.title.y = element_text(size = 8)
      )                  
  
   gg1 <- ggplotGrob(g1)
   gg2 <- ggplotGrob(g2) 
   grid.draw(rbind(gg1, gg2))

   ggsave('Outputs/Fig_1-testes-body_01.png',rbind(gg1,gg2, size = "last"), width = 7, height =10, units = 'cm')  

  # not used - images - ontop of the plot
   g0 = ggplot(t, aes(x = Morph, y = Gonadmass)) + 
    geom_boxplot() + 
    annotation_custom(gi, xmin=0.75, xmax=1.25, ymin = 4.665) + 
    annotation_custom(gs, xmin=1.75, xmax=2.25, ymin = 4.665)+#, ymin=4.5, ymax=5.1) + 
    annotation_custom(gfc, xmin=2.77, xmax=3.23,ymin=4.655)+#, ymin=4.68, ymax=4.82) +
    scale_y_continuous('Testes mass [g]', limits = c(4.7,4.8), breaks = c(4.7,4.8))+#scale_y_continuous(limits = c(4.9,5.1))+
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_text(size = 10, color = 'transparent'), 
          axis.text.y = element_text(color = 'transparent'),
          axis.ticks = element_blank(),
          plot.margin = unit(c(1,1.75,0,2), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background  = element_rect(fill='white', color = 'white'),
          legend.position = "none"
          )
   
   g1 = ggplot(t, aes(x = Morph, y = Gonadmass)) + 
    #annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=4.5, ymax=5.1) + 
    #annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=4.5, ymax=5.1) + 
    #annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=4.47, ymax=4.9) +
    geom_boxplot(col = 'grey50') + 
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'darkgrey', aes(fill =Morph), dotsize = 1.1)+
    scale_color_manual(values=c(ind,sat,fae))+
    scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Testes mass [g]', limits = c(1.93,4.5),expand = c(0, 0))+
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
      axis.title = element_text(size = 10),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = 'grey70'),
      #axis.title.y = element_text(size = 8),
      legend.position = "none")

   g2 = ggplot(t, aes(x = Morph, y = Bodymass)) + 
    #annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=139, ymax=160) + 
    #annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=139, ymax=160) + 
    #annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=140, ymax=151) + 
    geom_boxplot(col = 'grey50') +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(),  aes(fill =Morph, col = Morph))+
    
    scale_color_manual(values=c('black','darkgrey','#bf925a'))+ #scale_fill_viridis(discrete=TRUE)+
    scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Body mass [g]', limits = c(100,200), expand = c(0, 0))+
    theme_bw()+
    theme(axis.ticks = element_blank(),
      axis.title = element_text(size = 10),
      panel.border = element_rect(color = 'grey70'),
      legend.position = "none"
      #axis.title.y = element_text(size = 8)
      )     

   ggA = ggarrange(
    g0,#+theme(plot.margin = unit(c(0,1.75,0,2), "mm")),
    g1+theme(plot.margin = unit(c(0,1.75,0.3,2), "mm")),
    g2+theme(plot.margin = unit(c(1.5,1.75,0.3,2), "mm")),  
    nrow=3, heights=c(1.5, 4, 4.8),  align = 'v'
    )  
   ggA         
   ggsave('Outputs/Fig_1-testes-body_width-55mm_v5.png',ggA, width = 7, height =13, units = 'cm', bg="white", dpi = 600)
  # images - ontop of the plot, larger icons
   g0 = ggplot(t, aes(x = Morph, y = Gonadmass)) + 
    geom_boxplot() + 
    annotation_custom(gi, xmin=0.725, xmax=1.3, ymin = 4.665) + 
    annotation_custom(gs, xmin=1.725, xmax=2.3, ymin = 4.665)+#, ymin=4.5, ymax=5.1) + 
    annotation_custom(gfc, xmin=2.745, xmax=3.255,ymin=4.655)+#, ymin=4.68, ymax=4.82) +
    scale_y_continuous('Testes mass [g]', limits = c(4.7,4.8), breaks = c(4.7,4.8))+#scale_y_continuous(limits = c(4.9,5.1))+
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_text(size = 10, color = 'transparent'), 
          axis.text.y = element_text(color = 'transparent'),
          axis.ticks = element_blank(),
          plot.margin = unit(c(1,1.75,0,2), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background  = element_rect(fill='white', color = 'white'),
          legend.position = "none"
          )
   
   g1 = ggplot(t, aes(x = Morph, y = Gonadmass)) + 
    #annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=4.5, ymax=5.1) + 
    #annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=4.5, ymax=5.1) + 
    #annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=4.47, ymax=4.9) +
    geom_boxplot(col = 'grey50', outlier.shape = NA) + 
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(),  aes(col = Morph, fill =Morph), dotsize = 1.1)+ #col = 'darkgrey',
    scale_color_manual(values=c('black','darkgrey','#bf925a'))+ #scale_fill_viridis(discrete=TRUE)+
    scale_fill_manual(values=c(ind,sat,fae))+
    #scale_color_manual(values=c(ind,sat,fae))+
    #scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Testes mass [g]', limits = c(1.93,4.5),expand = c(0, 0))+
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
      axis.title = element_text(size = 10, , colour="grey10"),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = 'grey70'),
      #axis.title.y = element_text(size = 8),
      legend.position = "none")

   g2 = ggplot(t, aes(x = Morph, y = Bodymass)) + 
    #annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=139, ymax=160) + 
    #annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=139, ymax=160) + 
    #annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=140, ymax=151) + 
    geom_boxplot(col = 'grey50', outlier.shape = NA) +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(),  aes(fill =Morph, col = Morph))+
    
    scale_color_manual(values=c('black','darkgrey','#bf925a'))+ #scale_fill_viridis(discrete=TRUE)+
    scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Body mass [g]', limits = c(100,200), expand = c(0, 0))+
    theme_bw()+
    theme(axis.ticks = element_blank(),
      axis.title = element_text(size = 10, colour="grey10"),
      panel.border = element_rect(color = 'grey70'),
      legend.position = "none"
      #axis.title.y = element_text(size = 8)
      )     

   ggA = ggarrange(
    g0,#+theme(plot.margin = unit(c(0,1.75,0,2), "mm")),
    g1+theme(plot.margin = unit(c(0,1.75,0.3,2), "mm")),
    g2+theme(plot.margin = unit(c(1.5,1.75,0.3,2), "mm")),  
    nrow=3, heights=c(1.5, 4, 4.8),  align = 'v'
    )  
   ggA         
   ggsave('Outputs/Fig_1_width-50mm.png',ggA, width = 7, height =13, units = 'cm', bg="white", dpi = 600)

  # not used - images - ontop of the plot clean & 5cm
   g0 = ggplot(t, aes(x = Morph, y = round(Gonadmass,1))) + 
    geom_boxplot() + 
    annotation_custom(gi, xmin=0.75, xmax=1.25, ymin = 4.665) + 
    annotation_custom(gs, xmin=1.75, xmax=2.25, ymin = 4.665)+#, ymin=4.5, ymax=5.1) + 
    annotation_custom(gfc, xmin=2.77, xmax=3.23,ymin=4.655)+#, ymin=4.68, ymax=4.82) +
    scale_y_continuous('Testes mass [g]', limits = c(4.7,4.8), breaks = c(4.7,4.8))+#scale_y_continuous(limits = c(4.9,5.1))+
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_text(size=7, color = 'transparent'), 
          axis.text.y = element_text(size=6, color = 'transparent'), 
          axis.ticks = element_blank(),
          plot.margin = unit(c(1,1.75,0,2), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background  = element_rect(fill='white', color = 'white'),
          legend.position = "none"
          )
   
   g1 = ggplot(t, aes(x = Morph, y = Gonadmass)) + 
    #annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=4.5, ymax=5.1) + 
    #annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=4.5, ymax=5.1) + 
    #annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=4.47, ymax=4.9) +
    geom_boxplot(col = 'grey50', outlier.shape = NA) + 
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'darkgrey', aes(fill =Morph), dotsize = 1.1)+
    scale_color_manual(values=c(ind,sat,fae))+
    scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Testes mass [g]', limits = c(1.93,4.5),expand = c(0, 0))+
    theme_MB +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
      # axis.title.y = element_text(colour="grey10"),
      # axis.title = element_text(size = 10),
      axis.ticks = element_blank(),
      #panel.border = element_rect(color = 'grey70'),
      #axis.title.y = element_text(size = 8),
      legend.position = "none")

   g2 = ggplot(t, aes(x = Morph, y = Bodymass)) + 
    #annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=139, ymax=160) + 
    #annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=139, ymax=160) + 
    #annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=140, ymax=151) + 
    geom_boxplot(col = 'grey50',outlier.shape = NA) +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(),  aes(fill =Morph, col = Morph))+
    scale_color_manual(values=c('black','darkgrey','#bf925a'))+ #scale_fill_viridis(discrete=TRUE)+
    scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Body mass [g]', limits = c(100,200), expand = c(0, 0))+
    theme_MB+
    theme(axis.ticks = element_blank(),
      #axis.title = element_text(colour="grey10"),
      #axis.title = element_text(size = 10),
     # panel.border = element_rect(color = 'grey70'),
      legend.position = "none"
      #axis.title.y = element_text(size = 8)
      )     

   ggA = ggarrange(
    g0,#+theme(plot.margin = unit(c(0,1.75,0,2), "mm")),
    g1+theme(plot.margin = unit(c(0,1.75,0.3,2), "mm")),
    g2+theme(plot.margin = unit(c(1.5,1.75,0.3,2), "mm")),  
    nrow=3, heights=c(1.5, 4, 4.8),  align = 'v'
    )  

   ggA         
   ggsave('Outputs/Fig_1-testes-body_v10-small.png',ggA, width = 5, height =8.57, units = 'cm', bg="white", dpi = 600)
  # not used - images - ontop of the plot clean & 5cm & larger icons
   g0 = ggplot(t, aes(x = Morph, y = Gonadmass)) + 
    geom_boxplot() +
    annotation_custom(gi, xmin=0.725, xmax=1.3, ymin = 4.665) + 
    annotation_custom(gs, xmin=1.725, xmax=2.3, ymin = 4.665)+#, ymin=4.5, ymax=5.1) + 
    annotation_custom(gfc, xmin=2.7405, xmax=3.282,ymin=4.655)+#, ymin=4.68, ymax=4.82) +
    scale_y_continuous('Testes mass [g]', limits = c(4.7,4.8), breaks = c(4.7,4.8))+#scale_y_continuous(limits = c(4.9,5.1))+
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_text(size=7, color = 'transparent'), 
          axis.text.y = element_text(size=6, color = 'transparent'), 
          axis.ticks = element_blank(),
          plot.margin = unit(c(1,1.75,0,2), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background  = element_rect(fill='white', color = 'white'),
          legend.position = "none"
          )
   
   g1 = ggplot(t, aes(x = Morph, y = Gonadmass)) + 
    #annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=4.5, ymax=5.1) + 
    #annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=4.5, ymax=5.1) + 
    #annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=4.47, ymax=4.9) +
    geom_boxplot(col = 'grey50', outlier.shape = NA) + 
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'darkgrey', aes(fill =Morph), dotsize = 1.1)+
    scale_color_manual(values=c(ind,sat,fae))+
    scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Testes mass [g]', limits = c(1.93,4.5),expand = c(0, 0))+
    theme_MB +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
      # axis.title.y = element_text(colour="grey10"),
      # axis.title = element_text(size = 10),
      axis.ticks = element_blank(),
      #panel.border = element_rect(color = 'grey70'),
      #axis.title.y = element_text(size = 8),
      legend.position = "none")

   g2 = ggplot(t, aes(x = Morph, y = Bodymass)) + 
    #annotation_custom(gi, xmin=0.75, xmax=1.25, ymin=139, ymax=160) + 
    #annotation_custom(gs, xmin=1.75, xmax=2.25, ymin=139, ymax=160) + 
    #annotation_custom(gf, xmin=2.75, xmax=3.25, ymin=140, ymax=151) + 
    geom_boxplot(col = 'grey50',outlier.shape = NA) +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(),  aes(fill =Morph, col = Morph))+
    scale_color_manual(values=c('black','darkgrey','#bf925a'))+ #scale_fill_viridis(discrete=TRUE)+
    scale_fill_manual(values=c(ind,sat,fae))+ #scale_fill_viridis(discrete=TRUE)+
    scale_y_continuous('Body mass [g]', limits = c(100,200), expand = c(0, 0))+
    theme_MB+
    theme(axis.ticks = element_blank(),
      #axis.title = element_text(colour="grey10"),
      #axis.title = element_text(size = 10),
     # panel.border = element_rect(color = 'grey70'),
      legend.position = "none"
      #axis.title.y = element_text(size = 8)
      )     

   ggA = ggarrange(
    g0,#+theme(plot.margin = unit(c(0,1.75,0,2), "mm")),
    g1+theme(plot.margin = unit(c(0,1.75,0.3,2), "mm")),
    g2+theme(plot.margin = unit(c(1.5,1.75,0.3,2), "mm")),  
    nrow=3, heights=c(1.5, 3.95, 4.85),  align = 'v'
    )  

   ggA         
   ggsave('Outputs/Fig_1-testes-body_50mm_clean_v02.png',ggA, width = 5, height =9, units = 'cm', bg="white", dpi = 600)

  #gg0 <- ggplotGrob(g0)
  #gg1 <- ggplotGrob(g1)
  #gg2 <- ggplotGrob(g2) 
  
  #g12 <- rbind(gg1,gg2, size = "last") 
  #grid.draw(rbind(gg0, g12))
**Figure 1** | Morph-specific testes size.
   
# METHODS
## distribution of correlation coefficients for traits with inbreeding
    r1 = aml[, cor(HL,value), by = mot]
    names(r1)[1] = 'trait'
    r2 = a[, cor(HL,Length_avg), by = part]
    names(r2)[1] = 'trait'
    r3= cv_[, cor(HL,CV), by = part]
    names(r3)[1] = 'trait'
    summary(rbind(r1,r2,r3))

## Fig Shl
  # dummie dataset to create empty panels
    dum = aml[1:4] 
    dum[,mot2 := c('w','x','y','z')]
    aml[, mot2 := factor(mot,levels=c('Curvilinear','Straight line','Average path','w','x','y','z'))]
    amld = rbind(aml,dum)

  gv =
  ggplot(aml, aes(x = HL, y = value)) +
    stat_smooth(method = MASS::rlm, col ='grey30') +
    geom_point(pch = 21, alpha = 0.75, aes(fill = Morph, col = Morph))+
    #stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    facet_wrap(~mot2, scales = 'free_y', nrow = 1,drop=FALSE)+
  
    scale_color_manual(values = cols)+
    scale_fill_manual(values = fills)+
    
    xlab('Homozygousity by locus') +
    ylab('Motility [μm/s]') +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size=8),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 8,colour="grey10"),
      axis.title.x = element_blank(),
      axis.text = element_text(size=7), 
      axis.text.x = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_text(size = 7, color="grey20",  margin=margin(1,1,1,1,"mm")),
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
      
      panel.border = element_rect(color = 'grey70')
      )  

  ggv = ggplotGrob(gv)
  rm_grobs <- ggv$layout$name %in% c("panel-4-1", "panel-5-1","panel-6-1","panel-7-1", "strip-t-4-1", "strip-t-5-1","strip-t-6-1", "strip-t-7-1")

  ggv$grobs[rm_grobs] <- NULL
  ggv$layout <- ggv$layout[!rm_grobs, ]

  gm =
  ggplot(a, aes(x = HL, y = Length_avg)) +
    stat_smooth(method = MASS::rlm, col ='grey30') +
    geom_point(pch = 21, alpha = 0.75, aes(fill = Morph, col = Morph))+
    #stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    facet_wrap(~part, scales = 'free_y', nrow = 1,drop=FALSE)+
  
    scale_color_manual(values = cols)+
    scale_fill_manual(values = fills)+
    
    xlab('Homozygousity by locus') +
    ylab('Length [µm]') +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size=8),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 8,colour="grey10"),
      axis.title.x = element_blank(),
      axis.text = element_text(size=7), 
      axis.text.x = element_blank(),
      axis.text.y = element_text(margin = margin(r = -1)),
      
      strip.text = element_text(size = 7, color="grey20",  margin=margin(1,1,1,1,"mm")),
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
      
      panel.border = element_rect(color = 'grey70')
      )  
  
  gcv = 
  ggplot(cv_, aes(x = HL, y = CV)) +
    stat_smooth(method = MASS::rlm, col ='grey30') +
    geom_point(pch = 21, alpha = 0.75, aes(fill = Morph, col = Morph))+
    #stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
    stat_cor(method="pearson",size = 2, cor.coef.name = 'r',aes(label = ..r.label..)) +
    facet_wrap(~part, scales = 'free_y', nrow = 1,drop=FALSE)+
  
    scale_color_manual(values = cols)+
    scale_fill_manual(values = fills)+
    
    xlab('Homozygousity by locus') +
    ylab('Coefficient of variation') +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size=8),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 8,colour="grey10"),
      axis.text = element_text(size=7), 
      axis.text.x = element_text(margin = margin(t = -1)),
      axis.text.y = element_text(margin = margin(r = -1)),
      strip.text = element_blank(),
      strip.background = element_rect(fill=NA,colour=NA, size=0.25),
      
      panel.border = element_rect(color = 'grey70')
      )  

  ggHL = ggarrange(
    ggv,
    gm,
    gcv,  
    nrow=3, align = 'v' #heights=c(1.5, 4, 4.8),  
    )   
  
   # prepare legend
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

  g_exp =   
    ggHL + 
    annotation_custom(gp_ind_grob, xmin=.47, xmax=.52, ymin = 0.5)+
    annotation_custom(gi, xmin=.47, xmax=.52, ymin = 0.7) +
    annotation_custom(gp_sat_grob, xmin=.47+.06, xmax=.52+.06, ymin = 0.5) +
    annotation_custom(gs, xmin=.47+.06, xmax=.52+.06, ymin = 0.7) +
    annotation_custom(gp_fae_grob, xmin=.47+.12, xmax=.52+.12, ymin = 0.5) +
    annotation_custom(gf, xmin=.47+.12, xmax=.52+.12, ymin = 0.7) 

  ggsave('Outputs/Fig_Shl_140mm.png',g_exp, width = 14/(5/7), height =9, units = 'cm', bg="white", dpi = 600)


# RESULTS

bw[!duplicated(bird_ID),summary(Morph)] # N sampled males
nrow(ss[duplicated(bird_ID)]) # for 42 males velocity recorded twice - in May and June 
ss[bird_ID%in%ss[duplicated(bird_ID), bird_ID], summary(as.factor(month))]

# cor among motility estimates
  

# Fig Sm
  dl = data.table(melt(d[,.(bird_ID,month,Morph,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph"), variable.name = "mot"))
  dl[mot =='VAP', mot:='Average path']
  dl[mot =='VSL', mot:='Straight line']
  dl[mot =='VCL', mot:='Curvilinear']
  dl[, mot := factor(mot,levels=c('Curvilinear','Straight line','Average path'))]

  g = ggplot(dl, aes(x = Morph, y = value)) + 
      facet_wrap(~mot, scales = 'free_y', nrow = 3, strip.position="right") +
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(),  aes(fill =Morph, col = Morph))+
      geom_boxplot(col = 'grey50', fill = NA, outlier.shape = NA) +
      scale_color_manual(values=c(cols, '#a53708'))+
      scale_fill_manual(values=c(ind,sat,fae,'#f89f79'))+
      scale_y_continuous('Velocity [μm/s]', expand = c(0, 0))+
      coord_cartesian(clip = 'off') +
      theme_bw()+
      theme(
        axis.ticks = element_blank(),
        axis.title = element_text(size = 10, colour="grey10"),
        
        strip.text.y.right = element_text(color="grey20",  margin=margin(1,1,1,1,"mm"), angle=90),
        strip.background = element_rect(fill=NA,colour=NA, size=0.25),
        
        panel.border = element_rect(color = 'grey70'),
        plot.margin = margin(14,3,1,1, "mm"),

        legend.position = "none"
        #axis.title.y = element_text(size = 8)
        )  
  gg = ggarrange(g)   
  ggExp = gg + 
        annotation_custom(gi, xmin=0.125, xmax=0.275, ymin=0.91) + 
        annotation_custom(gs, xmin=0.125+0.2, xmax=0.275+0.2, ymin=0.91) + 
        annotation_custom(gf, xmin=0.125+0.4, xmax=0.275+0.4, ymin=0.906) +
        annotation_custom(gz, xmin=0.125+0.6, xmax=0.275+0.6, ymin=0.88) 
  ggsave('Outputs/Fig_Sm_zebra_v2_60mm.png',ggExp, width = 6/(5/7), height =13, units = 'cm', bg="white", dpi = 600)


# discriminant analyses
  m =  lda(Morph~scale(Nucleus)+scale(Midpiece)+scale(Tail) + scale(VSL), data = aw)    
  m =  lda(Morph~scale(Nucleus)+scale(Midpiece)+scale(Tail), data = aw)    
  m =  lda(Morph~scale(Nucleus)+scale(Midpiece), data = aw)    
  m =  lda(Morph~scale(Midpiece), data = aw)    
  m
  plot(m)
  predictions <- m %>% predict(aw)
  mean(predictions$class==aw$Morp)
  
  names(predictions)
  head(predictions$class, 6)
  head(predictions$posterior, 6) 
  head(predictions$x, 3) 

  lda.data <- cbind(aw, predict(m)$x)
  ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(color = Morph))

  mean(predictions$class==aw$Morp)

  ggplot(aw,aes(x = Midpiece, y = Tail, col = Morph)) + geom_point()

  as = a[, mean(Length_avg), by = list(part,Morph)]
    names(as)[3] = 'median'
    as$lwr = a[, mean(Length_avg)-sd(Length_avg), by = list(part,Morph)]$V1
    as$upr = a[, mean(Length_avg)+sd(Length_avg), by = list(part,Morph)]$V1

    as_n = as[part == 'Nucleus']
    as_n[,Midpiece_µm:= as[part == 'Midpiece',.(median)]]

    as_nm = as[part == 'Midpiece']
    as_nm[,Nucleus_µm:= as[part == 'Nucleus',.(median)]]

  g0a =
    ggplot(ha, aes(x = Nucleus_µm, y = Midpiece_µm, col = Morph, fill = Morph)) +
      geom_point(pch = 21)+
      geom_point(data = as_n, aes(x = median, y = Midpiece_µm), col = 'black', size = 2) +
      geom_segment(data = as_n, aes(x = lwr, xend = upr, y =Midpiece_µm, yend =Midpiece_µm), col = "black" ) +
      geom_segment(data = as_nm, aes(x = Nucleus_µm, xend = Nucleus_µm, y =lwr, yend =upr), col = "black" )+
      geom_point(data = as_n, aes(x = median, y = Midpiece_µm, col = Morph), size = 1) +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = aas_f, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
     #scale_colour_viridis(discrete=TRUE)+
      ggtitle('Male means') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

```{r Predictions,  warning = FALSE, message = FALSE, fig.align="center", fig.width=3, fig.height=3}
  set.seed(1)
  d = data.table(morph = factor(c( rep('independent', 60),
                                   rep('satelite', 30),
                                   rep('faeder', 15)
                                   ), 
                                levels = c('independent','satelite','faeder')
                                ),
                sperm_trait =  c(rnorm(60, mean = 1, sd = 0.5), 
                                 rnorm(30, mean = 2, sd = 0.5),
                                 rnorm(15, mean = 3, sd = 0.5) 
                                 ) 
                )
  #dev.new(width = 3.5, height = 3.5)
  ggplot(d, aes(x = morph, y = sperm_trait, col = morph)) +  
    geom_jitter() + theme_MB +  theme(legend.position="none") + xlab('Morph') + ylab('Sperm length or relative midpiece length') +
    scale_colour_manual(values = colors)
 
```
**Figure 2** | Predicted relationship between morph type and sperm traits.


## 2. Does sperm morphology predict sperm swimming speed?

Sperm morphology, specifically long sperm, sperm with short heads, and long midpiece (relative to overall sperm length), has been positively associated with sperm swimming speed [**CIT**]. We thus predict that longer sperm and/or sperm with longer head, and/or sperm with relatively longer midpiece will swim faster.

```{r Predictions2,  warning = FALSE, message = FALSE, fig.align="center", fig.width=4, fig.height=4}
  set.seed(1)
  d = data.table(morph = c( rep('independent', 60),
                                   rep('satelite', 30),
                                   rep('faeder', 15)
                                   ),
                midpiece =  c(rnorm(60, mean = 2, sd = 0.5), 
                                 rnorm(30, mean = 3, sd = 0.5),
                                 rnorm(15, mean = 4, sd = 0.5) 
                                 ),
                
                tail =  c(rnorm(60, mean = 2, sd = 0.5), 
                                 rnorm(30, mean = 3, sd = 0.5),
                                 rnorm(15, mean = 4, sd = 0.5) 
                                 ),
                velocity =  c(rnorm(60, mean = 2, sd = 0.5), 
                                 rnorm(30, mean = 3, sd = 0.5),
                                 rnorm(15, mean = 4, sd = 0.5) 
                                 )
                )
  d[,morph123 :=ifelse(morph == 'independent', 1, ifelse(morph == 'satelite', 2,3))]             

  colors_ <- colors[d$morph123]
  par(las = 1, cex.axis = 0.6, cex.lab = 0.8, cex.main = 0.8)
  s3d=scatterplot3d(d$midpiece, d$tail, d$velocity, pch = 16, type="h", 
              color=colors_, grid=TRUE, box=FALSE,
              xlab = "",
              ylab = "",
              zlab = "",
              x.ticklabs=c("short","","","","","long"),
              y.ticklabs=c("short","","","","long",""),
              z.ticklabs=c("slow","","","","","fast"),
              mar = c(3, 2, 0, 1.5)
              )     
  text(x = 7.5, y = 1, "Tail", srt = 0, cex = 0.8)
  text(x = 2.5, y = -0.5, "Mipiece", srt = 0,xpd = TRUE, cex = 0.8)
  text(x = -0.5, y = 2.5, "Velocity", srt = 90,xpd = TRUE, cex = 0.8)
  legend("bottom", legend = levels(factor(d$morph,levels = c('independent','satelite','faeder'))),
      col =  colors, pch = 16,xpd = TRUE, horiz = TRUE,inset = -0.125, bty = "n", cex = 0.7)
 
```
**Figure 3** | Predicted relationship between morph type and sperm traits.



# References