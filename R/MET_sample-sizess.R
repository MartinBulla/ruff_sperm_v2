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
  b[sample_ID==95, month:='June'] # for bird 1339 motility measured from May sample, but sperm from June, so metadata were changed to May sample to allow for correlation between the two, so here we change it back to get correct sample sizes\
  
# sample sizes
  # velocity
  v = unique(d[species!='zebrafinch', unique(bird_ID)])
  length(d[species=='zebrafinch', unique(bird_ID)]) # 5 zebra finch velocity measurements
  length(unique(d[species!='zebrafinch', unique(bird_ID)])) # N = 92 ruff males with velocity measurements
  nrow(d[species!='zebrafinch']) # 134 velocity measurements
  length(d[duplicated(bird_ID), bird_ID])  # 42 males with both June and May velocity measurements
  nrow(d[species!='zebrafinch' & month == 'June']) # 88 males recorded in June  
  nrow(d[species!='zebrafinch' & month == 'May']) # 46 males recorded in May  
  
  # morphology
  m = unique(b$bird_ID)
  length(m) # N = 92 males
  nrow(b[part == 'Acrosome']) # N = 920 sperm
  nrow(b) # N = 6440 measurements (92 males, 10 sperm/male, 7 sperm components)
  bb = b[!duplicated(sample_ID)]
  summary(factor(bb$month))  # 31 morpho samples from May, 61 from June used
  # END