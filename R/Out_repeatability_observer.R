  # TOOLS
    require(here)
    source(here::here('R/tools.R'))
    require(arm)
    require(effects)
    require(foreach)
    library(ggpubr)
    require(magrittr)
    require(PerformanceAnalytics)
    require(RColorBrewer)
    require(rptR) 
    require(stringi)
    

    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)
    # functions
        # for Repeatability output based on sim
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
        # mode output
         m_out = function(model = m, type = "mixed", 
            name = "define", dep = "define", fam = 'Gaussian',
            round_ = 3, nsim = 5000, aic = TRUE, save_sim = FALSE, N = NA, back_tran = FALSE, perc_ = 1){
              # perc_ 1 = proportion or 100%
            bsim = sim(model, n.sim=nsim)  
            
            if(save_sim!=FALSE){save(bsim, file = paste0(save_sim, name,'.RData'))}
           
            if(type != "mixed"){
             v = apply(bsim@coef, 2, quantile, prob=c(0.5))
             ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 

             if(back_tran == TRUE & fam == "binomial"){
              v = perc_*plogis(v)
              ci = perc_*plogis(ci)
             }
            if(back_tran == TRUE & fam == "binomial_logExp"){
                  v = perc_*(1-plogis(v))
                  ci = perc_*(1-plogis(ci))
                  ci = rbind(ci[2,],ci[1,])
                 }

             if(back_tran == TRUE & fam == "Poisson"){
              v = exp(v)
              ci = exp(ci)
             }

             oi=data.frame(type='fixed',effect=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,])
              rownames(oi) = NULL
              oi$estimate_r=round(oi$estimate,round_)
              oi$lwr_r=round(oi$lwr,round_)
              oi$upr_r=round(oi$upr,round_)
              if(perc_ == 100){
               oi$estimate_r = paste0(oi$estimate_r,"%")
               oi$lwr_r = paste0(oi$lwr_r,"%")
               oi$upr_r = paste0(oi$upr_r,"%")
              }
             x=data.table(oi[c('type',"effect", "estimate_r","lwr_r",'upr_r')]) 
           
            }else{
             v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
             ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 

             if(back_tran == TRUE & fam == "binomial"){
              v = perc_*plogis(v)
              ci = perc_*plogis(ci)
             }
            if(back_tran == TRUE & fam == "binomial_logExp"){
                  v = perc_*(1-plogis(v))
                  ci = perc_*(1-plogis(ci))
                  ci = rbind(ci[2,],ci[1,])
                 }

             if(back_tran == TRUE & fam == "Poisson"){
              v = exp(v)
              ci = exp(ci)
             }

             oi=data.frame(type='fixed',effect=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,])
                rownames(oi) = NULL
                oi$estimate_r=round(oi$estimate,round_)
                oi$lwr_r=round(oi$lwr,round_)
                oi$upr_r=round(oi$upr,round_)
                if(perc_ == 100){
                 oi$estimate_r = paste0(oi$estimate_r,"%")
                 oi$lwr_r = paste0(oi$lwr_r,"%")
                 oi$upr_r = paste0(oi$upr_r,"%")
                }
             oii=oi[c('type',"effect", "estimate_r","lwr_r",'upr_r')] 
            
             l=data.frame(summary(model)$varcor)
             l = l[is.na(l$var2),]
             l$var1 = ifelse(is.na(l$var1),"",l$var1)
             l$pred = paste(l$grp,l$var1)

             q050={}
             q025={}
             q975={}
             pred={}
             
             # variance of random effects
             for (ran in names(bsim@ranef)) {
               ran_type = l$var1[l$grp == ran]
               for(i in ran_type){
                q050=c(q050,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.5)))
                q025=c(q025,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.025)))
                q975=c(q975,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.975)))
                pred= c(pred,paste(ran, i))
                }
               }
             # residual variance
             q050=c(q050,quantile(bsim@sigma^2, prob=c(0.5)))
             q025=c(q025,quantile(bsim@sigma^2, prob=c(0.025)))
             q975=c(q975,quantile(bsim@sigma^2, prob=c(0.975)))
             pred= c(pred,'Residual')

             ri=data.frame(model = name,type='random %',effect=pred, estimate_r=round(100*q050/sum(q050)), lwr_r=round(100*q025/sum(q025)), upr_r=round(100*q975/sum(q975)))
               rx = ri[ri$effect == 'Residual',]
               if(rx$lwr_r>rx$upr_r){ri$lwr_r[ri$effect == 'Residual'] = rx$upr_r; ri$upr_r[ri$effect == 'Residual'] = rx$lwr_r}
               ri$estimate_r = paste0(ri$estimate_r,'%')
               ri$lwr_r = paste0(ri$lwr_r,'%')
               ri$upr_r = paste0(ri$upr_r,'%')
            
            x = data.table(rbind(oii,ri))
            }
            
            x[1, model := name]                                                                
            x[1, response := dep]                                                                
            x[1, error_structure := fam]                                                                
            x[1, N := N]                                                                

            x=x[ , c('model', 'response', 'error_structure', 'N', 'type',"effect", "estimate_r","lwr_r",'upr_r')] 

            if (aic == TRUE){   
                x[1, AIC := AIC(update(model,REML = FALSE))] 
                }
            if (aic == "AICc"){
                aicc = AICc(model)
                x[1, AICc := aicc] 
            }
            if(type == "mixed"){
              x[1, R2_mar := invisible({capture.output({r2_nakagawa(model)$R2_marginal})})]
              x[1, R2_con := invisible({capture.output({r2_nakagawa(model)$R2_conditional})})]
             }
            x[is.na(x)] = ""
            return(x)
          } 
  # DATA
    f = data.table(
      f = c(list.files(path = here::here('Data/'), pattern = 'test40', recursive = TRUE, full.names = TRUE)),
      f2 = c(list.files(path = here::here('Data/'), pattern = 'test40', recursive = TRUE, full.names = FALSE))
      )   

    ff = f[substr(f2,8,9) == 'MB']
    d2 = foreach(j = 1:nrow(ff), .combine = rbind) %do% {
      #j = 1
      k = ff[j,]
      x = fread(file = k[,f])
      x[, who := substr(k[,f2], 8, 9)]
      #x[, datetime_ := substr(j, 27, 45)]
      x[, id := substr(x$'File Name', 1, 2)]
      x[, pk := 1:nrow(x)]
      x[, pic := substr(k[,f2], 11, 13)]
      #x[, manip := substr(x$'File Name', 4, nchar(x$'File Name')-7-(nchar(pk)))]
      
      #table(x$id)
      xx = foreach(i = unique(x$id), .combine = rbind) %do% {
        #i = '01'
        xi = x[id == i]
        xi[, part := c('Acrosome','Nucleus','Midpiece','Tail')]
        return(xi)
        }    
      return(xx)
      }
      
    ff = f[substr(f2,8,9) == 'MC']
    d3 = foreach(j = 1:nrow(ff), .combine = rbind) %do% {
      #k = ff[5,]
      k = ff[j,]
      x = fread(file = k[,f])
      x[, who := substr(k[,f2], 8, 9)]
      #x[, datetime_ := substr(j, 27, 45)]
      x[, id := substr(x$'File Name', 1, 2)]
      x[, pk := 1:nrow(x)]
      x[, pic := substr(k[,f2], 13, 15)]
      xx = foreach(i = unique(x$id), .combine = rbind) %do% {
        #i = '02'
        xi = x[id == i]
        if(nrow(xi)==4){
         xi[Pixels == max(Pixels), part := 'Tail']
         xi[Pixels != max(Pixels), part := c('Acrosome','Nucleus','Midpiece')]
         } else{xi[, part := c('Acrosome','Nucleus','Midpiece')]}
         return(xi)
        }    
      return(xx)
      }
        
    d = rbind(d2,d3)  
    d[, id_part := paste(id, part)]
    d$pk = d$'File Name' = NULL
    # add composite measures
       d[, id_who_pic := paste(id, who, pic)]
       dt = data.table(ddply(d,.(who, id, pic, id_who_pic), summarise, part = 'Total', Pixels = sum(Pixels)))
       dt[,id_part := paste(id, part)]
       dt = dt[Pixels>1000] # excludes one trial with one sperm measured once
       
       dh = data.table(ddply(d[part %in% c('Acrosome','Nucleus')],.(who, id, pic, id_who_pic), summarise, part = 'Head', Pixels = sum(Pixels)))
       df = data.table(ddply(d[part %in% c('Midpiece','Tail')],.(who, id, pic, id_who_pic), summarise, part = 'Flagellum', Pixels = sum(Pixels)))
       df = df[Pixels>1000]

       d = rbind(d,dh,df,dt)

# within observer
  b = d[who=='MC']
  # estimate gre - des 
    lfrpt = list()
    for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
      part_ = i
      # part_ = "Acrosome"
      dd = b[part == part_ & pic %in% c('gre','des')]
      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(compar = 'MC gre-des', name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('compar','part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      RR[, pred := 100*R$R]
      RR[, lwr := 100*R$CI_emp[1]]
      RR[, upr := 100*R$CI_emp[2]]
      lfrpt[[i]] =  RR
      print(i)
    } 

    y1 = do.call(rbind,lfrpt)
  
  # estimate inv - des 
    lfrpt2 = list()
    for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
      part_ = i
      # part_ = "Acrosome"
      dd = b[part == part_ & pic %in% c('inv','des')]
      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(compar = 'MC inv-des', name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('compar','part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      RR[, pred := 100*R$R]
      RR[, lwr := 100*R$CI_emp[1]]
      RR[, upr := 100*R$CI_emp[2]]
      lfrpt2[[i]] =  RR
      print(i)
    } 

    y2 = do.call(rbind,lfrpt2)
  
  # estimate inv - gre  
    lfrpt3 = list()
    for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
      part_ = i
      # part_ = "Acrosome"
      dd = b[part == part_ & pic %in% c('inv','gre')]
      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(compar = 'MB inv-gre', name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('compar','part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      RR[, pred := 100*R$R]
      RR[, lwr := 100*R$CI_emp[1]]
      RR[, upr := 100*R$CI_emp[2]]
      lfrpt3[[i]] =  RR
      print(i)
    } 

    y3 = do.call(rbind,lfrpt3)
    
  # estimate all
    lfrpt4 = list()
    for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
      part_ = i
      # part_ = "Nucleus"
      dd = b[part == part_]
      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(compar = 'MC all', name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('compar','part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      RR[, pred := 100*R$R]
      RR[, lwr := 100*R$CI_emp[1]]
      RR[, upr := 100*R$CI_emp[2]]
      lfrpt4[[i]] =  RR
      print(i)
    } 

    y4 = do.call(rbind,lfrpt4)

y = rbind(y1,y2,y3,y4)

# between observers
    b = d[(who=='MC' & pic %in% c('des'))| who =='MB']
    l = list()
    for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
      part_ = i
      # part_ = "Acrosome"
      dd = b[part == part_]
      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(compar = 'MC-MB', name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('compar','part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      RR[, pred := 100*R$R]
      RR[, lwr := 100*R$CI_emp[1]]
      RR[, upr := 100*R$CI_emp[2]]
      l[[i]] =  RR
      print(i)
    } 
    yb = do.call(rbind,l)

    b = d[(who=='MC' & pic %in% c('inv'))| who =='MB']
    l = list()
    for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
      part_ = i
      # part_ = "Acrosome"
      dd = b[part == part_]
      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(compar = 'MC-MB inv', name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('compar','part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      RR[, pred := 100*R$R]
      RR[, lwr := 100*R$CI_emp[1]]
      RR[, upr := 100*R$CI_emp[2]]
      l[[i]] =  RR
      print(i)
    } 
    yc = do.call(rbind,l)

    b = d[(who=='MC' & pic %in% c('gre'))| who =='MB']
    l = list()
    for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
      part_ = i
      # part_ = "Acrosome"
      dd = b[part == part_]
      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(compar = 'MC-MB gre', name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('compar','part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      RR[, pred := 100*R$R]
      RR[, lwr := 100*R$CI_emp[1]]
      RR[, upr := 100*R$CI_emp[2]]
      l[[i]] =  RR
      print(i)
    } 
    yd = do.call(rbind,l)
  
#### from ruff_sperm

    # inverted only
      b = d[pic == "inv" & id_part %in% d[duplicated(id_part), id_part]] # only sperm ids and parts measured 
      
      # composite parts
           b[, id_who := paste(id, who)]
           bt = data.table(ddply(b,.(who, id, id_who), summarise, part = 'Total', Pixels = sum(Pixels)))
           bt = bt[Pixels>1000]
           bh = data.table(ddply(b[part %in% c('Acrosome','Nucleus')],.(who, id, id_who), summarise, part = 'Head', Pixels = sum(Pixels)))
           bf = data.table(ddply(b[part %in% c('Midpiece','Tail')],.(who, id, id_who), summarise, part = 'Flagellum', Pixels = sum(Pixels)))
           bf = bf[Pixels>1000]
           bc = rbind(bh,bf,bt)

      # wide format
           bw = reshape(b[,.(who,id,part,Pixels)], idvar = c('id','part'), timevar = 'who', direction = "wide")

           bcw = reshape(bc[,.(who,id,part,Pixels)], idvar = c('id','part'), timevar = 'who', direction = "wide")    
    # Maggie inverted vs desaturated vs gray scale
      # composite parts
           mm = d[who == 'MC']
           mm[, id_pic := paste(id, pic)]
           mt = data.table(ddply(mm,.(pic, id, id_pic), summarise, part = 'Total', Pixels = sum(Pixels)))
           mt = mt[Pixels>1000]
           mh = data.table(ddply(mm[part %in% c('Acrosome','Nucleus')],.(pic, id, id_pic), summarise, part = 'Head', Pixels = sum(Pixels)))
           mf = data.table(ddply(mm[part %in% c('Midpiece','Tail')],.(pic, id, id_pic), summarise, part = 'Flagellum', Pixels = sum(Pixels)))
           mf = mf[Pixels>1000]
           mc = rbind(mh,mf,mt)

      # wide format
           mw = reshape(mm[,.(pic,id,part,Pixels)], idvar = c('id','part'), timevar = 'pic', direction = "wide")

           mcw = reshape(mc[,.(pic,id,part,Pixels)], idvar = c('id','part'), timevar = 'pic', direction = "wide")    

  # Repeatability prep inv
    # all  
      # estimate BASIC 
              lfsim = list()
              lfrpt = list()
              for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
                part_ = i
                # part_ = "Acrosome"
                dd = b[part == part_]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
                RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
                RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
                lfrpt[[i]] =  RR[, method_CI := 'rpt package']
                #print(i)
              } 
              x = do.call(rbind,lfsim)
              names(x)[1] = "part"
              y = do.call(rbind,lfrpt)

              x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
              x[, lwr:= as.numeric(substr(CI,1,2))]
              x[, upr:= as.numeric(substr(CI,4,5))]

              y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
              y[nchar(CI) == 5, CI := paste0(0,CI) ]
              y[, lwr:= as.numeric(substr(CI,1,2))]
              y[, upr:= as.numeric(substr(CI,4,5))]
              names(y)[2] = tolower( names(y)[2])
              xy = rbind(x,y)
              xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
              xy[nchar(CI) == 7 & upr == 10, upr := 100]
      # estimate COMPOSITE
              lfsim = list()
              lfrpt = list()
              for(i in c("Head", "Flagellum","Total")){
                part_ = i
                # part_ = "Flagellum"
                dd = bc[part == part_]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
                RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
                RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
                lfrpt[[i]] =  RR[, method_CI := 'rpt package']
                #print(i)
              }
              
              x = do.call(rbind,lfsim)
              names(x)[1] = "part"
              y = do.call(rbind,lfrpt)

              x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
              x[, lwr:= as.numeric(substr(CI,1,2))]
              x[, upr:= as.numeric(substr(CI,4,5))]

              y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
              y[nchar(CI) == 5, CI := paste0(0,CI) ]
              y[, lwr:= as.numeric(substr(CI,1,2))]
              y[, upr:= as.numeric(substr(CI,4,5))]
              names(y)[2] = tolower( names(y)[2])
              xy2 = rbind(x,y)
              xy2[, part := factor(part, levels=c("Head","Flagellum","Total"))] 
              xy2[nchar(CI) == 7 & upr == 10, upr := 100]
    # KT, MB
      # estimate BASIC 
          lfsim = list()
          lfrpt = list()
          for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
            part_ = i
            # part_ = "Acrosome"
            dd = b[part == part_ & who %in% c('KT','MB')]
            m = lmer(Pixels ~ 1+(1|id), dd)
            Rf = R_out(part_)
            lfsim[[i]] = Rf[, method_CI:='arm package']

            R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
            RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
            RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
            lfrpt[[i]] =  RR[, method_CI := 'rpt package']
            print(i)
          } 
          x = do.call(rbind,lfsim)
          names(x)[1] = "part"
          y = do.call(rbind,lfrpt)

          x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
          x[, lwr:= as.numeric(substr(CI,1,2))]
          x[, upr:= as.numeric(substr(CI,4,5))]

          y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
          y[nchar(CI) == 5, CI := paste0(0,CI) ]
          y[, lwr:= as.numeric(substr(CI,1,2))]
          y[, upr:= as.numeric(substr(CI,4,5))]
          names(y)[2] = tolower( names(y)[2])
          xyKTMB = rbind(x,y)
          xyKTMB[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
          xyKTMB[nchar(CI) == 7 & upr == 10, upr := 100]
      # estimate COMPOSITE
          lfsim = list()
          lfrpt = list()
          for(i in c("Head", "Flagellum","Total")){
            part_ = i
            # part_ = "Flagellum"
            dd = bc[part == part_ & who %in% c('KT','MB')]
            m = lmer(Pixels ~ 1+(1|id), dd)
            Rf = R_out(part_)
            lfsim[[i]] = Rf[, method_CI:='arm package']

            R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
            RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
            RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
            lfrpt[[i]] =  RR[, method_CI := 'rpt package']
            print(i)
          }
          
          x = do.call(rbind,lfsim)
          names(x)[1] = "part"
          y = do.call(rbind,lfrpt)

          x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
          x[, lwr:= as.numeric(substr(CI,1,2))]
          x[, upr:= as.numeric(substr(CI,4,5))]

          y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
          y[nchar(CI) == 5, CI := paste0(0,CI) ]
          y[, lwr:= as.numeric(substr(CI,1,2))]
          y[, upr:= as.numeric(substr(CI,4,5))]
          names(y)[2] = tolower( names(y)[2])
          xy2KTMB = rbind(x,y)
          xy2KTMB[, part := factor(part, levels=c("Head","Flagellum","Total"))] 
          xy2KTMB[nchar(CI) == 7 & upr == 10, upr := 100]
    # MB, MC
      # estimate BASIC 
          lfsim = list()
          lfrpt = list()
          for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
            part_ = i
            # part_ = "Acrosome"
            dd = b[part == part_ & who %in% c('MC','MB')]
            m = lmer(Pixels ~ 1+(1|id), dd)
            Rf = R_out(part_)
            lfsim[[i]] = Rf[, method_CI:='arm package']

            R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
            RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
            RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
            lfrpt[[i]] =  RR[, method_CI := 'rpt package']
            print(i)
          } 
          x = do.call(rbind,lfsim)
          names(x)[1] = "part"
          y = do.call(rbind,lfrpt)

          x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
          x[, lwr:= as.numeric(substr(CI,1,2))]
          x[, upr:= as.numeric(substr(CI,4,5))]

          y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
          y[nchar(CI) == 5, CI := paste0(0,CI) ]
          y[, lwr:= as.numeric(substr(CI,1,2))]
          y[, upr:= as.numeric(substr(CI,4,5))]
          names(y)[2] = tolower( names(y)[2])
          xyMBMC = rbind(x,y)
          xyMBMC[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
          xyMBMC[nchar(CI) == 7 & upr == 10, upr := 100]
      # estimate COMPOSITE
          lfsim = list()
          lfrpt = list()
          for(i in c("Head", "Flagellum","Total")){
            part_ = i
            # part_ = "Flagellum"
            dd = bc[part == part_ & who %in% c('MC','MB')]
            m = lmer(Pixels ~ 1+(1|id), dd)
            Rf = R_out(part_)
            lfsim[[i]] = Rf[, method_CI:='arm package']

            R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
            RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
            RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
            lfrpt[[i]] =  RR[, method_CI := 'rpt package']
            print(i)
          }
          
          x = do.call(rbind,lfsim)
          names(x)[1] = "part"
          y = do.call(rbind,lfrpt)

          x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
          x[, lwr:= as.numeric(substr(CI,1,2))]
          x[, upr:= as.numeric(substr(CI,4,5))]

          y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
          y[nchar(CI) == 5, CI := paste0(0,CI) ]
          y[, lwr:= as.numeric(substr(CI,1,2))]
          y[, upr:= as.numeric(substr(CI,4,5))]
          names(y)[2] = tolower( names(y)[2])
          xy2MBMC = rbind(x,y)
          xy2MBMC[, part := factor(part, levels=c("Head","Flagellum","Total"))] 
          xy2MBMC[nchar(CI) == 7 & upr == 10, upr := 100]
    # KT, MC
      # estimate BASIC 
          lfsim = list()
          lfrpt = list()
          for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
            part_ = i
            # part_ = "Acrosome"
            dd = b[part == part_ & who %in% c('MC','KT')]
            m = lmer(Pixels ~ 1+(1|id), dd)
            Rf = R_out(part_)
            lfsim[[i]] = Rf[, method_CI:='arm package']

            R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
            RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
            RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
            lfrpt[[i]] =  RR[, method_CI := 'rpt package']
            print(i)
          } 
          x = do.call(rbind,lfsim)
          names(x)[1] = "part"
          y = do.call(rbind,lfrpt)

          x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
          x[, lwr:= as.numeric(substr(CI,1,2))]
          x[, upr:= as.numeric(substr(CI,4,5))]

          y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
          y[nchar(CI) == 5, CI := paste0(0,CI) ]
          y[, lwr:= as.numeric(substr(CI,1,2))]
          y[, upr:= as.numeric(substr(CI,4,5))]
          names(y)[2] = tolower( names(y)[2])
          xyKTMC = rbind(x,y)
          xyKTMC[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
          xyKTMC[nchar(CI) == 7 & upr == 10, upr := 100]
      # estimate COMPOSITE
          lfsim = list()
          lfrpt = list()
          for(i in c("Head", "Flagellum","Total")){
            part_ = i
            # part_ = "Flagellum"
            dd = bc[part == part_ & who %in% c('MC','KT')]
            m = lmer(Pixels ~ 1+(1|id), dd)
            Rf = R_out(part_)
            lfsim[[i]] = Rf[, method_CI:='arm package']

            R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
            RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
            RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
            lfrpt[[i]] =  RR[, method_CI := 'rpt package']
            print(i)
          }
          
          x = do.call(rbind,lfsim)
          names(x)[1] = "part"
          y = do.call(rbind,lfrpt)

          x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
          x[, lwr:= as.numeric(substr(CI,1,2))]
          x[, upr:= as.numeric(substr(CI,4,5))]

          y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
          y[nchar(CI) == 5, CI := paste0(0,CI) ]
          y[, lwr:= as.numeric(substr(CI,1,2))]
          y[, upr:= as.numeric(substr(CI,4,5))]
          names(y)[2] = tolower( names(y)[2])
          xy2KTMC = rbind(x,y)
          xy2KTMC[, part := factor(part, levels=c("Head","Flagellum","Total"))] 
          xy2KTMC[nchar(CI) == 7 & upr == 10, upr := 100]              

  # Repeatability prep MC acrosome
    ma = mm[part == 'Acrosome'] 
    # all  
      part_ = 'Acrosome'
      dd = ma[part == part_]
      m = lmer(Pixels ~ 1+(1|id), dd)
      Rf = R_out(part_)
      x = Rf[, method_CI:='arm package']
      names(x)[1] = "part"

      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      y=  RR[, method_CI := 'rpt package']
      #print(i)

      x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
      x[, lwr:= as.numeric(substr(CI,1,2))]
      x[, upr:= as.numeric(substr(CI,4,5))]

      y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
      y[nchar(CI) == 5, CI := paste0(0,CI) ]
      y[, lwr:= as.numeric(substr(CI,1,2))]
      y[, upr:= as.numeric(substr(CI,4,5))]
      names(y)[2] = tolower( names(y)[2])
      xy_ma = rbind(x,y)
      #xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
      xy_ma[nchar(CI) == 7 & upr == 10, upr := 100]     
    # inv, des
      part_ = 'Acrosome'
      dd = ma[part == part_ & pic %in%c('inv','des')]
      m = lmer(Pixels ~ 1+(1|id), dd)
      Rf = R_out(part_)
      x = Rf[, method_CI:='arm package']
      names(x)[1] = "part"

      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      y=  RR[, method_CI := 'rpt package']
      #print(i)
       
      x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
      x[, lwr:= as.numeric(substr(CI,1,2))]
      x[, upr:= as.numeric(substr(CI,4,5))]

      y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
      y[nchar(CI) == 5, CI := paste0(0,CI) ]
      y[, lwr:= as.numeric(substr(CI,1,2))]
      y[, upr:= as.numeric(substr(CI,4,5))]
      names(y)[2] = tolower( names(y)[2])
      xy_ma_id = rbind(x,y)
      #xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
      xy_ma_id[nchar(CI) == 7 & upr == 10, upr := 100]
    # inv, gre
      part_ = 'Acrosome'
      dd = ma[part == part_ & pic %in%c('inv','gre')]
      m = lmer(Pixels ~ 1+(1|id), dd)
      Rf = R_out(part_)
      x = Rf[, method_CI:='arm package']
      names(x)[1] = "part"

      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      y=  RR[, method_CI := 'rpt package']
      #print(i)
       
      x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
      x[, lwr:= as.numeric(substr(CI,1,2))]
      x[, upr:= as.numeric(substr(CI,4,5))]

      y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
      y[nchar(CI) == 5, CI := paste0(0,CI) ]
      y[, lwr:= as.numeric(substr(CI,1,2))]
      y[, upr:= as.numeric(substr(CI,4,5))]
      names(y)[2] = tolower( names(y)[2])
      xy_ma_ig = rbind(x,y)
      #xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
      xy_ma_ig[nchar(CI) == 7 & upr == 10, upr := 100]
    # des, gre
      part_ = 'Acrosome'
      dd = ma[part == part_ & pic %in%c('des','gre')]
      m = lmer(Pixels ~ 1+(1|id), dd)
      Rf = R_out(part_)
      x = Rf[, method_CI:='arm package']
      names(x)[1] = "part"

      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      y=  RR[, method_CI := 'rpt package']
      #print(i)
       
      x[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
      x[, lwr:= as.numeric(substr(CI,1,2))]
      x[, upr:= as.numeric(substr(CI,4,5))]

      y[, pred:= as.numeric(substr(Repeatability,1,nchar(Repeatability)-1))]
      y[nchar(CI) == 5, CI := paste0(0,CI) ]
      y[, lwr:= as.numeric(substr(CI,1,2))]
      y[, upr:= as.numeric(substr(CI,4,5))]
      names(y)[2] = tolower( names(y)[2])
      xy_ma_dg = rbind(x,y)
      #xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
      xy_ma_dg[nchar(CI) == 7 & upr == 10, upr := 100]
