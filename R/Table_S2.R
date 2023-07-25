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
       dh[,id_part := paste(id, part)]

       df = data.table(ddply(d[part %in% c('Midpiece','Tail')],.(who, id, pic, id_who_pic), summarise, part = 'Flagellum', Pixels = sum(Pixels)))
       df = df[Pixels>1000]
       df[,id_part := paste(id, part)]

       d = rbind(d,dh,df,dt)
       d = d[order(id_who_pic, id_part)]

# prepare repeatability estimates
  # within observer
    b = d[who=='MC']
    lw = list()
    for(i in unique(b$part)){
      part_ = i
      # part_ = "Acrosome"
      dd = b[part == part_ & pic %in% c('inv','des')]
      R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(compar = 'MC inv-des', name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('compar','part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      RR[, pred := 100*R$R]
      RR[, lwr := 100*R$CI_emp[1]]
      RR[, upr := 100*R$CI_emp[2]]
      lw[[i]] =  RR
      print(i)
    } 

    yw = do.call(rbind,lw)
  # between observers
      b = d[pic %in% c('inv')]
      lb = list()
      for(i in unique(b$part)){
        part_ = i
        # part_ = "Acrosome"
        dd = b[part == part_]
        R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
        RR = data.table(merge(data.frame(compar = 'MC-MB inv', name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('compar','part', 'Repeatability'))
        RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
        RR[, pred := 100*R$R]
        RR[, lwr := 100*R$CI_emp[1]]
        RR[, upr := 100*R$CI_emp[2]]
        lb[[i]] =  RR
        print(i)
      } 
      yb = do.call(rbind,lb)
  # combine
    yw[,type:='within observer']
    yb[,type:='between observer']
    y = rbind(yw,yb) 
    
# savee
  y[, Repeatability := format(round(pred), nsmall = 0)]
  y[,lwr_f := format(round(lwr,1), nsmall = 1)]
  y[,upr_f := format(round(upr,1), nsmall = 1)]
  y[, es_ci := paste0(Repeatability, '% (', lwr_f,' - ',upr_f,')')]
  y_w = reshape(y[,.(part, type, es_ci)], idvar = c('part'), timevar = 'type', direction = "wide") 
  fwrite(y_w, file = 'Outputs/Table_S2.csv')

# not in the MS: Figure SRo
y[, part := factor(part, levels = rev(c("Acrosome", "Nucleus", "Midpiece", "Tail", "Total", "Head", "Flagellum")))]
col_ <- c(
  viridis(1, alpha = 1, begin = 0.2, end = 0.2, direction = 1, option = "D"),
  viridis(1, alpha = 1, begin = 0.9, end = 0.9, direction = 1, option = "D")
)
tx <- 0.15
# show_col(col_)
gm <-
  ggplot(y, aes(y = part, x = pred, shape = type, col = type)) +
  geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, position = position_dodge(width = 0.6)) +
  # ggtitle ("Sim based")+
  geom_point(position = position_dodge(width = 0.6), bg = "white", size = 1.1) +
  # scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
  # scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) +
  scale_x_continuous(limits = c(80, 100), expand = c(0, 0)) +
  scale_shape_manual(guide = guide_legend(reverse = TRUE), values = c(23, 21)) +
  scale_color_manual(guide = guide_legend(reverse = TRUE), values = col_) +
  labs(y = NULL, x = "Repeatability [%]") +
  geom_text(x = 93, y = 4 + tx, label = "within observer", col = col_[2], size = 2.2, adj = 1) +
  geom_text(x = 93, y = 4 - tx, label = "between observer", col = col_[1], size = 2.2, adj = 1) +
  # coord_flip()+
  theme_bw() +
  theme(
    plot.subtitle = element_text(size = 9, color = "grey30"),
    plot.margin = margin(3, 3, 1, 1, "mm"),
    panel.border = element_rect(color = "grey70"),
    legend.position = "none",
    # legend.title = element_blank(),
    # legend.text=element_text(size=7.5, color = 'grey30'),
    # legend.key.height= unit(0.2,"line"),
    # legend.margin=margin(0,0,0,0),
    # legend.position=c(0.5,1.6),

    axis.ticks = element_blank(),
    axis.title = element_text(size = 10, color = "grey10")
  )

ggsave(here::here("Outputs/Fig_Ro_width-43mm.png"), gm, width = 4.3 / (5 / 7), height = 6, units = "cm", dpi = 600)

# END