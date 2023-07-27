  # =============================================================
  # ❗ Runs relative to the project's root directory and 
  # exports Table S4 into ./Outputs/
  # =============================================================

# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  require(ggnewscale)
  require(ggpubr) 
  require(ggsci)
  require(MuMIn)

  width_ = 0.6 # spacing between error bars

# DATA
  source(here::here('R/DAT_prepare.R'))       

  dd = d[!Morph%in%'Zebra finch']
    
  # long format
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

      aw_ = aw[, c('bird_ID','Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')]

      adw = merge(ddx, aw_, by = 'bird_ID')

      adwl = merge(ddxl, aw_, by = 'bird_ID')
      adwl[Motility == 'VAP' ,Motility:='Average path']
      adwl[Motility == 'VCL' ,Motility:='Curvilinear']
      adwl[Motility == 'VSL' ,Motility:='Straight line']
      adwl[, Motility := factor(Motility, levels=(c("Curvilinear", "Straight line", "Average path")))] 

  adw[,motileCount_ln_z := scale(log(motileCount))]
  adwl[,motileCount_ln_z := scale(log(motileCount))]
  ar[,motileCount_ln_z := scale(log(motileCount))]
 
# compare AIC for simple and quadratic
   l = list()
   for(i in unique(adwl$mot)){
        #i ='Straight line'
        ai = adwl[mot == i]
        
        ma = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Acrosome), ai)
        map = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Acrosome,2), ai)
        aic = data.table(AICc(ma, map))
        aic[, part:='Acrosome']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mn = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Nucleus), ai)
        mnp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Nucleus,2), ai)
        aic = data.table(AICc(mn, mnp))
        aic[, part:='Nucleus']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mm = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Midpiece), ai)
        mmp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Midpiece,2), ai)
        aic = data.table(AICc(mm, mmp))
        aic[, part:='Midpiece']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mt = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Tail), ai)
        mtp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Tail,2), ai)
        aic = data.table(AICc(mt, mtp))
        aic[, part:='Tail']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mo = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Total), ai)
        mop = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Total,2), ai)
        aic = data.table(AICc(mo, mop))
        aic[, part:='Total']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic
        
        mh = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Head), ai)
        mhp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Head,2), ai)
        aic = data.table(AICc(mh, mhp))
        aic[, part:='Head']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic
        

        mf = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Flagellum), ai)
        mfp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Flagellum,2), ai)
        aic = data.table(AICc(mf, mfp))
        aic[, part:='Flagellum']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mmr = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Midpiece_rel), ai)
        mmrp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Midpiece_rel,2), ai)
        aic = data.table(AICc(mmr, mmrp))
        aic[, part:='Midpiece_rel']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic

        mfr = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Flagellum_rel), ai)
        mfrp = lm(scale(value) ~ motileCount_ln_z+ Morph+poly(Flagellum_rel,2), ai)
        aic = data.table(AICc(mfr, mfrp))
        aic[, part:='Flagellum_rel']
        aic[, mot:=i]
        l[[paste(i, aic$part[1])]] = aic
       
        }    
   aic = do.call(rbind, l)
   #aic[df == 6, poly:='no']
   #aic[df == 7, poly:='yes']
   aic[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Midpiece","Tail","Total","Head","Flagellum","Midpiece_rel","Flagellum_rel")))] 
   aic[, mot := factor(mot, levels=rev(c("Curvilinear", "Straight line", "Average path")))] 
   aic = merge(aic[df==6], aic[df==7], by = c('part', 'mot'))
   setnames(aic, old = c('AICc.x', 'AICc.y'), new = c('simple', 'quadratic'))
   aic$df.x = aic$df.y = NULL
   aic = aic[,.(mot, part, simple, quadratic)]
   aic[,simple:=round(simple,2)]
   aic[,quadratic:=round(quadratic,2)]
   aic[, deltaAICc:=round(quadratic-simple,2)]
   aic[, prob := round(exp(-0.5*deltaAICc)/sum(exp(-0.5*deltaAICc)),2)]
   aic[, ER := round((1-prob)/prob, 2)]    
   aic = aic[order(mot,part, decreasing = TRUE)]   
   #aic_w = reshape(aic, idvar = c('part'), timevar = 'mot', direction = "wide")  
   fwrite(aic, file = 'Outputs/Table_S4.csv') 
   #fwrite(aic_w, file = 'Outputs/Table_SpolyAIC_w.csv')

# Not in the MS: compare AIC for simple and interactions
  l = list()
  for(i in unique(adwl$mot)){
    #i ='Curvilinear'
    ai = adwl[mot == i]
    ms = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Head)+scale(Midpiece)+scale(Tail), ai)
    mi = lm(scale(value) ~ motileCount_ln_z+ Morph+scale(Head)*scale(Midpiece)+scale(Head)*scale(Tail) + scale(Midpiece)*scale(Tail), ai)
    aic = data.table(AICc(ms,mi))
    aic[, mot:=i]
    l[[i]] = aic
  }
  aic = do.call(rbind, l)
  aic[, mot := factor(mot, levels=rev(c("Curvilinear", "Straight line", "Average path")))] 
  aic = merge(aic[df==8], aic[df==11], by = c('mot'))
  setnames(aic, old = c('AICc.x', 'AICc.y'), new = c('simple', 'interactions'))
  aic$df.x = aic$df.y = NULL
  aic[,simple:=round(simple,2)]
  aic[,interactions:=round(interactions,2)]
  aic[, deltaAICc:=round(interactions-simple,2)]
  aic[, prob := round(exp(-0.5*deltaAICc)/sum(exp(-0.5*deltaAICc)),2)]
  aic[, ER := round((1-prob)/prob, 2)]    
  aic = aic[order(mot, decreasing = TRUE)]   
  fwrite(aic, file = 'Outputs/Table_SintAIC.csv') 
    
# END