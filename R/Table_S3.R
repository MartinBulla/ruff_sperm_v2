  # =============================================================
  # ❗ Runs relative to the project's root directory, uses
  # data generated by DAT_relatedness_effects.R
  # and exports Table S3 into ./Outputs/
  # =============================================================
  
# Tools 
 require(here)
 require(data.table)

# Load summary data from simmulatins DAT_relatedness_effects.R
 w = fread("Data/DAT_resRel_morpho_avg_5000.csv")
 v = fread("Data/DAT_resRel_morpho_single_5000.csv")
 
 wcv = fread("Data/DAT_resRel_CV_5000.csv")
 
 wv <- fread("Data/DAT_resRel_motil_June_5000.csv")
  
 vv = fread("Data/DAT_resRel_motil_single_5000.csv")

 #load("Outputs/freeze/resRel_morpho_5000.Rdata")
 #load("Outputs/freeze/resRel_motil_5000.Rdata")
 #load("Outputs/freeze/resRel_CV_5000.Rdata")

# Prepare and save
  w[, data := "means"]
  w[, type := "morphology"]
  v[, data := "all"]
  v[, type := "morphology"]

  wcv[, data := "CV"]
  wcv[, type := "morphology"]

  wv[, data := 'June'] 
  wv[, type := 'motility'] 
  
  vv[, data := 'all']
  vv[, type := 'motility'] 

  x = rbind(v,vv,w,wv, wcv)
  cols_=c('Rhat_no', 'Rhat_yes')
  x[ , (cols_) := lapply(.SD, function(x){round(x,2)}), .SDcols = cols_]
  x[, prob:=round(prob,3)]
  x[, bf:=round(bf)]
  cols_=c('est', 'lwr','upr')
  #x[ , (cols_) := lapply(.SD, function(x){round(x*100,1)}), .SDcols = cols_]
  #x[substring(upr,1,1) != 0, upr := round(upr)]
  #x[substring(lwr,1,1) != 0, lwr := round(lwr)]
  #x[, est := paste0(est, '%')]
  #x[, CI := paste0(lwr, '-', upr, '%')]
  
  x[ , (cols_) := lapply(.SD, function(x){format(round(x*100,1), nsmall = 1)}), .SDcols = cols_]
  x[, est := paste0(est, '%')]
  x[, lwr := paste0(lwr,'%')]
  x[, upr := paste0(upr,'%')]

  y = x[,.(data, type, response, est, lwr, upr, bf, prob, Rhat_no, Rhat_yes)]
  y[, data := factor(data, levels = c("June", "means", "CV", "all"))]
  y[, type := factor(type, levels=c("motility","morphology"))] 
  y[, response := factor(response, levels=c("Curvilinear",'Straight line',"Average path","Acrosome", "Nucleus",  "Midpiece","Tail","Total", "Head","Flagellum"))] 

  #summary(x[data=='all',.(est, lwr, upr, bf, prob, Rhat_no, Rhat_yes)])
  #summary(x[data!='all',.(est, lwr, upr, bf, prob, Rhat_no, Rhat_yes)])

  fwrite(y[order(data, type,response)], file = paste0('Outputs/Table_S3',sample_,'.csv'))
 
# END