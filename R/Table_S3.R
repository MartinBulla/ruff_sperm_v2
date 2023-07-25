
# Tools 
 require(here)
 require(data.table)

# Load summary data from simmulatins DAT_relatedness_effects.R
 load('Outputs/resRel_morpho_5000.Rdata')
 load('Outputs/resRel_motil_5000.Rdata')
 load('Outputs/resRel_CV_5000.Rdata')

# Prepare and save
  wv[, data := 'June'] 
  wv[, type := 'motility'] 
  w[, data := 'means'] 
  w[, type := 'morphology'] 
  vv[, data := 'all']
  vv[, type := 'motility'] 
  v[, data := 'all']
  v[, type := 'morphology']
  wcv[, data := 'CV']
  wcv[, type := 'morphology'] 

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