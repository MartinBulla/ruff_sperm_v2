# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(arm)
  require(ggpubr) 
  require(ggsci)
# DATA
  source(here::here('R/DAT_prepare.R'))       

# Prepare estimates
  # motility - only one measure per bird (June)
    d[, motileCount_ln:=scale(log(motileCount))]
    dd = d[!Morph%in%'Zebra finch']
   
    # use June values and for 4 males without June, May
      dd1 = dd[month == 'June']
      dd2 = dd[month == 'May']
      ddx = rbind(dd1,dd2[!bird_ID%in%dd1$bird_ID])
  # extract p-values via multcomp package 
    l = list()
    # VAP
      m = lm(scale(VAP) ~ -1+ scale(log(motileCount)) + Morph, ddx)
      l[["VAP"]] = data.frame(response = "Average-path", p_value = as.numeric(summary(glht(m, linfct = c(
          "MorphIndependent - MorphSatellite = 0",
          "MorphIndependent - MorphFaeder = 0",
          "MorphSatellite - MorphFaeder = 0"
          )))$test$pvalues))
    # VSL
      m = lm(scale(VSL) ~ -1 + scale(log(motileCount)) + Morph, ddx)
      l[["VSL"]] = data.frame(response = "Straight line", p_value = as.numeric(summary(glht(m, linfct = c(
               "MorphIndependent - MorphSatellite = 0",
               "MorphIndependent - MorphFaeder = 0",
               "MorphSatellite - MorphFaeder = 0"
           )))$test$pvalues))
    # VCL
      m = lm(scale(VCL) ~ -1 + scale(log(motileCount)) + Morph, ddx)
      l[["VCL"]] <- data.frame(response = "Curvilinear", p_value = as.numeric(summary(glht(m, linfct = c(
          "MorphIndependent - MorphSatellite = 0",
          "MorphIndependent - MorphFaeder = 0",
          "MorphSatellite - MorphFaeder = 0"
      )))$test$pvalues))
  # morpho - averages
    for(i in unique(a$part)){
      #i ='Nucleus'
      m = lm(scale(Length_avg) ~ -1 + Morph, a[part == i])
      l[[i]] <- data.frame(response = i, p_value = as.numeric(summary(glht(m, linfct = c(
                "MorphIndependent - MorphSatellite = 0",
                "MorphIndependent - MorphFaeder = 0",
                "MorphSatellite - MorphFaeder = 0"
            )))$test$pvalues))
      print(i)     
      }          
    for(i in unique(ar$part)){
      if(i == 'Midpiece'){ii = 'Midpiece_rel'}
      if(i == 'Flagellum'){ii = 'Flagellum_rel'}
      m = lm(scale(Length_rel) ~ -1 + Morph, ar[part == i])
      l[[i]] <- data.frame(response = i, p_value = as.numeric(summary(glht(m, linfct = c(
          "MorphIndependent - MorphSatellite = 0",
          "MorphIndependent - MorphFaeder = 0",
          "MorphSatellite - MorphFaeder = 0"
      )))$test$pvalues))
      print(i)     
      }          
    for(i in unique(cv_$part)){
      #i ='Nucleus'
      m = lm(scale(CV) ~ -1 + Morph, cv_[part == i])
      l[[paste(i,'cv')]] <- data.frame(response = paste(i,'cv'), p_value = as.numeric(summary(glht(m, linfct = c(
          "MorphIndependent - MorphSatellite = 0",
          "MorphIndependent - MorphFaeder = 0",
          "MorphSatellite - MorphFaeder = 0"
      )))$test$pvalues))
      print(paste(i,'cv'))     
      }          
    
    t = data.table(do.call(rbind,l) ) 

# plot
  g = 
  ggplot(t, aes(x = p_value)) + geom_histogram()    
  ggsave(file = 'Outputs/hist_pvalues.png', width = 5, height = 5, unit = 'cm')
# END