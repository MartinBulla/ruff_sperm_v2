# =============================================================
# ❗ The script runs relative to the project's root directory,
# and shows sample sizes used in Methods of the paper
# =============================================================

# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(ggpubr)

# DATA
  source(here::here('R/DAT_prepare.R'))       
  b[sample_ID==95, month:='June'] # for bird 1339 motility measured from May sample, but sperm from June, so metadata were changed to May sample to allow for correlation between the two, so here we change it back to get correct sample sizes\
  
# sample sizes
  # velocity
  length(d[species=='zebrafinch', unique(bird_ID)]) # 5 zebra finch velocity measurements
  length(unique(d[species!='zebrafinch', unique(bird_ID)])) # N = 92 ruff males with velocity measurements
  nrow(d[species!='zebrafinch']) # 134 velocity measurements
  length(d[duplicated(bird_ID), bird_ID])  # 42 males with both June and May velocity measurements
  nrow(d[species!='zebrafinch' & month == 'June']) # 88 males recorded in June  
  nrow(d[species!='zebrafinch' & month == 'May']) # 46 males recorded in May  

  # number of recorded sperm
  v = unique(d[species!='zebrafinch'])
  summary(v$motileCount)
  summary(v$motileCount[v$Morph%in%'Independent'])
  length(v$motileCount[v$Morph%in%'Independent'])
  summary(v$motileCount[v$Morph%in%'Satellite'])
  length(v$motileCount[v$Morph%in%'Satellite'])
  summary(v$motileCount[v$Morph%in%'Faeder'])
  length(v$motileCount[v$Morph%in%'Faeder'])

  ggplot(v, aes(x = Morph, y = motileCount )) +
  geom_violin(alpha = 0.5) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5) +
  theme(legend.position = "none")
  
  # morphology
  m = unique(b$bird_ID)
  length(m) # N = 92 males
  nrow(b[part == 'Acrosome']) # N = 920 sperm
  nrow(b) # N = 6440 measurements (92 males, 10 sperm/male, 7 sperm components)
  bb = b[!duplicated(sample_ID)]
  summary(factor(bb$month))  # 31 morpho samples from May, 61 from June used
  

# CV
  source(here::here('R/DAT_prepare.R')) 

  b5 = b[,.SD[sample(.N, min(5,.N))],by = list(bird_ID, part)]
  cv_5 =  b5[, cv(Length_µm), by = list(bird_ID, part, Morph, loc_avi)]
  cv_5[ , Morph123 := as.numeric(Morph)]
  setnames(cv_5, 'V1', 'CV5')
  cv_5 = merge(cv_5,z[,.(sampleid, HL)], by.x = 'bird_ID', by.y = 'sampleid')

  cv2 = merge(cv_, cv_5[,.(bird_ID, part, CV5)], by = c('bird_ID', 'part'))

  ggplot(cv2, aes(x = CV5, y = CV)) + 
    facet_wrap(~part, ncol = 4) + 
    geom_point() + 
    stat_cor(cor.coef.name = "r", aes(label = after_stat(r.label)), r.accuracy = 0.1, cex = 3) +
    geom_abline(slope = 1,
              intercept = 0,
              lty = 3, col = 'red')+
    labs(x = 'Coefficient of varation based on five sperm', y = 'Coefficient of varation based on ten sperm')


# END