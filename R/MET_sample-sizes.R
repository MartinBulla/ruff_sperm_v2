# =============================================================
# ❗ The script runs relative to the project's root directory,
# and shows sample sizes used in Methods of the paper
# =============================================================

# TOOLS
  require(here)
  source(here::here('R/tools.R'))

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