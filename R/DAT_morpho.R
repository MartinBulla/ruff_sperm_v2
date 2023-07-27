# =============================================================
# ❗ The script runs relative to the project's root directory,
# links the randomized and blinded sperm measurements with 
# the individual males and exports the DAT_morphometrics.csv 
# dataset to ./DATA/
# =============================================================

# TOOLS & DATA
  require(here)
  require(data.table)
  source(here::here('R/tools.R'))

      packages = c("data.table", "ggplot2", "ggthemes", "glue", "grid", "gridExtra", "htmlTable", "lattice", "lubridate", "magrittr", "maptools", "plyr", "png", "raster", "RColorBrewer", "readxl", "scales", "scatterplot3d", "stringr", "zoo", "viridis", "writexl")
      sapply(packages, function(x) suppressPackageStartupMessages(using(x)))

  # CHECK WHICH package has %do%
  packages = c('anytime','data.table', 'DataEntry.validation', 'DT', 'foreach', 'ggplot2', 'ggthemes', 'glue','googledrive', 'googlesheets4', 'grid', 'htmlTable', 'lattice', 'lubridate', 'magrittr', 'maptools', 'openxlsx','plyr','raster','readxl','stringr','zoo')
> sapply(packages, function(x) suppressPackageStartupMessages(using(x)) )
   
# DATA 
  p = fread('Data/DAT_randomised_pics.csv') #p = fread(here::here('R/all_randomized_2022-03-21.csv')) # reads id as integer, instead of character
  p[, pic :=as.character(id)]
  p[nchar(pic)==1, pic := paste0('00',pic)]
  p[nchar(pic)==2, pic := paste0('0',pic)]
  p$id = NULL
  
  # extracted whether pics further manipulated before measuerement
  m1 = data.table(
    f = c(list.files(path = here::here('Data/sperm_morpho/'), pattern = '.erT', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("Data/sperm_morpho/"), pattern = ".erT", recursive = TRUE, full.names = FALSE))
    )
  m1[, pic := substr(m1$f2, 28, 30)]
  m2 = data.table(
    f = c(list.files(path = here::here('Data/sperm_morpho/'), pattern = '.erA', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("Data/sperm_morpho/"), pattern = ".erA", recursive = TRUE, full.names = FALSE))
    )
  m2[, pic := substr(m2$f2, 28, 30)]

  m3 = data.table(
    f = c(list.files(path = here::here('Data/sperm_morpho/'), pattern = '.erN', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("Data/sperm_morpho/"), pattern = ".erN", recursive = TRUE, full.names = FALSE))
    )
  m3[, pic := substr(m3$f2, 28, 30)]

  m4 = data.table(
    f = c(list.files(path = here::here('Data/sperm_morpho/'), pattern = '.erM', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("Data/sperm_morpho/"), pattern = ".erM", recursive = TRUE, full.names = FALSE))
    )
  m4[, pic := substr(m4$f2, 28, 30)]

# get and prepare measurements
  d = data.table(
    f = c(list.files(path = here::here('Data/sperm_morpho/'), pattern = '.csv', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("Data/sperm_morpho/"), pattern = ".csv", recursive = TRUE, full.names = FALSE))
    )
  d[ , measured := substring(d$f2,29,47)]

  b = foreach(j = 1:nrow(d), .combine = rbind) %do% {
       #j =1
       ff = d[j, ]
       
       x = fread(ff$f, stringsAsFactors = FALSE) #skip = skip, col.names = varnames, colClasses = colclass, 
       x[, measured := ff$measured]
       x[, pic := substr(x$'File Name', 1, 3)]
       if(j==8){
        x=x[!pic%in%363]
        x[, part := rep(c('Acrosome','Nucleus','Midpiece','Tail'),49)]
        }else if(j==9){x[, part := c('Acrosome','Nucleus','Midpiece','Tail')]
          }else if(j==20){x[, part := rep(c('Acrosome','Nucleus','Midpiece','Tail'),20)]
          }else{x[, part := rep(c('Acrosome','Nucleus','Midpiece','Tail'),50)]}
       x[pic %in% m1$pic & part == 'Tail', manip := 'erT']
       x[pic %in% m2$pic & part == 'Acrosome', manip := 'erA']
       x[pic %in% m3$pic & part == 'Nucleus', manip := 'erN']
       x[pic %in% m4$pic & part == 'Midpiece', manip := 'erM']
       print(j)
       return(x)
   }

  bp =  merge(b,p, all.x = TRUE)

  fwrite(bp, file = 'Data/DAT_morphometrics.csv')
# END