# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))

# April_MAY_2021 
  d = data.table(
    f = c(list.files(path = here::here('April_May_2021/Photos_PRE_sample_ID_MAY/'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here('April_May_2021/Photos_PRE_sample_ID_MAY/'), pattern = '.jpg', recursive = TRUE, full.names = FALSE))
    )
  d = d[!grepl('metadata', f2, fixed = TRUE)]
  d[ , sample_ID := substring(f2,6,8)]
  d[, file_name :=sub(".*/", "", f2)]
  #d[nchar(file_name)<32]
  d[sample_ID %in% c('400'), file_name := paste(substring(file_name,1,8), '2021-04-19', substring(file_name,10))]
  d[sample_ID %in% c('401'), file_name := paste(substring(file_name,1,9), '2021-04-19', substring(file_name,11))]
  d[sample_ID %in% c('402'), file_name := paste(substring(file_name,1,8), '2021-04-20', substring(file_name,11))]
  
  #paste(substring('Ruff 320 Snap-591.jpg',1,8), '2021-04-20', substring('Ruff 320 Snap-591.jpg',10))

  #dd = d[sample_ID %in% c('400')]
  d[, new_name := paste0("Ruff ", sample_ID, "_", substring(file_name, 6))]
  d[substring(new_name, nchar(new_name)-8, nchar(new_name)-7) == "p-", new_name := gsub("p-", "p-0", new_name)]
  
  #d$new_name
  
  # rename & copy
  for(i in 1:nrow(d)){
    #i = 1 
    file.copy(from = d$f[i], to = glue('all_photos_original/',d$new_name[i])) 
  }

# June_2021 
  d = data.table(
    f = c(list.files(path = here::here('June_2021/Photos_PRE_sample_ID_no'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here('June_2021/Photos_PRE_sample_ID_no'), pattern = '.jpg', recursive = TRUE, full.names = FALSE))
  )
  d = d[!grepl('metadata', f2, fixed = TRUE)]
  d[ , sample_ID := substring(f2,6,8)]
  d[, file_name :=sub(".*/", "", f2)]
  #d[nchar(file_name)<32]
  
  d[, new_name := paste0("Ruff ", sample_ID, "_", substring(file_name, 6))]
  #d[substring(new_name, nchar(new_name)-8, nchar(new_name)-7) == "p-", new_name := gsub("p-", "p-0", new_name)]
  
  # rename & copy
  for(i in 1:nrow(d)){
    #i = 1 
    file.copy(from = d$f[i], to = glue('all_photos_original/',d$new_name[i])) 
  }

 # test
  d = data.table(
    f = c(list.files(path = here::here('all_photos_original'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here('all_photos_original'), pattern = '.jpg', recursive = TRUE, full.names = FALSE))
  )
  d[ , sample_ID := substring(f2,6,8)]
  nrow(d)
  unique(d$sample_ID)
  length(unique(d$sample_ID))

# REST with correct names
  d = data.table(
    f = c(list.files(path = here::here('June_2021/Photos_WITH_sample_ID_no'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here('June_2021/Photos_WITH_sample_ID_no'), pattern = '.jpg', recursive = TRUE, full.names = FALSE))
  )
  d = d[!grepl('metadata', f2, fixed = TRUE)]
  d[ , sample_ID := substring(f2,6,8)]
  d[, file_name :=sub(".*/", "", f2)]
  #d[nchar(file_name)<32]
  
  # test whether abnormal sperm sample present
  d[grepl('abn', f2, fixed = TRUE)]

  #d[, new_name := paste0("Ruff ", sample_ID, "_", substring(file_name, 6))]
  #d[substring(new_name, nchar(new_name)-8, nchar(new_name)-7) == "p-", new_name := gsub("p-", "p-0", new_name)]
  
  
  # copy
  for(i in 1:nrow(d)){
    #i = 1 
    file.copy(from = d$f[i], to = glue('all_photos_original/',d$file_name[i])) 
  }

 # test
  d = data.table(
    f = c(list.files(path = here::here('all_photos_original'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("all_photos_original"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  d[ , sample_ID := substring(f2,6,8)]
  nrow(d) # 1990
  unique(d$sample_ID)
  length(unique(d$sample_ID))
  summary(factor(d$sample_ID))
  dd = d[ ,  .N , by = sample_ID ]
  dd[N<10]
  o = d

# test 2
n <- data.table(
  f = c(list.files(path = here::here("all_photos_inv/inverted_all"), pattern = ".jpg", recursive = TRUE, full.names = TRUE)),
  f2 = c(list.files(path = here::here("all_photos_inv/inverted_all"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
)

nrow(n) # 2066
n[, inv := substring(f2, nchar(f2) - 6, nchar(f2) - 4)]
n[inv != "inv", .(f2)]
n[inv != "inv", pic := substring(f2, 1, nchar(f2) - 4)]

n[, pic := substring(f2, 1, nchar(f2) - 8)]

o[, pic := substring(f2, 1, nchar(f2) - 4)]
n[!pic%in%o$pic,.(f2)]

# RANDOMIZE and RENAME

  d1 = data.table(
    f = c(list.files(path = here::here('all_photos_inv/1'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("all_photos_inv/1"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  d2 = data.table(
    f = c(list.files(path = here::here('all_photos_inv/2'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("all_photos_inv/2"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  d3 = data.table(
    f = c(list.files(path = here::here('all_photos_inv/3'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("all_photos_inv/3"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  d4 = data.table(
    f = c(list.files(path = here::here('all_photos_inv/4'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("all_photos_inv/4"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  d5 = data.table(
    f = c(list.files(path = here::here('all_photos_inv/5'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("all_photos_inv/5"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  d6 = data.table(
    f = c(list.files(path = here::here('all_photos_inv/6'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("all_photos_inv/6"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  d7 = data.table(
    f = c(list.files(path = here::here('all_photos_inv/7'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("all_photos_inv/7"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  d = rbind(d1,d2,d3,d4,d5,d6,d7)
  d[ , sample_ID := substring(f2,6,8)]
  d[, file_name :=sub(".*/", "", f2)]
  d[, bird_ID :=substring(f2, 10,nchar(f2)-29)]
  d[, bird_ID :=substring(f2, 10,nchar(f2)-29)]
  d[bird_ID=="", bird_ID :=substring(f2, 10,nchar(f2)-25)]
  #substring('Ruff 018_1307 2021-11-10 Snap-3277.jpg', 10,nchar('Ruff 018_1307 2021-11-10 Snap-3277.jpg')-25)
  d[bird_ID == 'G200067 2', bird_ID := 'G200067']
  length(unique(d$bird_ID)) # 1 - 14, 2 - 16, 3 - 26, 4 - 16, 5 - 23
  #d[nchar(file_name)<32]
  
  # TEST
    d[, inv := substring(f2, nchar(f2) - 6, nchar(f2) - 4)]
    d[inv != "inv", .(f2)]

    dd = d[ ,  .N , by = sample_ID ]
    #dd
    dd[N<10]

    dd = d[ ,  .N , by = bird_ID ]
    #dd
    dd[N<10]
    dd[N>10]
    d[bird_ID%in%dd$bird_ID[dd$N>10],.N, by = list(bird_ID,sample_ID)]

  # add randomization
    d[, id := as.character(sample.int(n = nrow(d), size = nrow(d)))]
    d[nchar(id)==1, id := paste0('00',id)]
    d[nchar(id)==2, id := paste0('0',id)]
    d =d[order(id)]
    fwrite(d, file = 'R/all_randomized_2022-03-21.csv')

    for(i in 1:nrow(d)){
    #i = 1 
    file.copy(from = d$f[i], to = glue('random_inv/',d$id[i], '.jpg'))
    print(i)    
    }