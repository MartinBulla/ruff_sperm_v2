# =============================================================
# ❗ Runs relative to the project's root directory,
# loads selected pictures for measurements in Sperm Sizer from
# "Data/sperm_pics/measured", randomizes the pictures and 
# copies the renamed pictures into "Data/sperm_pics/random_inv/"
# and "randomised_pics.txt" containing the links between blinded
# and original pictures to ./DATA/
# =============================================================

# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))

# LOAD picture names
  m <- data.table(
    f = c(list.files(path = here::here("Data/sperm_pics/measured"), pattern = ".jpg", recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("Data/sperm_pics/measured"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  
  m[, sample_ID := substring(f2, 6, 8)]
  
  nrow(m) # N = 920
  
  # check N sample/bird  
    mm <- m[, .N, by = sample_ID]
    mm[N < 10]
    length(unique(m$sample_ID)) # N birds = 92
    #summary(factor(m$sample_ID))

# RANDOMISE
 m[, file_name := sub(".*/", "", f2)]
 m[, bird_ID := substring(f2, 10, nchar(f2) - 29)]
 m[bird_ID == "", bird_ID := substring(f2, 10, nchar(f2) - 25)]
 # substring('Ruff 018_1307 2021-11-10 Snap-3277.jpg', 10,nchar('Ruff 018_1307 2021-11-10 Snap-3277.jpg')-25)
 m[bird_ID == "G200067 2", bird_ID := "G200067"]
 length(unique(m$bird_ID)) # 92
 # m[nchar(file_name)<32]

  # TEST
    dd <- m[, .N, by = bird_ID]
    dd[N < 10]
    dd[N > 10]
    #m[bird_ID %in% dd$bird_ID[dd$N > 10], .N, by = list(bird_ID, sample_ID)]

 m[, id := as.character(sample.int(n = nrow(m), size = nrow(m)))]
 m[nchar(id) == 1, id := paste0("00", id)]
 m[nchar(id) == 2, id := paste0("0", id)]
 m <- m[order(id)]
 fwrite(m, file = "Data/DAT_randomised_pics.csv")

# COPY randomised
 for (i in 1:nrow(m)) {
   # i = 1
   file.copy(from = m$f[i], to = glue("Data/sperm_pics/random_inv/", m$id[i], ".jpg"))
   print(i)
 }

# CHECKS
 
  # N randomized
  r <- data.table(
    f = c(list.files(path = here::here("Data/sperm_pics/random_inv"), pattern = ".jpg", recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("Data/sperm_pics/random_inv"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  nrow(r)

  # N of all photographs
  d = data.table(
    f = c(list.files(path = here::here("Data/sperm_pics/original"), pattern = ".jpg", recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("Data/sperm_pics/original"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  d[, sample_ID := substring(f2, 6, 8)]
  nrow(d) # 2065
  length(unique(d$sample_ID)) # unique(d$sample_ID)
  summary(factor(d$sample_ID))
  dd = d[, .N, by = sample_ID]
  dd[N < 10]

  # N of all inverted
  n <- data.table(
    f = c(list.files(path = here::here("Data/sperm_pics/inverted"), pattern = ".jpg", recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here("Data/sperm_pics/inverted"), pattern = ".jpg", recursive = TRUE, full.names = FALSE))
  )
  nrow(n) # 2066
  
# END