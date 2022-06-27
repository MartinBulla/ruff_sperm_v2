require(data.table)
require(gtools)
require(kingship2)
require(MasterBayes)
require(MCMCglmm)
require(pedigreemm)
require(readxl)

Ruffs2020_inbreeding
pedigree2014to2019


# load data
  a[, animal := bird_ID]
  b[, animal := bird_ID]

  p = fread('Data/Dat_parentage.txt',na.strings = c(NA_character_, ""))
  setnames(p, old = c('bird_ID', 'Mother', 'Father'), new = c('id', 'dam', 'sire'))
  p = p[!is.na(id),.(id,sire,dam)]
  pc = tidyped(p)

  visped(pc) 

# pedigreemm
  p2 <- with(data.frame(pc), pedigreemm::pedigree(sire=Sire, dam=Dam, label=Ind))

  i = 'Midpiece'
  # on single values
   prior1 <- list(
      G = list(G1 = list(V = 1, nu = 0.002)),
      R = list(V = 1, nu = 0.002)
    )

  bb = b[part == i]
  bb[, id := bird_ID]
  m = lmer(scale(Length_µm) ~ Morph + (1|bird_ID), bb)
  bb[, res := resid(m)]

  mp1 = pedigreemm(res ~ 1 + (1|id), data = bb, pedigree = list(id = p2))
  summary(mp1)

# MCMCglmm
  Ainv = inverseA(data.frame(pc[,.(Ind,Dam,Sire)]))$Ainv
  i = 'Midpiece'
  # on single values
   prior1 <- list(
      G = list(G1 = list(V = 1, nu = 0.002)),
      R = list(V = 1, nu = 0.002)
    )

   bb = b[part == i]
   bb[, animal := bird_ID]
   m = lmer(scale(Length_µm) ~ Morph + (1|bird_ID), bb)
   bb[, res := resid(m)]

   m1 <- MCMCglmm(res ~ 1, random = ~animal, ginv = list(animal = Ainv), data = bb, prior = prior1,
        nitt = 65000, thin = 50, burnin = 15000, verbose = FALSE)
   plot(m1$Sol)
   plot(m1$VCV)
   autocorr.diag(m1$VCV)
   posterior.mode(m1$VCV)
   HPDinterval(m1$VCV)
   
   ph1 <- m1$VCV[, "animal"] / (m1$VCV[, "animal"] + m1$VCV[, "units"])
   posterior.mode(ph1)
   HPDinterval(ph1, 0.95)

   p.var <- var(bb$res, na.rm = TRUE)
   prior2 <- list(
      G = list(G1 = list(V = matrix(p.var * 0.95), nu = 1)),
      R = list(V = matrix(p.var * 0.05), nu = 1)
      )

   m2 <- MCMCglmm(res ~ 1, random = ~animal, ginv = list(animal = Ainv), data = bb, prior = prior2,
        nitt = 65000, thin = 50, burnin = 15000, verbose = FALSE)
   
   plot(m2$Sol)
   plot(m2$VCV)
   autocorr.diag(m2$VCV)
   posterior.mode(m2$VCV)
   HPDinterval(m2$VCV)
  
   ph2 <- m2$VCV[, "animal"] / (m2$VCV[, "animal"] + m2$VCV[, "units"])
   posterior.mode(ph2)
   HPDinterval(ph2, 0.95)
  

  # on averages 
   prior1 <- list(
      G = list(G1 = list(V = 1, nu = 0.002)),
      R = list(V = 1, nu = 0.002)
    )

   i ='Nucleus'
   aa = a[part == i]
   aa[, animal := bird_ID]
   m = lm(scale(Length_avg) ~ Morph, aa)
   aa[, res := resid(m)]
   
   m1 <- MCMCglmm(res ~ 1, random = ~animal, ginv = list(animal = Ainv), data = aa, prior = prior1,
        nitt = 65000, thin = 50, burnin = 15000, verbose = FALSE)
   plot(m1$Sol)
   plot(m1$VCV)
   autocorr.diag(m1$VCV)
   posterior.mode(m1$VCV)
   HPDinterval(m1$VCV)
   
   ph1 <- m1$VCV[, "animal"] / (m1$VCV[, "animal"] + m1$VCV[, "units"])
   posterior.mode(ph1)
   HPDinterval(ph1, 0.95)

   p.var <- var(aa$res, na.rm = TRUE)
   prior2 <- list(
      G = list(G1 = list(V = matrix(p.var * 0.95), nu = 1)),
      R = list(V = matrix(p.var * 0.05), nu = 1)
      )

   m2 <- MCMCglmm(res ~ 1, random = ~animal, ginv = list(animal = Ainv), data = aa, prior = prior2,
        nitt = 65000, thin = 50, burnin = 15000, verbose = FALSE)
   
   posterior.mode(m1$VCV)
   posterior.mode(m2$VCV)

   ph2 <- m2$VCV[, "animal"] / (m2$VCV[, "animal"] + m2$VCV[, "units"])
   posterior.mode(ph2)
   HPDinterval(ph2, 0.95)


# end