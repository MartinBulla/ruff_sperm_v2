# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(brms)
  require(gtools)
  require(kingship2)
  require(MasterBayes)
  require(MCMCglmm)
  require(nadiv)
  require(pedigreemm)
  require(readxl)
  require(visPedigree)


# load data
  source(here::here('R/DAT_prepare.R'))      

  a[, animal := bird_ID]
  a[, id := bird_ID]
  b[, animal := bird_ID]
  b[, id := bird_ID]

  # prepare motility
    dd = d[!Morph%in%'Zebra finch']
    ddl = data.table(melt(dd[,.(bird_ID,month,Morph,age,motileCount,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph","age","motileCount"), variable.name = "Motility"))
    ddl[Motility == 'VAP' ,mot:='Average path']
    ddl[Motility == 'VCL' ,mot:='Curvilinear']
    ddl[Motility == 'VSL' ,mot:='Straight line']
    

    # use June values and for 4 males without June, May
      dd1 = dd[month == 'June']
      dd2 = dd[month == 'May']
      ddx = rbind(dd1,dd2[!bird_ID%in%dd1$bird_ID])
      ddxl = data.table(melt(ddx[,.(bird_ID,month,Morph,age,motileCount,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph","age","motileCount"), variable.name = "Motility"))
      ddxl[Motility == 'VAP' ,mot:='Average path']
      ddxl[Motility == 'VCL' ,mot:='Curvilinear']
      ddxl[Motility == 'VSL' ,mot:='Straight line']
  
  # pedigree
    p = fread('Data/Dat_parentage.txt',na.strings = c(NA_character_, ""))
    setnames(p, old = c('bird_ID', 'Mother', 'Father'), new = c('id', 'dam', 'sire'))
    p = p[!is.na(id),.(id,sire,dam)]
    pc = tidyped(p)
    #visped(pc) 
    p2 <- with(data.frame(pc), pedigreemm::pedigree(sire=Sire, dam=Dam, label=Ind))
# pedigreemm

  ls =list()
  la =list()
  for(i in unique(b$part)){
    #i = 'Midpiece'
    # on single values  
      bi = b[part == i]
      m = lmer(scale(Length_µm) ~ Morph + (1|bird_ID), bi)
      bi[, res := resid(m)]
      mp = pedigreemm(res ~ 1 + (1|id), data = bi, pedigree = list(id = p2))
      #bsim = sim(mp, n.sim=nsim) 
      ls[[i]] = m_out(name = "Table S test", dep = i,model = mp, nsim = 5000, save_sim = FALSE)
    
    # on averages values  
      ai = a[part == i]
      m = lm(scale(Length_avg) ~ Morph, ai)
      ai[, res := resid(m)]
      mp2 = pedigreemm(res ~ 1 + (1|id), data = ai, pedigree = list(id = p2))
      #bsim = sim(mp, n.sim=nsim) 
      ls[[i]] = m_out(name = "Table S test", dep = i,model = mp2, nsim = 5000, save_sim = FALSE)
        

      print(i)
  }

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

# brms
  pc_ = pc[,.(Ind,Dam,Sire)]
  setnames(pc_, old = 'Ind', new = 'ID')
  Amat <- as.matrix(nadiv::makeA(pc_))

  ls =list()
  la =list()
  for(i in unique(b$part)){
    #i = 'Midpiece'
    # on single values  
      bi = b[part == i]
      m = lmer(scale(Length_µm) ~ Morph + (1|bird_ID), bi)
      bi[, res := resid(m)]
    
      prior_no <- get_prior(res ~ 0 + Intercept + (1|bird_ID), data=bi)    
      mp_no = brm(res ~ 0 + Intercept + (1|bird_ID), data = bi, cores = 2, chains = 4, iter = 5000, thin = 5, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_no)
      
      prior_yes <- get_prior(res ~ 0 + Intercept + (1|bird_ID) + (1|gr(animal, cov = Amat)), data=bi, data2   = list(Amat = Amat)) 
      mp_yes = brm(res ~ 0 + Intercept  + (1|bird_ID) + (1 | gr(animal, cov = Amat)), data = bi,  data2 = list(Amat = Amat), cores = 2, chains = 4, iter = 5000, thin = 5, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = priors)
      
      save(mp_no, mp_yes, file = paste0('Data/sim/',i,'_res_sin.Rdata'))
      
      yy =bayes_factor(mp_no, mp_yes)
      zz = post_prob(mp_no, mp_yes)
      v_sc <- (VarCorr(mp_yes, summary = FALSE)$animal$sd)^2
      v_sp <- (VarCorr(mp_yes, summary = FALSE)$bird_ID$sd)^2
      v_r <- (VarCorr(mp_yes, summary = FALSE)$residual$sd)^2
      xx = summary(as.mcmc(v_sc / (v_sc + v_sp + v_r)))$quantiles
      ls[[i]] =  data.table(response = i, est = xx[3], lwr = xx[1], upr = xx[5], bf = yy$bf, prob=zz[1], Rhat_no = summary(mp_no)$spec_pars$Rhat, Rhat_yes =summary(mp_yes)$spec_pars$Rhat) #bf & prob in 

    # on averages values  
      ai = a[part == i]
      m = lm(scale(Length_avg) ~ Morph, ai)
      ai[, res := resid(m)]
    
      prior_no <- get_prior(res ~ 0 + Intercept, data=ai)    
      mp_no = brm(res ~ 0 + Intercept, data = ai, cores = 2, chains = 4, iter = 5000, thin = 5, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_no)
      
      prior_yes <- get_prior(res ~ 0 + Intercept + (1|gr(animal, cov = Amat)), data=ai, data2   = list(Amat = Amat)) 
      mp_yes = brm(res ~ 0 + Intercept  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat), cores = 2, chains = 4, iter = 5000, thin = 5, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = priors)
      
      save(mp_no, mp_yes, file = paste0('Data/sim/',i,'_res_avg.Rdata'))
      #summary(mp_yes)
      #plot(mp_yes)
      #mcmc_plot(mp_yes, type = "acf")

      #hypothesis(mp_yes, "sd_animal__Intercept^2 / (sd_animal__Intercept^2 + sigma^2) = 0", class = NULL)

      yy =bayes_factor(mp_no, mp_yes)
      zz = post_prob(mp_no, mp_yes)  
      v_sc <- (VarCorr(mp_yes, summary = FALSE)$animal$sd)^2
      #v_sp <- (VarCorr(mp2, summary = FALSE)$Species$sd)^2
      v_r <- (VarCorr(mp_yes, summary = FALSE)$residual$sd)^2
      xx = summary(as.mcmc(v_sc / (v_sc + v_r)))$quantiles
      la[[i]] =  data.table(response = i, est = xx[3], lwr = xx[1], upr = xx[5], bf = yy$bf, prob=zz[1], Rhat_no = summary(mp_no)$spec_pars$Rhat, Rhat_yes =summary(mp_yes)$spec_pars$Rhat ) #bf & prob in favor of NO pedigree model
      #bsim = sim(mp, n.sim=nsim) 
        
      print(i)
  }

  v = do.call(rbind,ls)
  w = do.call(rbind,la)

  vs =list()
  va =list()
  for(i in unique(ddl$mot)){
    #i = 'Average path'
    # on single values  
      bi = ddl[mot == i]
      m = lmer(scale(value) ~ scale(log(motileCount)) + Morph + (1|bird_ID), bi)
      bi[, res := resid(m)]
    
      prior_no <- get_prior(res ~ 0 + Intercept + (1|bird_ID), data=bi)    
      mp_no = brm(res ~ 0 + Intercept + (1|bird_ID), data = bi, cores = 2, chains = 4, iter = 5000, thin = 5, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_no)
      
      prior_yes <- get_prior(res ~ 0 + Intercept + (1|bird_ID) + (1|gr(animal, cov = Amat)), data=bi, data2   = list(Amat = Amat)) 
      mp_yes = brm(res ~ 0 + Intercept  + (1|bird_ID) + (1 | gr(animal, cov = Amat)), data = bi,  data2 = list(Amat = Amat), cores = 2, chains = 4, iter = 5000, thin = 5, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = priors)
      
      save(mp_no, mp_yes, file = paste0('Data/sim/',i,'_res_sin.Rdata'))
      
      yy =bayes_factor(mp_no, mp_yes)
      zz = post_prob(mp_no, mp_yes)
      v_sc <- (VarCorr(mp_yes, summary = FALSE)$animal$sd)^2
      v_sp <- (VarCorr(mp_yes, summary = FALSE)$bird_ID$sd)^2
      v_r <- (VarCorr(mp_yes, summary = FALSE)$residual$sd)^2
      xx = summary(as.mcmc(v_sc / (v_sc + v_sp + v_r)))$quantiles
      vs[[i]] =  data.table(response = i, est = xx[3], lwr = xx[1], upr = xx[5], bf = yy$bf, prob=zz[1], Rhat_no = summary(mp_no)$spec_pars$Rhat, Rhat_yes =summary(mp_yes)$spec_pars$Rhat) #bf & prob in 

    # on averages values  
      ai = ddxl[part == i]
      m = lm(scale(value) ~ scale(log(motileCount)) + Morph, ai)
      ai[, res := resid(m)]
    
      prior_no <- get_prior(res ~ 0 + Intercept, data=ai)    
      mp_no = brm(res ~ 0 + Intercept, data = ai, cores = 2, chains = 4, iter = 5000, thin = 5, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_no)
      
      prior_yes <- get_prior(res ~ 0 + Intercept + (1|gr(animal, cov = Amat)), data=ai, data2   = list(Amat = Amat)) 
      mp_yes = brm(res ~ 0 + Intercept  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat), cores = 2, chains = 4, iter = 5000, thin = 5, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = priors)
      
      save(mp_no, mp_yes, file = paste0('Data/sim/',i,'_res_avg.Rdata'))
      #summary(mp_yes)
      #plot(mp_yes)
      #mcmc_plot(mp_yes, type = "acf")

      #hypothesis(mp_yes, "sd_animal__Intercept^2 / (sd_animal__Intercept^2 + sigma^2) = 0", class = NULL)

      yy =bayes_factor(mp_no, mp_yes)
      zz = post_prob(mp_no, mp_yes)  
      v_sc <- (VarCorr(mp_yes, summary = FALSE)$animal$sd)^2
      #v_sp <- (VarCorr(mp2, summary = FALSE)$Species$sd)^2
      v_r <- (VarCorr(mp_yes, summary = FALSE)$residual$sd)^2
      xx = summary(as.mcmc(v_sc / (v_sc + v_r)))$quantiles
      va[[i]] =  data.table(response = i, est = xx[3], lwr = xx[1], upr = xx[5], bf = yy$bf, prob=zz[1], Rhat_no = summary(mp_no)$spec_pars$Rhat, Rhat_yes =summary(mp_yes)$spec_pars$Rhat ) #bf & prob in favor of NO pedigree model
      #bsim = sim(mp, n.sim=nsim) 
        
      print(i)
  }

  vv = do.call(rbind,vs)
  wv = do.call(rbind,va)


  save(v,w, file = 'Outputs/temp_resPed_test.Rdata')


# end