# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require('bayesplot')
  require(brms)
  require(gtools)
  require(kingship2)
  require(MasterBayes)
  require(MCMCglmm)
  require(nadiv)
  require(pedigreemm)
  require(readxl)
  require(visPedigree)

  # constants for MCMC
  cores_ = 4
  chains_ = 5
  iter_ = 40000
  thin_ = 20
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
    ddl[, animal := bird_ID]
    

    # use June values and for 4 males without June, May
      dd1 = dd[month == 'June']
      dd2 = dd[month == 'May']
      ddx = rbind(dd1,dd2[!bird_ID%in%dd1$bird_ID])
      ddxl = data.table(melt(ddx[,.(bird_ID,month,Morph,age,motileCount,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph","age","motileCount"), variable.name = "Motility"))
      ddxl[Motility == 'VAP' ,mot:='Average path']
      ddxl[Motility == 'VCL' ,mot:='Curvilinear']
      ddxl[Motility == 'VSL' ,mot:='Straight line']
      ddxl[, animal := bird_ID]
  
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
    #i = 'Acrosome'
    # on single values  
      bi = b[part == i]
      m = lmer(scale(Length_µm) ~ Morph + (1|bird_ID), bi)
      bi[, res := resid(m)]
    
      prior_no <- get_prior(res ~ 0 + Intercept + (1|bird_ID), data=bi)    
      mp_no = brm(res ~ 0 + Intercept + (1|bird_ID), data = bi, cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_no)
      
      prior_yes <- get_prior(res ~ 0 + Intercept + (1|bird_ID) + (1|gr(animal, cov = Amat)), data=bi, data2   = list(Amat = Amat)) 
      mp_yes = brm(res ~ 0 + Intercept  + (1|bird_ID) + (1 | gr(animal, cov = Amat)), data = bi,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_yes)
      
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
      mp_no = brm(res ~ 0 + Intercept, data = ai, cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_no)
      
      prior_yes <- get_prior(res ~ 0 + Intercept + (1|gr(animal, cov = Amat)), data=ai, data2   = list(Amat = Amat)) 
      mp_yes = brm(res ~ 0 + Intercept  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_yes)
      
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
  save(v,w, file = 'Outputs/temp_resPed_test_40000.Rdata')
  #load('Outputs/temp_resPed_test.Rdata')
  
  vs =list()
  va =list()
  


  for(i in unique(ddl$mot)){
    #i = 'Average path'
    # on single values  
      bi = ddl[mot == i]
      m = lmer(scale(value) ~ scale(log(motileCount)) + Morph + (1|bird_ID), bi)
      bi[, res := resid(m)]
    
      prior_no <- get_prior(res ~ 0 + Intercept + (1|bird_ID), data=bi)    
      mp_no = brm(res ~ 0 + Intercept + (1|bird_ID), data = bi, cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_no)
      
      prior_yes <- get_prior(res ~ 0 + Intercept + (1|bird_ID) + (1|gr(animal, cov = Amat)), data=bi, data2   = list(Amat = Amat)) 
      mp_yes = brm(res ~ 0 + Intercept  + (1|bird_ID) + (1 | gr(animal, cov = Amat)), data = bi,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_yes)
      
      save(mp_no, mp_yes, file = paste0('Data/sim/',i,'_res_sin.Rdata'))
      
      yy =bayes_factor(mp_no, mp_yes)
      zz = post_prob(mp_no, mp_yes)
      v_sc <- (VarCorr(mp_yes, summary = FALSE)$animal$sd)^2
      v_sp <- (VarCorr(mp_yes, summary = FALSE)$bird_ID$sd)^2
      v_r <- (VarCorr(mp_yes, summary = FALSE)$residual$sd)^2
      xx = summary(as.mcmc(v_sc / (v_sc + v_sp + v_r)))$quantiles
      vs[[i]] =  data.table(response = i, est = xx[3], lwr = xx[1], upr = xx[5], bf = yy$bf, prob=zz[1], Rhat_no = summary(mp_no)$spec_pars$Rhat, Rhat_yes =summary(mp_yes)$spec_pars$Rhat) #bf & prob in 

    # on averages values  
      ai = ddxl[mot == i]
      m = lm(scale(value) ~ scale(log(motileCount)) + Morph, ai)
      ai[, res := resid(m)]
    
      prior_no <- get_prior(res ~ 0 + Intercept, data=ai)    
      mp_no = brm(res ~ 0 + Intercept, data = ai, cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_no)
      
      prior_yes <- get_prior(res ~ 0 + Intercept + (1|gr(animal, cov = Amat)), data=ai, data2   = list(Amat = Amat)) 
      mp_yes = brm(res ~ 0 + Intercept  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_yes)
      
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

  save(vv,wv, file = 'Outputs/temp_resPed_test_mot_40000
    .Rdata')

arning messages:
1: There were 103 divergent transitions after warmup. See
https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
to find out why this is a problem and how to eliminate them.
2: There were 244 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
3: There were 5 chains where the estimated Bayesian Fraction of Missing Information was low. See
https://mc-stan.org/misc/warnings.html#bfmi-low
4: Examine the pairs() plot to diagnose sampling problems

5: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
Running the chains for more iterations may help. See
https://mc-stan.org/misc/warnings.html#bulk-ess
6: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
Running the chains for more iterations may help. See
https://mc-stan.org/misc/warnings.html#tail-ess
7: There were 103 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
8: There were 2 divergent transitions after warmup. See
https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
to find out why this is a problem and how to eliminate them.
9: Examine the pairs() plot to diagnose sampling problems

10: There were 2 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
11: There were 4 divergent transitions after warmup. See
https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
to find out why this is a problem and how to eliminate them.
12: There were 27 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
13: There were 5 chains where the estimated Bayesian Fraction of Missing Information was low. See
https://mc-stan.org/misc/warnings.html#bfmi-low
14: Examine the pairs() plot to diagnose sampling problems

15: There were 4 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
16: There were 70 divergent transitions after warmup. See
https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
to find out why this is a problem and how to eliminate them.
17: There were 5 chains where the estimated Bayesian Fraction of Missing Information was low. See
https://mc-stan.org/misc/warnings.html#bfmi-low
18: Examine the pairs() plot to diagnose sampling problems

19: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
Running the chains for more iterations may help. See
https://mc-stan.org/misc/warnings.html#tail-ess
20: There were 70 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
21: There were 55 divergent transitions after warmup. See
https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
to find out why this is a problem and how to eliminate them.
22: There were 5 chains where the estimated Bayesian Fraction of Missing Information was low. See
https://mc-stan.org/misc/warnings.html#bfmi-low
23: Examine the pairs() plot to diagnose sampling problems

24: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
Running the chains for more iterations may help. See
https://mc-stan.org/misc/warnings.html#tail-ess
25: There were 55 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
26: There were 8 divergent transitions after warmup. See
https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
to find out why this is a problem and how to eliminate them.
27: Examine the pairs() plot to diagnose sampling problems

28: There were 8 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
29: There were 3 divergent transitions after warmup. See
https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
to find out why this is a problem and how to eliminate them.
30: Examine the pairs() plot to diagnose sampling problems

31: There were 3 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

 load('Outputs/temp_resPed_test_mot_40000.Rdata')
 load('Outputs/temp_resPed_test_40000.Rdata')

wv[, data := 'June'] 
wv[, type := 'motility'] 
w[, data := 'means'] 
w[, type := 'morphology'] 
vv[, data := 'all']
vv[, type := 'motility'] 
v[, data := 'all']
v[, type := 'morphology'] 

x = rbind(v,vv,w,wv)
cols_=c('Rhat_no', 'Rhat_yes')
x[ , (cols_) := lapply(.SD, function(x){round(x,2)}), .SDcols = cols_]
x[, prob:=round(prob,3)]
x[, bf:=round(bf)]
cols_=c('est', 'lwr','upr')
x[ , (cols_) := lapply(.SD, function(x){round(x*100,1)}), .SDcols = cols_]
x[substring(upr,1,1) != 0, upr := round(upr)]
x[substring(lwr,1,1) != 0, lwr := round(lwr)]
x[, est := paste0(est, '%')]
x[, CI := paste0(lwr, '-', upr, '%')]


summary(x[data=='all',.(est, lwr, upr, bf, prob, Rhat_no, Rhat_yes)])
summary(x[data!='all',.(est, lwr, upr, bf, prob, Rhat_no, Rhat_yes)])

y = x[,.(type, data, response, est, CI, bf, prob, Rhat_no, Rhat_yes)]
y[, data := factor(data, levels=c("all","June","means"))] 
y[, type := factor(type, levels=c("motility","morphology"))] 
y[, response := factor(response, levels=c("Curvilinear","Average path",'Straight line',"Acrosome", "Nucleus",  "Midpiece","Tail","Total", "Head","Flagellum"))] 


fwrite(y[order(data, type,response)], file = 'Outputs/Table_Phylo.csv')


ggplot(b, aes(x = Length_µm))+facet_wrap(~part, scales = 'free') + geom_density()
ggplot(b, aes(x = Length_µm))+facet_wrap(~part, scales = 'free') + geom_histogram()
ggplot(a, aes(x = Length_avg))+facet_wrap(~part, scales = 'free') + geom_density()
ggplot(a, aes(x = Length_avg))+facet_wrap(~part, scales = 'free') + geom_histogram()

load('Data/sim/Acrosome_res_avg.RData')
model = mp_yes
posterior_cp <-  as.array(model)
lp_cp <- log_posterior(model)
np = nuts_params(model) 

summary(model)
plot(model, ask = FALSE)
pp_check(model, ndraws = 100) # same as function ppc_dens_overlay - see below
pp_check(model, ndraws = 4, type ='hist') # same as function ppc_hist - see below
pp_check(model, ndraws = 100, type ='stat') # same as function ppc_hist - see below
mcmc_plot(model, type = "acf")


mcmc_parcoord(posterior_cp, np = np)
mcmc_parcoord(lp_cp, np = np)
mcmc_parcoord(model, np = np)
mcmc_parcoord(model, transformations = 'log', np = np)

mcmc_plot(model, type = "neff")
mcmc_plot(model, type = "neff_hist")
mcmc_plot(model, type = "rhat")
mcmc_plot(model, type = "rhat_hist")
     
mcmc_plot(model, type = "nuts_acceptance")
mcmc_plot(model, type = "nuts_divergence")

mcmc_plot(model, type = "nuts_energy")
mcmc_plot(model, type = "trace_highlight")



#mcmc_pairs(posterior_cp, np = np)
mcmc_trace(posterior_cp, np = np, pars = c("b_Intercept",'sd_animal__Intercept','sigma')) + xlab("Post-warmup iteration")
pairs(model, np = np)
mcmc_scatter(posterior_cp, pars = c("b_Intercept", "sigma"),  size = 1,
  transform = list(sigma = "log"), 
  np = np) #mcmc_plot(model, type = 'scatter', pars=c("b_Intercept",'sd_animal__Intercept'))
mcmc_scatter(posterior_cp, pars = c("sd_animal__Intercept", "sigma"),  size = 1,
  #transform = list(sigma = "log"), 
  np = np)


# posterior predictive checks - easier to use pp_check 
  res = ai$res
  res_rep = posterior_predict(model, draws = 100)
 
  #ppc_stat() 
  ppc_dens_overlay(y = as.numeric(res), res_rep)     

     



  type: The type of the plot. Supported types are (as names) ‘hist’,
          ‘dens’, ‘hist_by_chain’, ‘dens_overlay’, ‘violin’,
          ‘intervals’, ‘areas’, ‘acf’, ‘acf_bar’,‘trace’,
          ‘trace_highlight’, ‘scatter’, ‘rhat’, ‘rhat_hist’, ‘neff’,
          ‘neff_hist’ ‘nuts_acceptance’, ‘nuts_divergence’,
          ‘nuts_stepsize’, ‘nuts_treedepth’, and ‘nuts_energy’. For an
          overview on the various plot types see ‘MCMC-overview’.

# end