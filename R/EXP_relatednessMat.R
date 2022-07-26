# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require('bayesplot')
  require(brms)
  require(gtools)
  require(MasterBayes)
  require(MCMCglmm)
  require(nadiv)
  require(pedigreemm)
  require(readxl)
  require(visPedigree)

  require(admisc)
  require(gap)
  require(lme4qtl)
  require(related)
  require(reshape2)
  require(matrixcalc)

  # constants for MCMC
    cores_ = 20
    chains_ = 4
    iter_ = 50000
    thin_ = 20
    sample_ = chains_*(iter_/2)/thin_
# load data
  load(file = 'Data/DAT_rel-mat.RData') # m - created using DAT_pedigree.R
  m0 = pmax(m,0) # in case it is needed, change negative estimates to zero
  Amat = m0
  diag(Amat) = diag(Amat)+0.1 # making matrix positive def https://discourse.mc-stan.org/t/positive-definite-vs-positive-semidefinite-required-for-phylogenetic-regression-error-message-in-brms/17049

  source(here::here('R/DAT_prepare.R'))      

  a[, animal := bird_ID]
  a[, id := bird_ID]
  b[, animal := bird_ID]
  b[, id := bird_ID]
  cv_[, animal := bird_ID]
  cv_[, id := bird_ID]
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
  
# brms
  # morpho
    ls =list()
    la =list()
    for(i in unique(b$part)){
      #i = 'Acrosome'
      # on single values  
        bi = b[part == i]
        m = lmer(scale(Length_µm) ~ Morph + (1|bird_ID), bi)
        bi[, res := resid(m)]
      
        #prior_no <- get_prior(res ~ 0 + Intercept + (1|bird_ID), data=bi)    
        mp_no = brm(res ~ 0 + Intercept + (1|bird_ID), data = bi, cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))#, prior   = prior_no)
        
        #prior_yes <- get_prior(res ~ 0 + Intercept + (1|bird_ID) + (1|gr(animal, cov = Amat)), data=bi, data2   = list(Amat = Amat)) 
        mp_yes = brm(res ~ 0 + Intercept  + (1|bird_ID) + (1 | gr(animal, cov = Amat)), data = bi,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))#, prior   = prior_yes)
        
        save(mp_no, mp_yes, file = paste0('Data/sim/',i,'_res_sin_relatedness_',sample_,'.Rdata'))
        
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
      
        #prior_no <- get_prior(res ~ 0 + Intercept, data=ai)    
        mp_no = brm(res ~ 0 + Intercept, data = ai, cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))#, prior   = prior_no)
        
        #prior_yes <- get_prior(res ~ 0 + Intercept + (1|gr(animal, cov = Amat)), data=ai, data2   = list(Amat = Amat)) 
        mp_yes_ped = brm(res ~ 0 + Intercept  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))#, prior   = prior_yes)
        
        save(mp_no, mp_yes, file = paste0('Data/sim/',i,'_res_avg_relatedness_',sample_,'.Rdata'))
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
    save(v,w, file = paste0('Outputs/resRel_morpho_',sample_,'.Rdata'))
    #load('Outputs/temp_resPed_test.Rdata')
  # motil
    vs =list()
    va =list()

    for(i in unique(ddl$mot)){
      #i = 'Average path'
      # on single values  
        bi = ddl[mot == i]
        m = lmer(scale(value) ~ scale(log(motileCount)) + Morph + (1|bird_ID), bi)
        bi[, res := resid(m)]
      
        #prior_no <- get_prior(res ~ 0 + Intercept + (1|bird_ID), data=bi)    
        mp_no = brm(res ~ 0 + Intercept + (1|bird_ID), data = bi, cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))#, prior   = prior_no)
        
        #prior_yes <- get_prior(res ~ 0 + Intercept + (1|bird_ID) + (1|gr(animal, cov = Amat)), data=bi, data2   = list(Amat = Amat)) 
        mp_yes = brm(res ~ 0 + Intercept  + (1|bird_ID) + (1 | gr(animal, cov = Amat)), data = bi,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))#, prior   = prior_yes)
        
        save(mp_no, mp_yes, file = paste0('Data/sim/',i,'_res_sin_relatedness_',sample_,'.Rdata'))
        
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
      
        #prior_no <- get_prior(res ~ 0 + Intercept, data=ai)    
        mp_no = brm(res ~ 0 + Intercept, data = ai, cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))#, prior   = prior_no)
        
        #prior_yes <- get_prior(res ~ 0 + Intercept + (1|gr(animal, cov = Amat)), data=ai, data2   = list(Amat = Amat)) 
        mp_yes = brm(res ~ 0 + Intercept  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))#, prior   = prior_yes)
        
        save(mp_no, mp_yes, file = paste0('Data/sim/',i,'_res_avg_relatedness_',sample_,'.Rdata'))
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
    save(vv,wv, file = paste0('Outputs/resRel_motil_',sample_,'.Rdata'))
  # CV
    vcv_ =list()
    for(i in unique(b$part)){
      #i = 'Tail'
        ai = cv_[part == i]
        m = lm(scale(CV) ~ Morph, ai)
        ai[, res := resid(m)]

        adapt_d = 0.999
      
        #prior_no <- get_prior(res ~ 0 + Intercept, data=ai)    
        mp_no = brm(res ~ 0 + Intercept, data = ai, cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))#, prior   = prior_no)
        
        #prior_yes <- get_prior(res ~ 0 + Intercept + (1|gr(animal, cov = Amat)), data=ai, data2   = list(Amat = Amat)) 
        mp_yes = brm(res ~ 0 + Intercept  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = adapt_d), sample_prior="yes",save_pars = save_pars(all = TRUE))#, prior   = prior_yes)
        
        save(mp_no, mp_yes, file = paste0('Data/sim/',i,'_res_CV_relatedness',sample_,'.Rdata'))
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
        vcv_[[i]] =  data.table(response = i, est = xx[3], lwr = xx[1], upr = xx[5], bf = yy$bf, prob=zz[1], Rhat_no = summary(mp_no)$spec_pars$Rhat, Rhat_yes =summary(mp_yes)$spec_pars$Rhat ) #bf & prob in favor of NO pedigree model
        #bsim = sim(mp, n.sim=nsim) 
          
        print(i)
    }
   
    wcv = do.call(rbind,vcv_)
    save(wcv, file = paste0('Outputs/resRel_CV_',sample_,'.Rdata'))
    
# Export summary
 #load('Outputs/temp_resRel_morpho.Rdata')
 #load('Outputs/temp_resRel_motil.Rdata')
 #load('Outputs/temp_resRel_CV.Rdata')

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
  x[ , (cols_) := lapply(.SD, function(x){round(x*100,1)}), .SDcols = cols_]
  x[substring(upr,1,1) != 0, upr := round(upr)]
  x[substring(lwr,1,1) != 0, lwr := round(lwr)]
  x[, est := paste0(est, '%')]
  x[, CI := paste0(lwr, '-', upr, '%')]


  summary(x[data=='all',.(est, lwr, upr, bf, prob, Rhat_no, Rhat_yes)])
  summary(x[data!='all',.(est, lwr, upr, bf, prob, Rhat_no, Rhat_yes)])

  y = x[,.(type, data, response, est, CI, bf, prob, Rhat_no, Rhat_yes)]
  y[, data := factor(data, levels=c("all","June","means", 'CV'))] 
  y[, type := factor(type, levels=c("motility","morphology"))] 
  y[, response := factor(response, levels=c("Curvilinear","Average path",'Straight line',"Acrosome", "Nucleus",  "Midpiece","Tail","Total", "Head","Flagellum"))] 


  fwrite(y[order(data, type,response)], file = paste0('Outputs/Table_Srelate_',sample_,'.csv'))

# explore chains
  ggplot(b, aes(x = Length_µm))+facet_wrap(~part, scales = 'free') + geom_density()
  ggplot(b, aes(x = Length_µm))+facet_wrap(~part, scales = 'free') + geom_histogram()
  ggplot(a, aes(x = Length_avg))+facet_wrap(~part, scales = 'free') + geom_density()
  ggplot(a, aes(x = Length_avg))+facet_wrap(~part, scales = 'free') + geom_histogram()

# yes CV Acrosome 2 div - is ok
# yes CV Nucleus 1 div - is ok
# yes CV Tail 12 div - looks okaisch

load('Data/sim/Nucleus_res_CV_relatedness.RData')
load('Data/sim/Straight line_res_avg_relatedness.RData')
load('Data/sim/Total_res_avg_relatedness.RData')
load('Data/sim/Acrosome_res_sin_relatedness.RData')
load('Data/sim/Curvilinear_res_sin_relatedness.RData')
model = mp_yes
posterior_cp <-  as.array(model)
lp_cp <- log_posterior(model)
np = nuts_params(model) 
summary(model)


f = c(list.files(path = here::here('Data/sim/'), pattern = 'relatedness', recursive = TRUE, full.names = TRUE))
  for(i in f){
    #i=f[1]
    load(i)
    print(paste('no',i))
    print(tryCatchWEM(summary(mp_no)))
    print(paste('yes',i))
    print(tryCatchWEM(summary(mp_yes)))
  }  


### START HEREE
# 1 flagellum
# 1 nucleus
# Increasing adapt_delta above 0.999 & explore whether ok or not - as of now looks okeisch

"yes /ds/grpkempenaers/Martin/ruff_sperm_v2/Data/sim//Flagellum_res_CV_relatedness5000.Rdata"
$warning
[1] "There were 13 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup"

"yes /ds/grpkempenaers/Martin/ruff_sperm_v2/Data/sim//Midpiece_res_CV_relatedness5000.Rdata"
$warning
[1] "There were 4 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup"

load(here::here('Data/sim/Flagellum_res_CV_relatedness5000.Rdata'))

model = mp_yes
posterior_cp <-  as.array(model)
lp_cp <- log_posterior(model)
np = nuts_params(model) 

plot(model, ask = FALSE)
pp_check(model, ndraws = 100) # same as function ppc_dens_overlay - see below
pp_check(model, ndraws = 4, type ='hist') # same as function ppc_hist - see below
pp_check(model, ndraws = 100, type ='stat') # same as function ppc_hist - see below
mcmc_plot(model, type = "acf")


mcmc_trace(posterior_cp, np = 
  np, pars = c("b_Intercept",'sd_animal__Intercept','sigma')) + xlab("Post-warmup iteration")
pairs(model, np = np)
mcmc_scatter(posterior_cp, pars = c("b_Intercept", "sigma"),  size = 1,
  transform = list(sigma = "log"), 
  np = np) #mcmc_plot(model, type = 'scatter', pars=c("b_Intercept",'sd_animal__Intercept'))
mcmc_scatter(posterior_cp, pars = c("sd_animal__Intercept", "sigma"),  size = 1,
  #transform = list(sigma = "log"), 
  np = np)



mcmc_parcoord(lp_cp, np = np)
mcmc_parcoord(model, np = np) 
mcmc_parcoord(posterior_cp, np = np)
mcmc_parcoord(model, transformations = 'log', np = np)

mcmc_plot(model, type = "neff")
mcmc_plot(model, type = "neff_hist")
mcmc_plot(model, type = "rhat")
mcmc_plot(model, type = "rhat_hist")
     
mcmc_plot(model, type = "nuts_acceptance")
mcmc_plot(model, type = "nuts_divergence")

mcmc_plot(model, type = "nuts_energy")
#mcmc_plot(model, type = "trace_highlight")


#mcmc_pairs(posterior_cp, np = np)

# posterior predictive checks - easier to use pp_check 
  #res = ai$res
  #res_rep = posterior_predict(model, draws = 100)
 
  #ppc_stat() 
  #ppc_dens_overlay(y = as.numeric(res), res_rep)     

     



  type: The type of the plot. Supported types are (as names) ‘hist’,
          ‘dens’, ‘hist_by_chain’, ‘dens_overlay’, ‘violin’,
          ‘intervals’, ‘areas’, ‘acf’, ‘acf_bar’,‘trace’,
          ‘trace_highlight’, ‘scatter’, ‘rhat’, ‘rhat_hist’, ‘neff’,
          ‘neff_hist’ ‘nuts_acceptance’, ‘nuts_divergence’,
          ‘nuts_stepsize’, ‘nuts_treedepth’, and ‘nuts_energy’. For an
          overview on the various plot types see ‘MCMC-overview’.

# end