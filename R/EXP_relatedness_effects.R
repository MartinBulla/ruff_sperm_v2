# TOOLS
require(here)
source(here::here("R/tools.R"))
require(admisc)
require(bayesplot)
require(brms)
require(gtools)
require(MasterBayes)

# EXPLORE DIVERGENCE - Increasing adapt_delta above 0.999 makes all look okeisch
f <- c(list.files(path = here::here("Data/sim/"), pattern = "relatedness", recursive = TRUE, full.names = TRUE))
for (i in f) {
  # i=f[1]
  load(i)
  print(paste("no", i))
  print(tryCatchWEM(summary(mp_no)))
  print(paste("yes", i))
  print(tryCatchWEM(summary(mp_yes)))
}

# EXPLORE CHAINS 
load(here::here("Data/sim/Flagellum_res_CV_relatedness5000.Rdata")) # choose output to load and explore
#load("Data/sim/Nucleus_res_CV_relatedness.RData")
#load("Data/sim/Straight line_res_avg_relatedness.RData")
#load("Data/sim/Total_res_avg_relatedness.RData")
#load("Data/sim/Acrosome_res_sin_relatedness.RData")
#load("Data/sim/Curvilinear_res_sin_relatedness.RData")

model <- mp_yes # choose model to use and explore mp_yes (controlled for relatedness) or mp_no (not)
posterior_cp <- as.array(model)
lp_cp <- log_posterior(model)
np <- nuts_params(model)
summary(model)

plot(model, ask = FALSE)
pp_check(model, ndraws = 100) # same as function ppc_dens_overlay - see below
pp_check(model, ndraws = 4, type ='hist') # same as function ppc_hist - see below
pp_check(model, ndraws = 100, type ='stat') # same as function ppc_hist - see below
mcmc_plot(model, type = "acf")

mcmc_trace(posterior_cp, np = np, pars = c("b_Intercept",'sd_animal__Intercept','sigma')) + xlab("Post-warmup iteration")
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


# INFO
  # yes CV Acrosome 2 div - is ok
  # yes CV Nucleus 1 div - is ok
  # yes CV Tail 12 div - looks okaisch

  # "yes /ds/grpkempenaers/Martin/ruff_sperm_v2/Data/sim//Flagellum_res_CV_relatedness5000.Rdata"$warning
  # [1] "There were 13 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup"

  # "yes /ds/grpkempenaers/Martin/ruff_sperm_v2/Data/sim//Midpiece_res_CV_relatedness5000.Rdata"
  # $warning
  # [1] "There were 4 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup"

# END