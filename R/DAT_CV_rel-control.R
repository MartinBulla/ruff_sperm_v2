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
  require(standist)

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
  
# relatedness controlled models for CV - cauchy prior
  effects_ = c('Satellite relative to independent','Faeder relative to independent', 'Faeder relative to satellite')
  # constants for MCMC
    cores_ = 60
    chains_ = 4
    iter_ = 50000
    thin_ = 20
    sample_ = chains_*(iter_/2)/thin_
    adapt_d = 0.999

    #cores_ = 1#20
    #chains_ = 1
    #iter_ = 4000
    #thin_ = 2
    #sample_ = chains_*(iter_/2)/thin_
    #adapt_d = 0.999
  
  prior_x = c(
    prior(normal(0, 5), "Intercept"), #prior(normal(0, 10), class = Intercept),
    prior(normal(0, 5), class = b, coef = MorphFaeder), #prior(normal(0, 10), class = b, coef = gender),
    prior(normal(0, 5), class = b, coef = MorphSatellite),
    prior(cauchy(0, 5), class = sd),
    prior(cauchy(0, 5), "sigma")
  )

  #visualize("normal(0, 3)", "normal(0, 5)", "normal(0, 10)",'normal(1,2)','student_t(3, 0, 20)', xlim = c(-20,20))
  #visualize("cauchy(0, 5)",xlim = c(0,20))


  lmi =list()
  lco =list()
  for(i in unique(b$part)){
    #i = 'Tail'
    ai = cv_[part == i]
    
    #get_prior(scale(CV)  ~ Morph  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat))

    mib = brm(scale(CV)  ~ Morph  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = adapt_d, max_treedepth = 15), sample_prior="yes",save_pars = save_pars(all = TRUE), prior   = prior_x)

    mi_d1 = data.table(summary(mi)$fixed)
    mi_d1[, param:= rownames(summary(mi)$fixed)]
    mi_d1[,type :='fixed']

    mi_d2 = data.table(summary(mi)$random$animal)
    mi_d2[, param := 'relatedness_sd']
    mi_d2[, type :='random']

    mi_d3 = data.table(summary(mi)$spec_pars)
    mi_d3[, param := 'residual_sd']
    mi_d3[, type :='random']

    mi_d=rbind(mi_d1,mi_d2,mi_d3)
    mi_d[,response:=paste('CV', i)]
    lmi[[i]] = mi_d


    mi_co = data.table(
            rbind(hypothesis(mi, "MorphSatellite = 0")$hypothesis,
                  hypothesis(mi, "MorphFaeder = 0")$hypothesis,
                  hypothesis(mi, "MorphFaeder - MorphSatellite = 0")$hypothesis
                  ))
    mi_co[, effect := effects_]
     mi_co[,response:=paste('CV', i)]
    lco[[i]] = mi_co

    save(mi, file = paste0('Data/sim/CV_',i,'_relatedness-control_',sample_,'.Rdata'))
    print(i)
  }
  mi_ = do.call(rbind,lmi)
  mi_co_ = do.call(rbind,lco)
  save(mi_,mi_co_,file = paste0('Outputs/CV_rel_control_cauchy_', sample_,'.Rdata'))
  

  # check for warnings - no divergent transitions
    load(file = 'Outputs/CV_rel_control_cauchy_5000.Rdata')
    f = c(list.files(path = here::here('Data/sim/'), pattern = 'relatedness-control_5000', recursive = TRUE, full.names = TRUE))
    for(i in f){
      #i=f[1]
      load(i)
      print(paste('no',i))
      print(tryCatchWEM(summary(mp_no)))
      print(paste('yes',i))
      print(tryCatchWEM(summary(mp_yes)))
    }  
# relatedness controlled models for CV - default
  effects_ = c('Satellite relative to independent','Faeder relative to independent', 'Faeder relative to satellite')
  # constants for MCMC
    cores_ = 60
    chains_ = 4
    iter_ = 50000
    thin_ = 20
    sample_ = chains_*(iter_/2)/thin_
    adapt_d = 0.999

    #cores_ = 1#20
    #chains_ = 1
    #iter_ = 4000
    #thin_ = 2
    #sample_ = chains_*(iter_/2)/thin_
    #adapt_d = 0.999
  

  lmi =list()
  lco =list()
  for(i in unique(b$part)){
    #i = 'Tail'
    ai = cv_[part == i]
    
    #get_prior(scale(CV)  ~ Morph  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat))

    mib = brm(scale(CV)  ~ Morph  + (1 | gr(animal, cov = Amat)), data = ai,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = adapt_d, max_treedepth = 15), sample_prior="yes",save_pars = save_pars(all = TRUE))

    mi_d1 = data.table(summary(mi)$fixed)
    mi_d1[, param:= rownames(summary(mi)$fixed)]
    mi_d1[,type :='fixed']

    mi_d2 = data.table(summary(mi)$random$animal)
    mi_d2[, param := 'relatedness_sd']
    mi_d2[, type :='random']

    mi_d3 = data.table(summary(mi)$spec_pars)
    mi_d3[, param := 'residual_sd']
    mi_d3[, type :='random']

    mi_d=rbind(mi_d1,mi_d2,mi_d3)
    mi_d[,response:=paste('CV', i)]
    lmi[[i]] = mi_d


    mi_co = data.table(
            rbind(hypothesis(mi, "MorphSatellite = 0")$hypothesis,
                  hypothesis(mi, "MorphFaeder = 0")$hypothesis,
                  hypothesis(mi, "MorphFaeder - MorphSatellite = 0")$hypothesis
                  ))
    mi_co[, effect := effects_]
     mi_co[,response:=paste('CV', i)]
    lco[[i]] = mi_co

    save(mi, file = paste0('Data/sim/CV_',i,'_relatedness-control_default_',sample_,'.Rdata'))
    print(i)
  }
  mi_ = do.call(rbind,lmi)
  mi_co_ = do.call(rbind,lco)
  save(mi_,mi_co_,file = paste0('Outputs/CV_rel_control_default_', sample_,'.Rdata'))
  

  # check for warnings - no divergent transitions
    load(file = 'Outputs/CV_rel_control_default_5000.Rdata')
    f = c(list.files(path = here::here('Data/sim/'), pattern = 'relatedness-control_default_5000', recursive = TRUE, full.names = TRUE))
    for(i in f){
      #i=f[1]
      load(i)
      print(paste('no',i))
      print(tryCatchWEM(summary(mp_no)))
      print(paste('yes',i))
      print(tryCatchWEM(summary(mp_yes)))
    }  

# END