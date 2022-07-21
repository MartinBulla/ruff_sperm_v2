# packages
  require(arm)
  require(here)
  require(lme4qtl)
  require(related)
  require(reshape2)
  require(data.table)
  require(gap)
  require(matrixcalc)

# no need to run - makes kinship matrix
    input <- readgenotypedata('Data/Dat_GenotypesRuff_sampled.txt')  
    rel <- coancestry(input$gdata, lynchli=1)
    r = data.table(rel$relatedness)
    #r[ind1.id == 'AA_n001aripo' & ind2.id == 'AA_n001aripo']   

    rr = r[,.(pair.no,ind1.id,ind2.id,lynchli)]

    rx = data.table(pair.no = (nrow(r)+1):(nrow(r)+length(unique(c(r$ind2.id,r$ind1.id)))), ind1.id = unique(c(r$ind2.id,r$ind1.id)),ind2.id = unique(c(r$ind2.id,r$ind1.id)), lynchli = 1)

    ra = rbind(rr,rx)
    ra = ra[order(ind1.id, ind2.id)]

    m = acast(ra,ind1.id~ind2.id, value.var = 'lynchli')   

    for(i in 1:length(unique(ra$pair.no))) {  
      #i = 2 
      ToDo <- ra[pair.no==unique(ra$pair.no)[i]]
      if(ToDo$ind1.id == ToDo$ind2.id){ next } else { m[which(rownames(m)==ToDo$ind2.id),which(colnames(m)==ToDo$ind1.id)] <- ToDo$lynchli
          }
          
      }    
    save(m, file = 'DATA/DAT_rel-mat.RData')

# load data
  load(file = 'Data/data_for_A.Rdata') # bb
  load(file = 'Data/DAT_rel-mat.RData') # m
  #str(m)
  
  # various matrix setups
    mi = solve(m) # inverted the matrix (Jarrod's advice)
    m0 = pmax(m,0) # in case it is needed, change negative estimates to zero
    m0i = solve(m0) # invert the matrix (Jarrod's advice)    

# extract model residuals 
  
  bb = b[part == 'Midpiece'] 
  mm = lmer(scale(Length_µm) ~ Morph + (1|bird_ID), bb)
  bb[, res := resid(mm)]
  bb[, animal := bird_ID]

# check whether residuals confounded by relatedness   

  mo = relmatLmer(res ~ 1 + (1|bird_ID), data = bb, relmat = list(myID = m0i))
  summary(mo)
  # check whether all individuals within bird_ID present within the matrix
    n_ = unique(bb$bird_ID)
    n_[!n_%in%colnames(m)]
    colnames(m)[!colnames(m)%in%n_]    

# brms
  cores_ = 2
  chains_ = 2
  iter_ = 5000
  thin_ = 10
  
  mp_no = brm(res ~ 0 + Intercept + (1|bird_ID), data = bb, cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))
    
  
  Amat = m0
  diag(Amat) = diag(Amat)+0.2 # making matrix positive def https://discourse.mc-stan.org/t/positive-definite-vs-positive-semidefinite-required-for-phylogenetic-regression-error-message-in-brms/17049

  mp_m= brm(res ~ 0 + Intercept  + (1|bird_ID) + (1 | gr(animal, cov = Amat)), data = bb,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))
    
  # mp_yes +0,1
  # mp_yes2 +0.00001  

  #mp_m
  #mp_m0i
  #mp_m

  Amat = m0i
  diag(Amat) = diag(Amat)+0.1 
  mp_m= brm(res ~ 0 + Intercept  + (1|bird_ID) + (1 | gr(animal, cov = Amat)), data = bb,  data2 = list(Amat = Amat), cores = cores_, chains = chains_, iter = iter_, thin = thin_, seed = 5,  control = list(adapt_delta = 0.99), sample_prior="yes",save_pars = save_pars(all = TRUE))
  

  mdat <- matrix(c(1,0,1, 0,1,0, 1,0,1), nrow = 3, ncol = 3, byrow = TRUE)
  mdat = jitter(mdat,0.1)
  diag(mdat) =1
  mdat[upper.tri(mdat)] <- mdat[upper.tri(mdat)]

  M1<-matrix(1:25,ncol=5)
  diag(M1) =1
  M1[upper.tri(M1)]<-t(M1)[upper.tri(M1)]
  M1
# END