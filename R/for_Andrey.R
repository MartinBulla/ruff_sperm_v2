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

# check whether residuals confounded by relatedness   

  mo = relmatLmer(res ~ 1 + (1|bird_ID), data = bb, relmat = list(myID = m0i))
  summary(mo)
  # check whether all individuals within bird_ID present within the matrix
    n_ = unique(bb$bird_ID)
    n_[!n_%in%colnames(m)]
    colnames(m)[!colnames(m)%in%n_]    

# END