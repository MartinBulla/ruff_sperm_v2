# packages
  require(here)
  require(lme4qtl)
  require(related)
  require(reshape2)
  require(data.table)
  require(MCMCglmm)
  require(gap)
  require(matrixcalc)
  require(pedigreemm)

# function from https://stats.stackexchange.com/questions/18563/estimating-random-effects-and-applying-user-defined-correlation-covariance-struc#
  relmatmm <- function (formula, data, family = NULL, REML = TRUE, relmat = list(), 
    control = list(), start = NULL, verbose = FALSE, subset, 
    weights, na.action, offset, contrasts = NULL, model = TRUE, 
    x = TRUE, ...) {
    mc <- match.call()
    lmerc <- mc
    lmerc[[1]] <- as.name("lmer")
    lmerc$relmat <- NULL
    if (!length(relmat)) 
        return(eval.parent(lmerc))
    stopifnot(is.list(relmat), length(names(relmat)) == length(relmat))
    #lmerc$doFit <- FALSE
    lmf <- eval(lmerc, parent.frame())
    relfac <- relmat
    relnms <- names(relmat)
    stopifnot(all(relnms %in% names(lmf$FL$fl)))
    asgn <- attr(lmf$FL$fl, "assign")
    for (i in seq_along(relmat)) {
        tn <- which(match(relnms[i], names(lmf$FL$fl)) == asgn)
        if (length(tn) > 1) 
            stop("a relationship matrix must be associated with only one random effects term")
        Zt <- lmf$FL$trms[[tn]]$Zt
        relmat[[i]] <- Matrix(relmat[[i]][rownames(Zt), rownames(Zt)], 
            sparse = TRUE)
        relfac[[i]] <- chol(relmat[[i]])
        lmf$FL$trms[[tn]]$Zt <- lmf$FL$trms[[tn]]$A <- relfac[[i]] %*% Zt
    }
    ans <- do.call(if (!is.null(lmf$glmFit)) 
        lme4:::glmer_finalize
    else lme4:::lmer_finalize, lmf)
    ans <- new("pedigreemm", relfac = relfac, ans)
    ans@call <- match.call()
    ans
}
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
  load(file = 'Data/data_for_J.Rdata') # aa
  load(file = 'Data/DAT_rel-mat.RData') # m
  str(m)
  
  # various matrix setups, but none works
    mi = solve(m) # invert the matrix (Jarrod's advice)
    mip = Matrix(mi, sparse = TRUE)

    m0 = pmax(m,0) # in case it is needed, change negative estimates to 
    
    mp = Matrix(m, sparse = TRUE)  #dsCMatrix
    str(mp)

    m0i = solve(m0) # invert the matrix (Jarrod's advice)
    m0ip = Matrix(m0i, sparse = TRUE) 
      
    mp0 = pmax(mp,0) # in case it is needed, change negative estimates to 
    

# extract model residuals 
  
  mm = lm(scale(Length_avg) ~ Morph, aa)
  aa[, res := resid(mm)]

  bb = b[part == 'Nucleus'] 
  mm = lmer(scale(Length_µm) ~ Morph + (1|bird_ID), bb)
  bb[, res := resid(mm)]

# check whether residuals confounded by relatedness   
  # lme4qtl
    n_ = unique(bb$bird_ID)
    n_[!n_%in%colnames(m)]
    colnames(m)[!colnames(m)%in%n_]

    mo = relmatLmer(res ~ 1 + (1|bird_ID), data = bb, relmat = list(myID = m0i))
    mo = relmatLmer(res ~ 1 + (1|bird_ID), data = bb, relmat = list(myID = m0i))
  
  # rematmm hack - https://stats.stackexchange.com/questions/18563/estimating-random-effects-and-applying-user-defined-correlation-covariance-struc#
    m <- relmatmm(res ~ 1 + (1|bird_ID), relmat=list(gen=m0ip), data=bb)
  # MCMCglmm
    prior1 <- list(
          G = list(G1 = list(V = 1, nu = 0.002)),
          R = list(V = 1, nu = 0.002)
          )

    model1 <- MCMCglmm(res ~ 1,
        random = ~animal, ginv = list(animal = mip),
        data = aa, prior = prior1
      )

    p.var <- var(aa$res, na.rm = TRUE)
    prior2 <- list(
          G = list(G1 = list(V = matrix(p.var * 0.95), nu = 1)),
          R = list(V = matrix(p.var * 0.05), nu = 1)
        )

    model2 <- MCMCglmm(scale(res) ~ 1,
        random = ~animal, ginv = list(animal = mp),
        data = data.frame(aa), prior = prior2
      )

  # gap 
    m1 = MCMCgrm(res ~ 1, prior= prior1.1, data =aa, GRM = m)
    m1 = MCMCgrm(res ~ 1, prior= prior1.1, data =aa, GRM = mp)
    m1 = MCMCgrm(res ~ 1, prior= prior1.1, data =aa, GRM = m0ip)