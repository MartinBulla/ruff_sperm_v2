# ==========================================================================
# ❗ The script runs relative to the project's root directory,
# uses Dat_GenotypesRuff_sampled.txt to create and export relatedness matrix
# into ./DATA/
# ==========================================================================

# TOOLS
  require(data.table)
  require(here)
  require(related)
  require(reshape2)

# load data
  input <- readgenotypedata('Data/Dat_GenotypesRuff_sampled.txt')  

# DONE test estimators
  #compareestimators(input, 100) # we compared Correlation Coefficients Between Observed & Expected Values for each estimator. The lynchli estimator correlates best with the expected values (although by small margin) and hence we use this estimator. 

  #for all:
  #Correlation Coefficients Between Observed & Expected Values:
   #     wang        0.885598
   #     lynchli     0.889864
   #     lynchrd     0.788119
   #     quellergt   0.880375
  
  #for our sample
  #Correlation Coefficients Between Observed & Expected Values:
  #  wang        0.878313
  #  lynchli     0.893336
  #  lynchrd     0.831443
  #  quellergt   0.883814

# make kingship matrix
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

save(m, file = 'DATA/DAT_rel-mat.RData') # based on saampled

# END