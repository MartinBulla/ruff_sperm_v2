require(related)
require(reshape2)

# load data
  # merge datasets
    n1 = data.table(read_excel(here::here("Data/Ruffs_microsat_2020.xlsx"), sheet = "allBins", na = "NA"))
    n1 = n1[Sex == 1]
    n1 = n1[!is.na(OriginalRing)]
    n1 = n1[!OriginalRing %in% '7 -04 - 105']
    n1 = n1[,c('sort','Population',"newLAB#","eggID","Sex","Generation","Type", "Location","DNAconc", "DateGenotyping", "notes") := NULL]
    n2 = data.table(read_excel(here::here("Data/Ruffs_microsat_2010-19.xlsx"), sheet = "allAlleles", na = "NA", col_type = 'text'))
    n2 = n2[sex == 1]
    n2 = n2[!is.na(OriginalRing)]
    n2 = n2[,c('year',"eggID","sex","Lab_ID","DNA_notes", "notes") := NULL]
    setcolorder(n2, names(n1))
    n2 = n2[!OriginalRing%in%n1$OriginalRing]
    n = rbind(n1,n2)
    n = n[,c('Mother',"Father") := NULL]
    #n[, c('Ruff8a', 'Ruff8b') := NULL]

    change_columns = names(n)[2:length(names(n))]
    n[ ,(change_columns) := lapply(.SD, as.numeric),
           .SDcols = change_columns]
    nn = n

  # adjust bird IDs
    nn[ OriginalRing== '7 -04 - 105', OriginalRing := 704105]
    nn[ OriginalRing== 'AIF AO - 15 - 11', OriginalRing := 'AIFAO15-11']
    nn[ OriginalRing== 'A 8209 AIF AO 16 507', OriginalRing := 'AIFAO16507']
    nn[ OriginalRing== 'AO 2712-16-49-B 6.5', OriginalRing := 'AO27121649']
    nn[ OriginalRing== 'AO 3999-18-3-NL 6.5', OriginalRing := 'AO3999183']
    nn[ OriginalRing== 'AO 414-17-3-NL 6.5', OriginalRing := 'AO414-17-3-NL']
    nn[ OriginalRing== 'AO 414-17-5-NL 6.5', OriginalRing := 'AO414-17-5']
    nn[ OriginalRing== 'AO 5938-18-30-NL 6.5', OriginalRing := 'AO59381830']
    nn[ OriginalRing== 'AO 7422-18-1-NL 6.5', OriginalRing := 'AO7422181NL']
    nn[ OriginalRing== 'AO 7422-18-1-NL 6.5', OriginalRing := 'AO7422181NL']
    nn[ OriginalRing== 'AO 7956-17-18-NL 6.5', OriginalRing := 'AO7956-17-18']
    nn[ OriginalRing== 'AO 7956-16-54-NL 6.5', OriginalRing := 'AO79561654']
    nn[ OriginalRing== 'AO 7956-16-56-NL 6.5', OriginalRing := 'AO79561656']
    nn[ OriginalRing== 'CZ 6.5 005', OriginalRing := 'CZ005']
    nn[ OriginalRing== 'G 5.0 - 14 - 223', OriginalRing := 'G14223']
    nn[ OriginalRing== 'G 5.0 - 17 - 267', OriginalRing := 'G17267']
    nn[substr(OriginalRing, 1,3) == 'G20', OriginalRing := paste0('G20', substring(OriginalRing, 5,8))]
    nn[ OriginalRing== 'G 5.0 - 14 -221', OriginalRing := 'G1221'] 
    nn[ OriginalRing== 'AO 3999-18-4-NL 6.5', OriginalRing := 'AO3999-18-4-NL'] 

  fwrite(nn,file = 'Data/Dat_GenotypesRuff.txt', col.names = FALSE,sep="\t")
  input <- readgenotypedata('Data/Dat_GenotypesRuff.txt')  

# test estimators
  compareestimators(input, 100) # we compared Correlation Coefficients Between Observed & Expected Values for each estimator. The lynchli estimator correlates best with the expected values (although by small margin) and hence we use this estimator. 

  Correlation Coefficients Between Observed & Expected Values:
        wang        0.885598
        lynchli     0.889864
        lynchrd     0.788119
        quellergt   0.880375

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
#fwrite(file="Data/DAT_rel-mat.txt", m)
save(m, file = 'DATA/DAT_rel-mat.RData')
#write.table(m, file="Data/DAT_rel-mat_test.txt", row.names=TRUE, col.names=TRUE)

#read(file = 'DATA/DAT_rel-mat.RData')

# TESTNING
 input <- readgenotypedata('Data/GenotypeData.txt')
        compareestimators(input, 100)

        #---Estimate relatedness---#
        rel <- coancestry(input$gdata, lynchli=1, lynchrd=1, quellergt=1, wang=1)
        str(rel$relatedness)
        #---Create simulated individuals of known relatedness---#
        sim <- familysim(input$freqs, 100)

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
fwrite(file="Data/DAT_rel-mat.txt", m)
fwrite(file="Data/DAT_rel-mat.txt", m)
save(m, file = 'DATA/DAT_rel-mat.RData')
write.table(m, file="Data/DAT_rel-mat_test.txt", row.names=TRUE, col.names=TRUE)


mm = fread(file="Data/DAT_rel-mat.txt")
mm = as.matrix(read.table("Data/DAT_rel-mat.txt"))

x = data.table(relatedness = c(1,0,0,1,0.5, 1), ind1 = c('a','a','a','b','b','c'), ind2 = c('a','b','c','b','c','c')) 


x = data.table(relatedness = c(0,0,0.5), ind1 = c('a','a','b'), ind2 = c('b','c','c'))
x2 = data.table(relatedness = 1, ind1 = unique(c(x$ind1, x$ind2)),ind2 = unique(c(x$ind1, x$ind2)))
x = rbind(x, x2)        
x[,pair:=paste(ind1,ind2)]
m = acast(x,ind1~ind2, value.var = 'relatedness')   

for(i in 1:length(unique(x$pair))) {  
    #i = 2 
    ToDo <- x[pair==unique(x$pair)[i]]
    if(ToDo$ind1 == ToDo$ind2){ next } else { m[which(rownames(m)==ToDo$ind2),which(colnames(m)==ToDo$ind1)] <- ToDo$relatedness
        }
    }


   relatedness ind1 ind2 
a1 = acast(x,ind1~ind2, value.var = 'relatedness')
a2 = acast(x,ind2~ind1, value.var = 'relatedness')

dat1 <- read.table(textConnection(a1), header = TRUE)
dat2 <- read.table(textConnection(a2), header = TRUE)


dat <- data.frame(relatedness = c(0,0,0.5), ind1 = c('a','a','b'), ind2 = c('b','c','c'))   
dat$pair <- paste(dat$ind1,dat$ind2, sep="__")

inds <- sort(unique(c(dat$ind1,dat$ind2)))
rel.mat <- data.frame(matrix(rep(1,length(inds)*length(inds)),ncol=length(inds)))
colnames(rel.mat) <- inds ; rownames(rel.mat) <- inds

for(i in 1:length(unique(dat$pair))) {
    ToDo <- subset(dat, pair==unique(dat$pair)[i])
    ind1 <- strsplit(ToDo$pair[1],"__",fixed=TRUE)[[1]][1]
    ind2 <- strsplit(ToDo$pair[1],"__",fixed=TRUE)[[1]][2]
    if(ind1==ind2) { next } else { rel.mat[which(rownames(rel.mat)==ind1),which(colnames(rel.mat)==ind2)] <- ToDo$relatedness }
    if(ind1==ind2) { next } else { rel.mat[which(rownames(rel.mat)==ind2),which(colnames(rel.mat)==ind1)] <- ToDo$relatedness }
 }




matrix(data = x$m, dimnames = list(c('a','b','c'), c('a','b','c')))
matrix(dimnames = list(c('a','b','c'), c('a','b','c')))

 matrix(, nrow = c('a','b','c'), ncol = c('a','b','c'))
 x = matrix(, nrow = 3, ncol = 3)
 rownames(x) <- c('a','b','c')
 colnames(x) <- c('a','b','c')