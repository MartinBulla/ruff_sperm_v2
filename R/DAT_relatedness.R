require(related)
require(reshape2)

readgenotypedata()


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

# TESTNING
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