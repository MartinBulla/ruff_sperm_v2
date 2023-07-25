# =======================================================================================
# ❗ The script runs relative to the project's root directory,
# generates homozygosity by locus using 
#   the microsat data (see DAT_genotypes.R) and
#   GENHETv3.1 function (Coulon 2010, https://doi.org/10.1111/j.1755-0998.2009.02731.x) 
# and exports the dataset to ./DATA/
# =======================================================================================

# tools
    require(data.table)
    require(gtools)
    require(readxl)

    source('R/GENHETv3.1.R')

# load data
    n = fread("Data/Dat_GenotypesRuff.txt")
    b = fread(here::here("Data/DAT_morphometrics.csv"))[, unique(bird_ID)] # bird IDs in the sperm dataset

# identify all loci 
    lo = names(n)[2:length(names(n))]
    lo = substring(lo,1,nchar(lo)-1)
    lo = lo[!duplicated(lo)]

# run GENHET function from "GENHETv3.1.R"
    
    Htest = GENHET(dat = data.frame(n), estimfreq = "T", locname = lo) # all data
    df=data.table(Htest)

    Htest2=GENHET(dat=data.frame(n[OriginalRing%in%b]),estimfreq="T",locname=lo) # sperm-sampled individuals only
    df2=data.table(Htest2)

    df1 = df[,.(sampleid,HL)]
    setnames(df1, old = 'HL', new = 'HL_all')
    df3 = merge(df2,df1, all.x = TRUE)

    df3[, HL := as.numeric(HL)]
    df3[, HL_all := as.numeric(HL_all)]
    
    require(ggplot2); ggplot(df3, aes(x = HL_all, y = HL)) + geom_jitter(alpha =0.5) + geom_abline(slope =1) # estimates from full and sperm sampled only data are similar

# export the estimates
    fwrite(file="Data/DAT_HL.txt", df3)

# END