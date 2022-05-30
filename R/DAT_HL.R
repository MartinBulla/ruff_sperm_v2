# generates homozygosity by locus from the microsat data and using GENHETv3.1 function (Coulon 2010;  https://doi.org/10.1111/j.1755-0998.2009.02731.x) 

# tools
    require(data.table)
    require(gtools)
    require(readxl)

    source('R/GENHETv3.1.R')

# merge datasets
    n1 = data.table(read_excel(here::here("Data/Ruffs_microsat_2020.xlsx"), sheet = "allBins", na = "NA"))
    n1 = n1[Sex == 1]
    n1 = n1[!is.na(OriginalRing)]
    n1 = n1[!OriginalRing %in% '7 -04 - 105']
    n1 = n1[,c('sort','Population',"newLAB#","eggID","Sex","Generation","Type", "Location","Mother","Father","DNAconc", "DateGenotyping", "notes") := NULL]
    n2 = data.table(read_excel(here::here("Data/Ruffs_microsat_2010-19.xlsx"), sheet = "allAlleles", na = "NA", col_type = 'text'))
    n2 = n2[sex == 1]
    n2 = n2[!is.na(OriginalRing)]
    n2 = n2[,c('year',"eggID","sex","Lab_ID","Mother","Father","DNA_notes", "notes") := NULL]
    setcolorder(n2, names(n1))
    n2 = n2[!OriginalRing%in%n1$OriginalRing]
    n = rbind(n1,n2)
    n[, c('Ruff8a', 'Ruff8b') := NULL]

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

# identify all loci 
    lo = names(n)[2:length(names(n))]
    lo = substring(lo,1,nchar(lo)-1)
    lo = lo[!duplicated(lo)]

# run GENHET function from "GENHETv3.1.R" (
    Htest=GENHET(dat=data.frame(nn),estimfreq="T",locname=lo)
    df=data.table(Htest)

    Htest2=GENHET(dat=data.frame(nn[OriginalRing%in%aw$bird_ID]),estimfreq="T",locname=lo)
    df2=data.table(Htest2)

    df1 = df[,.(sampleid,HL)]
    setnames(df1, old = 'HL', new = 'HL_all')
    df3 = merge(df2,df1, all.x = TRUE)

    df3[, HL := as.numeric(HL)]
    df3[, HL_all := as.numeric(HL_all)]
    ggplot(df3, aes(x = HL_all, y = HL)) + geom_point() + geom_abline(slope =1 )

# export the estimates
    fwrite(file="Data/DAT_HL.txt", df3)

# end    