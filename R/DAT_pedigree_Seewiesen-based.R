require(data.table)
require(gtools)
require(readxl)

# load data
  # merge datasets & fix issues
    n1 = data.table(read_excel(here::here("Data/Ruffs_microsat_2020.xlsx"), sheet = "allBins", na = "NA"))
    #n1 = n1[Sex == 1]
    n1 = n1[!is.na(OriginalRing)]
    n1 = n1[,c('Population',"OriginalRing","Sex",'Mother', 'Father', "notes")]
    n2 = data.table(read_excel(here::here("Data/Ruffs_microsat_2010-19.xlsx"), sheet = "allAlleles", na = "NA", col_type = 'text'))
    n2 = n2[!is.na(OriginalRing)]
    n2 = n2[,c('year',"OriginalRing","sex",'Mother', 'Father', "notes")]
    setnames(n2, old = 'sex', new = 'Sex')
    
    n = merge(n1,n2, all=TRUE) #
    nn=n[!is.na(Mother) | !is.na(Father)]
    #nn[OriginalRing==261] 
    nn[OriginalRing==1005 & Father == 301, notes:= 'father previously 288'] 
    nn = nn[!which(OriginalRing==1005 & Father == 288)] 
    nn[OriginalRing==1310 &  Father == 264, notes:= 'father previously 333'] 
    nn = nn[!which(OriginalRing==1310 &  Father == 333)] 
    nn[OriginalRing == 1332, Father :=302]

    nn = nn[Mother != "561.80600000000004"]
    nn = nn[OriginalRing!=1374]

 # adjust bird IDs
    nn[ OriginalRing== '7 -04 - 105', OriginalRing := 704105]
    nn[ Father== '7 -04 - 105', Father := 704105]
    nn[ OriginalRing== 'AIF AO - 15 - 11', OriginalRing := 'AIFAO15-11']
    nn[ OriginalRing== 'A 8209 AIF AO 16 507', OriginalRing := 'AIFAO16507']
    nn[ Father== 'A 8209 AIF AO 16 507', Father := 'AIFAO16507']
    nn[ OriginalRing== 'AO 2712-16-49-B 6.5', OriginalRing := 'AO27121649']
    nn[ Father== 'AO 2712-16-49-B 6.5', Father := 'AO27121649']
    nn[ OriginalRing== 'AO 3999-18-3-NL 6.5', OriginalRing := 'AO3999183']
    nn[ Father== 'AO 3999-18-3-NL 6.5', Father := 'AO3999183']
    nn[ OriginalRing== 'AO 414-17-3-NL 6.5', OriginalRing := 'AO414-17-3-NL']
    nn[ OriginalRing== 'AO 414-17-5-NL 6.5', OriginalRing := 'AO414-17-5']
    nn[ OriginalRing== 'AO 5938-18-30-NL 6.5', OriginalRing := 'AO59381830']
    nn[ OriginalRing== 'AO 7422-18-1-NL 6.5', OriginalRing := 'AO7422181NL']
    nn[ Father== 'AO 7422-18-1-NL 6.5', Father := 'AO7422181NL']
    nn[ OriginalRing== 'AO 7956-17-18-NL 6.5', OriginalRing := 'AO7956-17-18']
    nn[ OriginalRing== 'AO 7956-16-54-NL 6.5', OriginalRing := 'AO79561654']
    nn[ OriginalRing== 'AO 7956-16-56-NL 6.5', OriginalRing := 'AO79561656']
    nn[ OriginalRing== 'CZ 6.5 005', OriginalRing := 'CZ005']
    nn[ Father== 'CZ 6.5 005', Father := 'CZ005']
    nn[ OriginalRing== 'G 5.0 - 14 - 223', OriginalRing := 'G14223']
    nn[ Father== 'G 5.0 - 14 - 223', Father := 'G14223']
    nn[ OriginalRing== 'G 5.0 - 17 - 267', OriginalRing := 'G17267']
    nn[substr(OriginalRing, 1,3) == 'G20', OriginalRing := paste0('G20', substring(OriginalRing, 5,8))]
    nn[ OriginalRing== 'G 5.0 - 14 -221', OriginalRing := 'G1221'] 
    nn[ OriginalRing== 'AO 3999-18-4-NL 6.5', OriginalRing := 'AO3999-18-4-NL'] 

    setnames(nn, old='OriginalRing', new = 'bird_ID')  

  # DONE check whether all in based on a different dataset - all in and ok
      p = fread('Data/pedigree2014to2019.txt')
      setnames(p, old='ind', new = 'bird_ID') 
      p[,bird_ID:=as.character(bird_ID)]
  
      p[!bird_ID%in%nn$bird_ID]

      np = merge(nn,p, by = 'bird_ID', all=TRUE)
      np[Mother.x!=Mother.y]
      np[Father.x!=Father.y]

   # adjust as needed
      nrow(nn[!Mother%in%nn$bird__ID])
      nrow(nn[!Father%in%nn$bird__ID])

      x = nn$bird_ID
      f = nn$Father
      m = nn$Mother
      #length(f[!f%in%x])
      #length(m[!m%in%x])

      ff = data.table(bird_ID = unique(f[!f%in%x]), Sex = 1, Mother = NA, Father = NA, notes = NA, Population = NA, year = NA) 
      mm = data.table(bird_ID = unique(m[!m%in%x]), Sex = 0, Mother = NA, Father = NA, notes = NA, Population = NA, year = NA) 
      yy =  data.table(bird_ID = c("1681","AIFAO15-11","AO414-17-3-NL", "AO414-17-5","AO59381830","AO7956-17-18", "AO79561654","AO79561656","G17267"), Sex = 0, Mother = NA, Father = NA, notes = NA, Population = NA, year = NA) 
      fmn=rbind(ff,mm,yy,nn)

      #x = fmn$bird_ID
      #y = unique(b$bird_ID)
      #y[!y%in%x]

  fwrite(fmn[,.(bird_ID, Mother, Father)],file = 'Data/Dat_parentage.txt')