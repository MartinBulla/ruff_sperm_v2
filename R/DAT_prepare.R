# =============================================================
# ❗ The script runs relative to the project's root directory,
# requires DAT_morphometrics.csv, sperm_sampling.xlsx, 
# ruff_males_Seewiesen.xlsx, motility.xlsx, and DAT_HL.txt,
# is sourced via other scripts (i.e. doesn't run otherwise as
# it requires r-packages, makes composite morhpology measures,
# and prepares datasets on individual morphological or 
# velocity measurements, as well as datasets with a single
# morpho and velocity value per male, including CV
# =============================================================

# load data 
  x = fread(here::here('Data/DAT_morphometrics.csv')) 
  setnames(x,old = 'pic', new = 'sperm_ID')
  x[, sample_ID:=as.character(sample_ID)]

  # for each bird 10 sperm measured from one of the two sampling occasions
    #bb = b[part == 'Tail']
    #bbx = data.table(table(bb$bird_ID,bb$month))
    #bbx[!N%in%c(0,10), unique(V1)]

  # mneasurements from manipulated pictures vs rest are the same - so no need to control for
    #ggplot(x[part == 'Tail'], aes(x = manip, y = Pixels)) + geom_boxplot()

# make composite measures      
  d1 = x[,.(bird_ID, sample_ID, sperm_ID, part, Pixels)]
  
  d2 = d1[, sum(Pixels), by = list(bird_ID, sample_ID, sperm_ID)]
  names(d2)[4] = 'Pixels'
  d2[,part := 'Total']

  d3 = d1[part%in%c('Midpiece', 'Tail'), sum(Pixels), by = list(bird_ID, sample_ID, sperm_ID)]
  names(d3)[4] = 'Pixels'
  d3[,part := 'Flagellum']

  d4 = d1[part%in%c('Acrosome', 'Nucleus'), sum(Pixels), by = list(bird_ID, sample_ID, sperm_ID)]
  names(d4)[4] = 'Pixels'
  d4[,part := 'Head']

  b = merge(d1,d2, all = TRUE)
  b = merge(b,d3, all = TRUE)
  b = merge(b,d4, all = TRUE)
  b[, Length_µm:=Pixels*0.078]

# add metadata 
  b = merge(b, x[,.(bird_ID, sample_ID, sperm_ID, part, manip)], all.x = TRUE, by=c('bird_ID', 'sample_ID', 'sperm_ID', 'part'))
  b[manip%in%"", manip:=NA]
  # adjust multiple samples from the same time to the one recorded for motility
    b[sample_ID==51, sample_ID:=52]
    b[sample_ID==3, sample_ID:=4]
  s = data.table(read_excel(here::here('Data/sperm_sampling.xlsx'), sheet = 1))
  s = s[!is.na(recording)]
  m = data.table(read_excel(here::here('Data/ruff_males_Seewiesen.xlsx'), sheet = 1))#, range = "A1:G161"))
  m[, hatch_year:=as.numeric(substr(Hatch_date,1,4)) ]
  m[, age := 2021-hatch_year]

  s = merge(s,m[,.(Ind_ID, Morph, age)], by.x = 'DB_ID', by.y = 'Ind_ID', all.x = TRUE)
  # adjust IDs (where missed typed) for smooth merging
    b[bird_ID=='A027121649', bird_ID:='AO27121649']
    b[bird_ID=='A03999183', bird_ID:='AO3999183']
    b[bird_ID=='A0414173NL', bird_ID:='AO414-17-3NL']
    b[bird_ID=='A0414175', bird_ID:='AO414-17-5']
    b[bird_ID=='A059381830', bird_ID:='AO59381830']
    b[bird_ID=='A07422181', bird_ID:='AO7422181NL']
    b[bird_ID=='A079561654', bird_ID:='AO79561654']
    b[bird_ID=='A079561656', bird_ID:='AO79561656']
    b[bird_ID=='A079561718', bird_ID:='AO7956-17-18']
    b[bird_ID=='AIFA-016507', bird_ID:='AIFAO16507']
    b[bird_ID=='AIFA01511', bird_ID:='AIFAO15-11']
    b[bird_ID=='A079561718', bird_ID:='AO79561718']
    b[bird_ID=='G20005', bird_ID:='G200055']
    b[bird_ID=='cz005', bird_ID:='CZ005']
    b[bird_ID=='AO414-17-3NL', bird_ID:='AO414-17-3-NL']

  b = merge(b, s[,.(sample_ID, bird_ID, Morph, age, datetime, type, sperm, recording, rec_measured, month, location, avi)], by.x = c('sample_ID', 'bird_ID'), by.y = c('sample_ID', 'bird_ID'), all.x = TRUE)
  b[, loc_avi := paste(location, avi)]
  b[location == 9 & avi == 3, treat := "prison"]
  b[is.na(treat), treat := "mixed"]
  
  # fwrite(ss[!duplicated(bird_ID), .(DB_ID, bird_ID)], file = 'Data/used_males_for_Clemens.csv')
  #b = (b[!duplicated(b)]) # shouldn't be necessary,

# add motility on to individual morpho measurements (only one motility value per individual)
  d = data.table(read_excel(here::here('Data/motility.xlsx'), sheet = 1))
  d[, motileCount_ln := log(motileCount)]
  setnames(d, old=c('date','ID'), new = c('month','bird_ID'))

  # for bird 1339 motility measured from May sample, but sperm from June, so meta data changed to May sample 
     b[sample_ID==183, month:='May']
     b[sample_ID==183, datetime:=as.POSIXct('2021-05-04 12:33')]
     b[sample_ID==183, sample_ID:=95]
  
  b = merge(b, d, by = c('bird_ID', 'month'), all.x = TRUE)

  #b[is.na(Morph), unique(bird_ID)]
 
  b[Morph == 'F', Morph := 'Faeder']
  b[Morph == 'I', Morph := 'Independent']
  b[Morph == 'S', Morph := 'Satellite']
  b[is.na(issues), issues := 'zero']

  b[, part := factor(part, levels = c('Acrosome', 'Nucleus', 'Midpiece', 'Tail','Head','Flagellum', 'Total'))]
  b[part %in% c('Acrosome', 'Nucleus', 'Midpiece', 'Tail'), measure := 'original']
  b[part %in% c('Head','Flagellum', 'Total'), measure := 'composite']

  #nrow(d) # N 139

# add metadata to motility dataset
  # check for multiple samples per individuals at a single sampling period & keep only one for merging
    ajm = data.table(table(s$bird_ID,s$month))
    #s[bird_ID%in%ajm[N>1, unique(V1)]]
    ss = s[rec_measured %in% 'yes']    
  d = merge(d, ss[,.(sample_ID, bird_ID, Morph, age, datetime, type, sperm, recording, month, location, avi)], by = c('month','bird_ID'), all.x = TRUE)
   
  d[is.na(Morph), Morph := 'Zebra finch']
  d[Morph == 'F', Morph := 'Faeder']
  d[Morph == 'I', Morph := 'Independent']
  d[Morph == 'S', Morph := 'Satellite']
  d[is.na(issues), issues := 'zero']

  d[, loc_avi := paste(location, avi)]
  d[location == 9 & avi == 3, treat := "prison"]
  d[is.na(treat), treat := "mixed"]

# prepare for correlations and repeatability
  bw = reshape(b[, .(bird_ID, Morph, age, datetime, month, location, avi, loc_avi, sample_ID, sperm_ID, VAP, VSL, VCL, motileCount, prop_motile, part, Length_µm)], idvar = c("bird_ID", "Morph", "age", "datetime", "month", "location", "avi", "loc_avi","sample_ID", "sperm_ID", "VAP", "VSL", "VCL", "motileCount", "prop_motile"), timevar = "part", direction = "wide")
  setnames(bw,old = c('Length_µm.Acrosome', 'Length_µm.Nucleus','Length_µm.Head','Length_µm.Midpiece','Length_µm.Tail','Length_µm.Flagellum', 'Length_µm.Total'), new = c('Acrosome', 'Nucleus', 'Head','Midpiece', 'Tail','Flagellum','Total'))
  # add relative measures
    bw[, Midpiece_rel := Midpiece/Total]
    bw[, Flagellum_rel := Flagellum/Total]

 dw = reshape(d[bird_ID %in% d[duplicated(bird_ID), bird_ID], .(bird_ID, species, Morph, age, month, location, avi, loc_avi, VAP, VSL, VCL, motileCount, motileCount_ln,prop_motile)], idvar = c("bird_ID", "species", "Morph", "age"), timevar = "month", direction = "wide")

# mean/male dataset
  a = b[, list(mean(Length_µm), mean(VAP), mean(VSL), mean(VCL), mean(motileCount), mean(prop_motile)), by = list(month, location, avi, loc_avi, bird_ID, Morph, age, part)]
   setnames(a, old = c('V1', 'V2','V3','V4','V5','V6'), new = c('Length_avg', 'VAP', 'VSL', 'VCL', 'motileCount','prop_motile'))
  a1 = data.table(bird_ID = unique(a$bird_ID), blind = c(rep(c('Independent','Satellite', 'Faeder'), floor(length(unique(a$bird_ID))/3)), 'Independent','Satellite'))
  a = merge(a,a1, all.x = TRUE)

  aw = reshape(a, idvar = c("month", "location", "avi", "loc_avi","bird_ID", "Morph", "blind", "age", "VAP", "VSL", "VCL", "motileCount",'prop_motile'), timevar = "part", direction = "wide")
  names(aw) = c("bird_ID", "month", "location", "avi", "loc_avi", "Morph", "age", "VAP", "VSL", "VCL", "motileCount",'prop_motile', "blind", as.character(unique(a$part)))
  # add relative measures
  aw[, Midpiece_rel := Midpiece/Total]
  aw[, Flagellum_rel := Flagellum/Total]

# add HL    
    # source(here::here('R/DAT_HL.R'))
    z = fread("Data/DAT_HL.txt")
    aw = merge(aw,z[,.(sampleid, HL)], by.x = 'bird_ID', by.y = 'sampleid')
    a = merge(a,z[,.(sampleid, HL)], by.x = 'bird_ID', by.y = 'sampleid')

# morph as factor
  b[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
  bw[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
  a[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
  aw[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
  
  d[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder", "Zebra finch"))] 
  dw[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder", "Zebra finch"))] 
  #ggplot(dw, aes(x = VAP.May, y = VAP.June))+geom_point() + stat_smooth(method = 'lm') + geom_abline(intercept = 0, slope = 1, lty = 3, col = 'red')

# CV dataset
    cv_ =  b[, cv(Length_µm), by = list(bird_ID, part, Morph, loc_avi)]
    cv_[ , Morph123 := as.numeric(Morph)]
    setnames(cv_, 'V1', 'CV')
    cv_ = merge(cv_,z[,.(sampleid, HL)], by.x = 'bird_ID', by.y = 'sampleid')

# dataset for relative measurements
  bt = b[part == 'Total',.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount, motileCount_ln, prop_motile)]
  
  b1 = b[part%in%c('Midpiece'),.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount, motileCount_ln, prop_motile, location, avi, loc_avi)]
  b1[, Length_rel := Length_µm/bt$Length_µm]
 
  b2 = b[part %in% c("Flagellum"), .(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount, motileCount_ln, prop_motile, location, avi, loc_avi)]
  b2[, Length_rel := Length_µm/bt$Length_µm]
  
  br = rbind(b1,b2)
  br$Length_µm = NULL
  br[, loc_avi := paste(location, avi)]
  br[location == 9 & avi == 3, treat := "prison"]
  br[is.na(treat), treat := "mixed"]

  at = a[part == "Total", .(Morph, bird_ID, part, Length_avg, VAP, VSL, VCL, motileCount, prop_motile)]
  
  a1 = a[part %in% c("Midpiece"), .(Morph, bird_ID, part, Length_avg, VAP, VSL, VCL, motileCount, prop_motile, location, avi, loc_avi)]
  a1[, Length_rel := Length_avg/at$Length_avg]
 
  a2 = a[part %in% c("Flagellum"), .(Morph, bird_ID, part, Length_avg, VAP, VSL, VCL, motileCount, prop_motile, location, avi, loc_avi)]
  a2[, Length_rel := Length_avg/at$Length_avg]
  
  ar = rbind(a1,a2)
  ar$Length_avg = NULL
  ar = merge(ar,z[,.(sampleid, HL)], by.x = 'bird_ID', by.y = 'sampleid')

# dataset for correlations
  h = b[part == 'Head']
  setnames(h, old = "Length_µm", new="Head_µm")
  h[, Acrosome_µm := b[part == 'Acrosome',.(Length_µm)]]
  h[, Nucleus_µm := b[part == 'Nucleus',.(Length_µm)]]
  h[, Midpiece_µm := b[part == 'Midpiece',.(Length_µm)]]
  h[, Tail_µm := b[part == 'Tail',.(Length_µm)]]
  h[, Flagellum_µm := b[part == 'Flagellum',.(Length_µm)]]
  h[, Total_µm := b[part == 'Total',.(Length_µm)]]
  h[, Midpiece_rel := br[part == 'Midpiece',.(Length_rel)]]
  h[, Flagellum_rel := br[part == 'Flagellum',.(Length_rel)]]

  ha = a[part == 'Head']
  setnames(ha, old = "Length_avg", new="Head_µm")
  ha[, Acrosome_µm := a[part == 'Acrosome',.(Length_avg)]]
  ha[, Nucleus_µm := a[part == 'Nucleus',.(Length_avg)]]
  ha[, Midpiece_µm := a[part == 'Midpiece',.(Length_avg)]]
  ha[, Tail_µm := a[part == 'Tail',.(Length_avg)]]
  ha[, Flagellum_µm := a[part == 'Flagellum',.(Length_avg)]]
  ha[, Total_µm := a[part == 'Total',.(Length_avg)]]
  ha[, Midpiece_rel := ar[part == 'Midpiece',.(Length_rel)]]
  ha[, Flagellum_rel := ar[part == 'Flagellum',.(Length_rel)]]

# export for Mihai
  a$CV = cv_$CV[match(a$bird_ID, cv_$bird_ID)]
  #save(b,br,a,ar, file = 'Data/ruff_sperm_for_Mihai.RData')
