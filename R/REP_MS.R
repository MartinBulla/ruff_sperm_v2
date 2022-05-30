#' ---
#' title: "Lack of major differences in sperm traits of highly divergent male reproductive morphs in the ruff"
#' author: "Martin Bulla"
#' date: "`r Sys.time()`"
#' output: 
#'     html_document:
#'         toc: true
#'         toc_float: true
#'         toc_depth: 5
#'         code_folding: hide
#'         bibliography: ruff_sperm.xml
#'         link-citations: yes
#' ---

#+ r setup, include=FALSE 
knitr::opts_chunk$set(message = FALSE, warning = FALSE, cache = TRUE)



#' # Code to load tools and data
  # TOOLS 
    require(here)
    source(here::here('R/tools.R'))
    
    require(arm) 
    require(effects)
    require(ggpubr)
    require(ggsci) 
    require(grid)
    require(gridExtra)
    require(magrittr)
    require(multcomp)
    require(PerformanceAnalytics)
    require(rptR) 
    require(stringi)
    require(viridis)
    require(readxl)


    # constants
      round_ = 3 # number of decimal places to round model coefficients
      nsim = 5000 # number of simulations to extract estimates and 95%CrI
      ax_lines = "grey60" # defines color of the axis lines
      colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)
    
    # functions
      getime = function (x) {ifelse(is.na(x), as.numeric(NA), as.numeric(difftime(x, trunc(x,"day"), units = "hours")))}
      
      getDay = function (x) {as.Date(trunc(x, "day"))}
     
      # for R output based on sim
        R_out = function(name = "define", model = m, nsim = 5000){
         bsim <- sim(model, n.sim=nsim)  
         l=data.frame(summary(model)$varcor)
         l = l[is.na(l$var2),]
         l$var1 = ifelse(is.na(l$var1),"",l$var1)
         l$pred = paste(l$grp,l$var1)

         q50={}
         q025={}
         q975={}
         pred={}
         
         # variance of random effects
         for (ran in names(bsim@ranef)) {
           #ran =names(bsim@ranef)[1]
           ran_type = l$var1[l$grp == ran]
           for(i in ran_type){
              # i = ran_type[2]
            q50=c(q50,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.5)))
            q025=c(q025,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.025)))
            q975=c(q975,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.975)))
            pred= c(pred,paste(ran, i))
            }
           }
         # residual variance
         q50=c(q50,quantile(bsim@sigma^2, prob=c(0.5)))
         q025=c(q025,quantile(bsim@sigma^2, prob=c(0.025)))
         q975=c(q975,quantile(bsim@sigma^2, prob=c(0.975)))
         pred= c(pred,'Residual')

         ci = c(round(100*q025/sum(q025))[1], round(100*q975/sum(q975))[1])
         ci = ci[order(ci)]
         
         ri=data.table(model = name, repeatability=paste0(round(100*q50/sum(q50)),'%')[1], CI = paste0(paste(ci[1], ci[2], sep ="-"), '%'))
         
         
         return(ri)
         }
  # DATA 
    # composite measures
      x = fread(here::here('Data/morphometrics.csv')) 
      setnames(x,old = 'pic', new = 'sperm_ID')
      x[, sample_ID:=as.character(sample_ID)]
      
      # for each bird 10 sperm measured from one of the two sampling occasions
        #bb = b[part == 'Tail']
        #bbx = data.table(table(bb$bird_ID,bb$month))
        #bbx[!N%in%c(0,10), unique(V1)]

        #  mneasurements from minip vs rest are the same - so no need to control for
        #ggplot(x[part == 'Tail'], aes(x = manip, y = Pixels)) + geom_boxplot()
      
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

      b = merge(b, s[,.(sample_ID, bird_ID, Morph, age, datetime, type, sperm, recording, rec_measured, month)], by.x = c('sample_ID', 'bird_ID'), by.y = c('sample_ID', 'bird_ID'), all.x = TRUE)
      
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
      d = merge(d, ss[,.(sample_ID, bird_ID, Morph, age, datetime, type, sperm, recording, month)], by = c('month','bird_ID'), all.x = TRUE)
     
      d[is.na(Morph), Morph := 'Zebra finch']
      d[Morph == 'F', Morph := 'Faeder']
      d[Morph == 'I', Morph := 'Independent']
      d[Morph == 'S', Morph := 'Satellite']
      d[is.na(issues), issues := 'zero']

    # prepare for correlations and repeatability
      bw = reshape(b[,.(bird_ID,Morph, age, datetime, month, sample_ID, sperm_ID, VAP,VSL,VCL, motileCount, part, Length_µm)], idvar = c('bird_ID','Morph','age','datetime', 'month', 'sample_ID', 'sperm_ID','VAP','VSL','VCL', 'motileCount'), timevar = 'part', direction = "wide")  
      setnames(bw,old = c('Length_µm.Acrosome', 'Length_µm.Nucleus','Length_µm.Head','Length_µm.Midpiece','Length_µm.Tail','Length_µm.Flagellum', 'Length_µm.Total'), new = c('Acrosome', 'Nucleus', 'Head','Midpiece', 'Tail','Flagellum','Total'))
      # add relative measures
        bw[, Midpiece_rel := Midpiece/Total]
        bw[, Flagellum_rel := Flagellum/Total]

     dw = reshape(d[bird_ID%in%d[duplicated(bird_ID), bird_ID],.(bird_ID,species, Morph, age, month, VAP,VSL,VCL, motileCount, motileCount_ln)], idvar = c('bird_ID','species','Morph','age'), timevar = 'month', direction = "wide")  
   
    # mean/male dataset
      a = b[, list(mean(Length_µm), mean(VAP),mean(VSL), mean(VCL), mean(motileCount)), by = list(month, bird_ID, Morph, age, part)]
       setnames(a, old = c('V1', 'V2','V3','V4','V5'), new = c('Length_avg', 'VAP', 'VSL', 'VCL', 'motileCount'))
      a1 = data.table(bird_ID = unique(a$bird_ID), blind = c(rep(c('Independent','Satellite', 'Faeder'), floor(length(unique(a$bird_ID))/3)), 'Independent','Satellite'))
      a = merge(a,a1, all.x = TRUE)

      aw = reshape(a, idvar = c('month','bird_ID','Morph', 'blind','age','VAP','VSL','VCL', 'motileCount'), timevar = "part", direction = "wide")
      names(aw) = c('bird_ID','month','Morph', 'age','VAP','VSL','VCL', 'motileCount', 'blind', as.character(unique(a$part)))
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
        cv_ =  b[, cv(Length_µm), by = list(bird_ID, part, Morph)]
        cv_[ , Morph123 := as.numeric(Morph)]
        names(cv_) [4]='CV'
        cv_ = merge(cv_,z[,.(sampleid, HL)], by.x = 'bird_ID', by.y = 'sampleid')

    # dataset for relative measurements
      bt = b[part == 'Total',.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount, motileCount_ln)]
      
      b1 = b[part%in%c('Midpiece'),.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount, motileCount_ln)]
      b1[, Length_rel := Length_µm/bt$Length_µm]
     
      b2 = b[part%in%c('Flagellum'),.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount, motileCount_ln)]
      b2[, Length_rel := Length_µm/bt$Length_µm]
      
      br = rbind(b1,b2)
      br$Length_µm = NULL

      at = a[part == 'Total',.(Morph, bird_ID,  part, Length_avg, VAP, VSL, VCL, motileCount)]
      
      a1 = a[part%in%c('Midpiece'),.(Morph, bird_ID, part, Length_avg, VAP, VSL, VCL, motileCount)]
      a1[, Length_rel := Length_avg/at$Length_avg]
     
      a2 = a[part%in%c('Flagellum'),.(Morph, bird_ID, part, Length_avg, VAP, VSL, VCL, motileCount)]
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

#' ## To decide
#' - who shall be a co-author apart from us two, Kim, Tomas & Michael? Surely guess Dove and Wolfgang, perhaps Clemens, Jasmin, Katrin?

***
# '## Abstract
#' Ruff *Calidris pugnax* is a Palearctic lekking shorebird with three strikingly different mating morphs: aggressive ‘independents’, semi-cooperative ‘satellites’ and female-mimic ‘faeders’[@Hogan-Warburg1966; @vanRhijn1991; @Widemo1998; @Jukema2006]. An autosomal inversion [@Lank1995; @Kupper2015; @Lamichhaney2015] drives the developments of the major morph-differences in mating behaviours, ornamentation, body and relative testis size. However, whether the morphs also differ in sperm traits has not been investigated. Here, we used a captive-breeding population of ruffs to reveal that in comparison to independents, satellites have longer tails and hence sperm (because tails determine the total sperm length) and perhaps shorter nucleus and hence heads,and faeders had longer midpiece. In addition, there was a tendency for shorter nucleaus (and hence heads) in satellites  investigate morph-differences in  sperm morphology and velocity. We found

 
# Background 
Aggressive 'independents' constitute 80-95% of male ruffs, show spectacular diversity of predominantly dark ornamental plumage, and are dominant holders of display sites on leks [@Hogan-Warburg1966; @vanRhijn1991; @Widemo1998]. Submissive and slightly smaller 'satellites' constitute 5-20% of males, show predominantly white ornamental plumage and do not defend display sites [@Hogan-Warburg1966; @Höglund1989; @vanRhijn1991; @Widemo1998]. 'Independents' behave differently with 'satellites' than they do with other 'independents'. Presence of satellites on the lek assists independents with female attraction and allows satellites to 'steal' copulations [@Widemo1998]. Rare 'faeders' (<1%) mimic females in their plumage and smaller size and attempt rapid copulations (i.e. 'steal' copulations) when females solicit matings from ornamented displaying males [@Jukema2006]. However, how 'faeders' behave in the the wild is largely unknown. On the leg, the dominant 'independents' get multiple copulations and multiple females [@Lank2002; @Vervoort2019]. The likelihood of copulation for other ‘independents’ and ‘satellites’ is driven by the amount of time an individual spends on the lek [@Vervoort2019], whereas 'faeders' have the least copulation opportunity. 

An autosomal inversion is associated with the major differences in body and relative testis size, ornamentation, and mating behaviours of 'satellites' and 'faeders' [@Lank1995; @Kupper2015; @Lamichhaney2015]. In other words, the inversion is linked to how the three morphs invest in the pre-copulatory male-male competition. However, whether the morphs also differ in post-copulatory male-male competition, such as ejaculate volume and sperm density and/or sperm morphology and velocity, is unknown.

We investigated differences in sperm morphology and velocity among the three male ruff morphs. We sampled the sperm from a captive breeding population, maintained since 1985. Adult males were phenotyped for mating behavior ('independent', 'satellite' or 'faeder') and the presence of ornamental plumage (ornamented versus faeder). Some adult females were also phenotyped, with small size indicating carriers of Faeder [@Lank2013] and aggressive behavior following testosterone implantation permitting classification as independents or satellites[@Lank1999]. I IMAGINE THAT RECENTLY WE USE GENETICAL MORPH DETERMINATION, RIGHT?


# Research questions 
1. Do the three ruff morphs differ in the sperm morphology and velocity, as well as in the within- and between-individual variability of the sperm traits?  - see figure 1 [here](https://www.nature.com/articles/s41559-017-0236-1?platform=hootsuite)
2. Does sperm morphology predict sperm swimming speed? - see figure 2 [here](https://www.nature.com/articles/s41559-017-0236-1?platform=hootsuite)

# Predictions

## 1. Do the three ruff morphs differ in the sperm morphology and velocity?  

Mating system of the ruffs appears to be well studied [@Hogan-Warburg1966; @Höglund1989; @vanRhijn1991; @Lank1995; @Kupper2015; @Lamichhaney2015; @Vervoort2019], but rare 'faeders' (female mimics) are tough to observe in the wild, and hence we know little about their behaviour in the wild. For predictions about sperm morphology and velocity based on morph-specific behaviour, we need to know whether the frequency of copulation on the lek is linked to the siring success. Whereas we know that dominant independents get most of the copulations on the lek [@Vervoort2019], we do not know whether females that copulate with faeders also always copulate with other males. In other words, we need to know  whether faeder sperm may experience different level of competition than sperm of other two morphs. Consequently, instead of clear predictions, we provide XX scenarios that we envision for sperm traits among the three morphs.

First, we do not expect between-morph differences in total sperm length because in birds sperm length strongly coevolved with the length of female sperm storage tubules [@Briskie1993] and ruffs show no hint of assortative mating by morph. Since males of all three morphs attempt to mate with the same pool of females (**unpublished data**), they should have sperm of similar length. Moreover, all three morphs are under selection to optimize sperm, and a steeper selection gradient in 'faeders' due to its low chance of copulating (**cit**) does not mean that a beneficial mutation enhancing sperm speed will occur there rather than elsewhere. This predictions of “no difference” are mirrored in a recent review (Kustra and Alonzo 2020). Across taxa there are no consistent differences in sperm traits between dominant and sneaker males. Moreover, despite the general expectation that sneaker males should produce sperm that are more competitive (e.g. higher quality or performance), the existing theory does not predict explicitly how sneaker males should differ in sperm traits.


Second, the three morphs may differ in sperm traits, if genes regulating their development reside on the autosomal inversion that is behind the phenotypes of 'satelites' and 'faeders'. Indeed, the inversion contains a gene GAS8 that has been associated with sperm velocity, at least in mice [@Yeh2002]. The gene GAS8 lies in a region with the strongest differences between 'independents' and the other two morphs. Hence, we expect 'independents' to have a slower sperm than 'satellites' and 'faeders'.

Third,the morphs may also differ in sperm traits if investing in pre-copulatory male-male competition trade-offs with investing in post-copulatory male-male competition. We see gradient in pre-copulatory investment with 'independents' investing the most on the lek and 'faeders' the least, i.e. 'faeders' having the lowest frequency and chance of copulations [@Hogan-Warburg1966; @Höglund1989; @vanRhijn1991; @Jukema2006]. This gradient in pre-copulatory investment is perhaps mirrored in the post-copulatory investment because independents have the smallest and faeders the largest relative (but not absolute) testes size [@Jukema2006; @Kupper2015; and Figure 1 based on Loveland et al in prep]. Given these differences in frequency and chance of copulation, 'faeders' are expected to deliver the largest ejaculates with faster sperm and sperm optimized for speed. Specifically, we predict that 'faeders' have the longest sperm, midpiece and tail or midpiece and  flagellum length relative to the whole sperm length and have sperm that is fastest of the three morphs (Figure 2 & 3). As 'faeders' invest little in pre-copulatory male-male competintion, we also expect 'faeders' sperm to differ the most from the sperm of 'independents', while we expect 'satellites', who partly invest in pre-copulatory male-male competition, but to a lesser extant than 'independents' [@Vervoort2019], to have intermediate sperm between the two (Figure 2 & 3). 

```{r Testes,  warning = FALSE, message = FALSE, fig.align="center", fig.width=3, fig.height=4}
  d = data.table(read_excel(here::here('Data/GSI_for_Liam.xlsx'), sheet = 1))#, range = "A1:G161"))
  d[, Morph := factor(Morph, levels=c("Res", "Sat", "Faed"))] 
  d[Morph == 'Res', Morph := 'Independent']
  d[Morph == 'Sat', Morph := 'Satellite']
  d[Morph == 'Faed', Morph := 'Faeder']

# plot
  g1 = ggplot(d, aes(x = Morph, y = Gonadmass)) + 
    geom_boxplot() + 
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'darkgrey', aes(fill =Morph))+
    scale_fill_viridis(discrete=TRUE)+
    ylab('Gonad Mass (g)') +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")

  g2 = ggplot(d, aes(x = Morph, y = GSI)) + 
    geom_boxplot() +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'darkgrey', aes(fill =Morph))+
    scale_fill_viridis(discrete=TRUE)+
    ylab('Gonadosomatic Index') +
    theme_bw()+theme(legend.position = "none")  

  grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))
 
```
**Figure 1** | Morph-specific testes size.

```{r Predictions,  warning = FALSE, message = FALSE, fig.align="center", fig.width=3, fig.height=3}
  set.seed(1)
  d = data.table(morph = factor(c( rep('independent', 60),
                                   rep('satelite', 30),
                                   rep('faeder', 15)
                                   ), 
                                levels = c('independent','satelite','faeder')
                                ),
                sperm_trait =  c(rnorm(60, mean = 1, sd = 0.5), 
                                 rnorm(30, mean = 2, sd = 0.5),
                                 rnorm(15, mean = 3, sd = 0.5) 
                                 ) 
                )
  #dev.new(width = 3.5, height = 3.5)
  ggplot(d, aes(x = morph, y = sperm_trait, col = morph)) +  
    geom_jitter() + theme_MB +  theme(legend.position="none") + xlab('Morph') + ylab('Sperm length or relative midpiece length') +
    scale_colour_manual(values = colors)
 
```
**Figure 2** | Predicted relationship between morph type and sperm traits.

Fourth, the morphs may also differ in sperm traits if the morphs differ in the intensity of sperm competition. Between species, sperm length, coefficient of variation in sperm length, and swimming speed as well as investment in sperm numbers (testes size) is generally related to the intensity of sperm competition [@Kleven2008; @Lifjeld2010] ADD CIT. Whether intensity of sperm competition may differ among the three morphs of a single species is unclear. We could imagine a scenario where females that copulate with 'faeders' always copulate with other males (including males of other morphs) whereas females that copulate with 'independent' tend to copulate only with that male (or morph). Whereas all individuals (regardless of morph) will benefit from having competitive sperm, 'faeders' sperm might be under stronger selection. Also, given the rarity of 'faeders' and of their copulations, the sperm of 'independents' and 'satellites' may rarely compete with sperm of 'faeders'. Thus, if 'faeder' sperm would evolve a competitive advantage, the other two morphs may not "follow" upon such advantage. If so, we predict that 'faeders' will have the fastest sperm and sperm morphology optimized for speed (the longest sperm with longest midpiece, tail and hence flagellum). Consequently, we also expect the smallest coefficient of variation in the sperm morphometry as within- and between-individual variability in sperm traits is lower in species with higher sperm competition [@Kleven2008; @Lifjeld2010]. This expectation is based on the assumption that the sperm production is better optimized and less error prone in 'faeders' than in the other two morphs. Such optimizations shall increase 'faeders' chances of inseminating a female despite 'faeders' only sporadic chances to copulate. 

## 2. Does sperm morphology predict sperm swimming speed?

Sperm morphology, specifically long sperm, sperm with short heads, and long midpiece (relative to overall sperm length), has been positively associated with sperm swimming speed [**CIT**]. We thus predict that longer sperm and/or sperm with longer head, and/or sperm with relatively longer midpiece will swim faster.

```{r Predictions2,  warning = FALSE, message = FALSE, fig.align="center", fig.width=4, fig.height=4}
  set.seed(1)
  d = data.table(morph = c( rep('independent', 60),
                                   rep('satelite', 30),
                                   rep('faeder', 15)
                                   ),
                midpiece =  c(rnorm(60, mean = 2, sd = 0.5), 
                                 rnorm(30, mean = 3, sd = 0.5),
                                 rnorm(15, mean = 4, sd = 0.5) 
                                 ),
                
                tail =  c(rnorm(60, mean = 2, sd = 0.5), 
                                 rnorm(30, mean = 3, sd = 0.5),
                                 rnorm(15, mean = 4, sd = 0.5) 
                                 ),
                velocity =  c(rnorm(60, mean = 2, sd = 0.5), 
                                 rnorm(30, mean = 3, sd = 0.5),
                                 rnorm(15, mean = 4, sd = 0.5) 
                                 )
                )
  d[,morph123 :=ifelse(morph == 'independent', 1, ifelse(morph == 'satelite', 2,3))]             

  colors_ <- colors[d$morph123]
  par(las = 1, cex.axis = 0.6, cex.lab = 0.8, cex.main = 0.8)
  s3d=scatterplot3d(d$midpiece, d$tail, d$velocity, pch = 16, type="h", 
              color=colors_, grid=TRUE, box=FALSE,
              xlab = "",
              ylab = "",
              zlab = "",
              x.ticklabs=c("short","","","","","long"),
              y.ticklabs=c("short","","","","long",""),
              z.ticklabs=c("slow","","","","","fast"),
              mar = c(3, 2, 0, 1.5)
              )     
  text(x = 7.5, y = 1, "Tail", srt = 0, cex = 0.8)
  text(x = 2.5, y = -0.5, "Mipiece", srt = 0,xpd = TRUE, cex = 0.8)
  text(x = -0.5, y = 2.5, "Velocity", srt = 90,xpd = TRUE, cex = 0.8)
  legend("bottom", legend = levels(factor(d$morph,levels = c('independent','satelite','faeder'))),
      col =  colors, pch = 16,xpd = TRUE, horiz = TRUE,inset = -0.125, bty = "n", cex = 0.7)
 
```
**Figure 3** | Predicted relationship between morph type and sperm traits.


# Methods #
## Individuals 
In 2018, we collected and stored in ~5% formalin the vas deferens from 15 sacrificed males (5 per morph) from the **Simon Fraser** University **colony**. The colony was founded with 110 ruffs hatched from wild eggs collected in Finland in 1985, 1989 and 1990 plus two 'faeder' males from the Netherlands in 2006 [@Lank1995; @Lank2013]. 

In 2021, we plan to collect sperm of all males from the **Max Planck colony** (`r nrow(s[morph == 'I'])` 'independents', `r nrow(s[morph == 'S'])` 'satelites', `r nrow(s[morph == 'F'])+4` 'faeders'; of which `r nrow(sv)` males were sampled already in Vancouver: `r nrow(sv[morph == 'r'])`  'independents', `r nrow(sv[morph == 's'])` 'satelites' and `r nrow(sv[morph == 'f'])` 'faeders'). The colony was founded in 2018 from the Simon Fraser colony. 

## Housing - Needed?

## Sperm sampling
Compared to passerines, shorebirds do not have a cloaca protuberance. Thus, abdominal massage and cloacal lavage method of sperm collection (e.g. @Knief2017 ; see the [specific protocol](https://raw.githack.com/MartinBulla/ruff_sperm/main/Protocols/protocol_sperm.html)) often lead to unclean samples contaminated by excrements. Such samples are then unsuitable for measuring sperm velocity. Hence, in 2021 along with abdominal massage and cloacal lavage, we plan to collect sperm by electro‐stimulation [@Lierz2013]. The length and diameter of the electro-stimulation probe, as well as the electric current and the number of electric impulses will be adapted to the sampled individuals. The electro-stimulation should result in clean sperm samples suitable for velocity measurements. 

The morphs will be sampled in a haphazard order to avoid temporal trends in morph sampling. In other words, all males of one morph will not be sampled at the same time. Sperm will be pipetted from the cloaca (~0.5–3μl) and immediately diluted in a preheated (40°C) Dulbecco’s Modified Eagle’s Medium (Advanced D-MEM, Invitrogen, USA). For the velocity measurements, an aliquot will be pipetted onto a standard 20μm two-chamber count slide (Leja, The Netherlands) placed on a thermal plate (Tokai Hit, Tokai Hit Co., LtD., Japan) kept at 40 °C. For the morphology measurements, an aliquot will be pipetted onto a microscopy slide and the rest of the sperm sample will be fixed in 100μl ~5% formalin solution.

## Sperm morphometry
From each sperm sample, we will pipet ~10μl onto a microscopy slide, let it dry at room temperature, and - in case of formalin-fixed samples - wash the slide gently under  distilled water to wash away the dried formalin and phosphate-buffered saline solution. Sperm from the vas deferens will be first extracted by cutting a piece from the middle of the tubules, plucking apart the tissue with curved extraction forceps to release the sperm, and washing away the sperm with the phosphate phosphate-buffered saline solution. 

We will inspect the slides with a light microscope ('Zeiss Axio Imager.M2') under 200x magnification.  To aid distinction of sperm parts, dried slides will be stained with Hoechst 33342 (which identifies the sperm nucleus) and Mitotracker Green FM (which identifies sperm midpiece) not earlier than 48h before photographing. For step by step slide preparation protocol see [here](https://raw.githack.com/MartinBulla/ruff_sperm/main/Protocols/protocol_sperm.html).

```{r Pic1, echo=FALSE, fig.cap="**Picture 1** | Ruff sperm blue-stained nucleus by Hoechst 33342 and green-stained midpiece by Mitotracker Green FM", echo=FALSE}
knitr::include_graphics(here::here("Pictures/Pic1_staining_example.jpg"))
```
For each sample, we will photograph at least 10 intact normal-looking spermatozoa, with a 12 megapixel (4250 × 2838) digital camera ('Zeiss Axiocam 512 color' with pixel size of 3.1μm × 3.1μm). Ten sperm per male were previously used to investigate the relationship between coefficient of variation in sperm trait and sperm competition [@Kleven2008; @Lifjeld2010]. For each sperm, we will measure the length of the acrosome, the nucleus, the midpiece and the tail to the nearest 0.1μm using the open source software [ImageJ](https://imagej.nih.gov/ij/) (for a detailed sperm measuring protocol see [here](https://www.dropbox.com/s/ouqfqih4gzx6rq5/Protocol_imageJ_sperm_measuring_manual.docx?dl=0)). To minimize observer error, all measurements will be taken by one person (KT), who will be blind to the individual morph and will measure the sperm from all morphs and samples in a random order (NOTE THAT THIS IS NOT THE CASE FOR THE VAS DEFERENS).

We will calculate total sperm length as the sum of all parts, head length as the sum of acrosome and nucleus length, and flagellum length as the sum of midpiece and tail length. We will compute coefficients of variation for each trait (CV = [SD/mean]) both within males and between males of each morph. If we are unable to measure 10 sperm for all males, we will correct the coefficients of variation for variation in sample size (CVadj = [1 + 1/(4n)]*CV; @Sokal1981). 

Each measured sperm and sperm part will be numbered and referenced in the database (see [sperm measuring protocol](https://www.dropbox.com/s/ouqfqih4gzx6rq5/Protocol_imageJ_sperm_measuring_manual.docx?dl=0) and [database](https://www.dropbox.com/s/9uk7bcxmd4o1prx/ruff_sperm_morphology.xlsx?dl=0)). This ensures transparency and allows re-measurement of the same sperm by the same or different person.

In addition, we will scan the slides for abnormal sperm by counting 100 sperms and noting how many are abnormal and in what manner

## Sperm velocity
For each sperm sample we will record sperm velocity for approximately 45s in eight different fields of the slide under a 100x magnification using a phase contrast and a digital camera (UI-1540-C, Olympus) mounted on a microscope (CX41, Olympus) fitted with a thermal plate (Tokai Hit, Tokai Hit Co., LtD., Japan) kept at a constant temperature of 40°C. The videos will be recorded XX SECONDS AFTER SAMPLE COLLECTION at 25 frames per second. A single person will analyze each recorded field using the CEROS computer-assisted sperm analysis system (Hamilton Thorne Inc., Beverly, Massa- chusetts, USA), visually inspect the tracked objects and exclude non-sperm objects and static spermatozoa from the analysis [@Laskemoen2010; @Cramer2016; @Opatová2016]. The dilution medium does not contain any spermatozoa attractants to guide the spermatozoa towards one direction; thus, we will use  curvilinear velocity rather than straight-line velocity as our measurement of sperm swimming speed [@Laskemoen2010]. We will report how many sperm cells per sperm sample will be recorded (median, 95%CI, range) and - if necessary - control for the (log-transformed) number of measured sperm cells per sample. 

We will test the method beginning of May and adjust it if necessary. Passerine sperm move straight, by turning like a screw around their axis (with their screw shaped heads). It is uncertain how sperm of ruffs move. If we find the above method problematic, we adjust the method, e.g.  use different (deeper) size of Leja slides or sperm attractants. 

## Inbreeding
Because our sperm samples came and will come from males bred in captivity, we expect higher levels of inbreeding compared with males from wild populations. In birds and mammals (but not in insects) inbred males have a higher proportion of abnormal sperm and lower sperm velocity than outbred males [@Gomendio2000; @Heber2013; @Ala-Honkola2013; @Opatová2016]. Thus, our sperm velocity measurements may not reflect sperm velocity of wild ruffs. However, there is no evidence that the  morphology  of normal‐looking sperm (e.g., length, coefficient of variation) differs between inbred and outbred males [@Mehlis2012; @Ala-Honkola2013; @Opatová2016]. Based on these studies, we assume that our measurements reflect the variation in sperm morphology observed in wild ruffs. Nevertheless, we know the relatedness of individuals within our population (pedigree) and can estimate inbreeding from microsatellite markers. Hence, we will investigate whether variation in inbreeding and genetic diversity is linked to the sperm traits and if so we will attempt to control for it in our models. Similarly, if find between morph differences in sperm traits, we will investigate whether such differences might have occurred by chance due to population bottlenecks.

## Analysis plan 
Sample sizes will reflect the maximum available data. No data selection will be done conditional on the outcome of statistical tests. To avoid unconscious data dragging, we will first analyze the data on sperm traits blind to the males morphs. Specifically, we will first randomly assign male morph to the measurements, explore such data and finalize the analyze using such data. Only then we will run the analyses with true assignment of the morphs. We will report all results, all data exclusions, all manipulations and all measures in the study at the [GitHub repository](https://github.com/MartinBulla/ruff_sperm). 

## Statistical analyses
All analyses will be in R [@R-Core-Team2020] using lmer() and glmer() functions of the lme4 package to fit linear and generalized linear mixed-effects models, and - if necessary - one of the package for fitting the pedigree structure as a random effect (e.g. [brms](https://mikheyev.github.io/brms-wam/), MCMCglmm, pedigreeMM). We will estimate the variance explained by the fixed effects of our mixed-effects models as marginal R2-values, using the r.squaredGLMM() function of the MuMIn package (v1.15.6). Model fits will be visually validated (by qq-plots of residuals and plots of residuals against fitted values). 

In general, we will fit two sets of linear models. A first set of models will have as a response variable the raw sperm trait measurements and will be controlled for multiple measurements per male by fitting male identity as a random effect. A second set of models will have average sperm trait per male as a response and - if necessary - will be controlled for the number of sampled sperm per individual (N) by specifying a weights term as sqrt(N-3)(*cit needed*). Note that we aim to keep the N constant. For each morph, we will exclude individuals with <10 sperm, as long as these represent <5% of individuals sampled within each morph. MB THE FIRST SET OF MODLES IS RATHER COMPLEX (if pedigree is used), THAT IS WHY I PROPOSE ALSO THE SIMPLER AVERAGE VARIANT. THE AVERAGE VARIANT WILL ALSO CORRESPOND WITH MODELS FOR VELOCITY and CV. HOWEVER, IF THIS MAKES THINGS TO COMPLICATED, WE CAN RUN JUST THE FIRST SET.

To investigate whether sperm traits differ between the morphs, we will fit a model for each sperm trait separately, with the male morph (three-level factor) as an explanatory variable, controlling for the aviary, in which an individual is housed (aviary ID), and - if necessary - controlling for the pedigree structure.

To investigate whether sperm morphology explains variation in sperm velocity, we will first fit a linear model with average velocity per male as the dependent variable, and with male-average total sperm length, head (*acrosome + nucleus*), midpiece and tail length (**?or rather flagelum - midpiece+tail?**), as well as their squared terms (after mean-centering) and all two-way interactions between the three linear terms, as explanatory variables. If necessary, we will control for a pedigree structure. In case, some of the sperm traits substantially correlate (Pearson's r>0.6), we will use only one of the correlated traits. We will then investigate whether predicted sperm velocity corresponds with observed sperm velocities. 

We will investigate whether the sperm traits are representative of each male and morph by calculating the repeatability of sperm measurements per male and morph, which will be obtained through 1,000 parametric bootstrap iterations [@Stoffel2017].

<span style="color:grey">**THE FOLLOWING IS LIKELY NOT NEEDED:** To evaluate multicollinearity between all main effect predictors, we will estimate their variance inflation factor using the corvif() function [@Zuur2009] in R. A general guideline is that a variance inflation factor larger than 5 or 10 is large, indicating that the model has problems estimating the coefficient. However, this in general does not degrade the quality of predictions. If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response. Thus, if VIF is larger than 5 we can use the model for predictions.</span>

# Proposed time line
- ongoing - measuring sperm morphology from vas deferens samples
- May - data collection
- June - July - photographing and measuring sperm samples, quantifying sperm velocity
- August - September - analyses and 1st draft

# TO DECIDE
- Does it matter whether we let the microscope slides dry at room temperature or whether we heat fix them?
- What sperm traits to use? Total length with/without acrosome, acrosome, head (acrosome + nucleus), midpiece, tail, flagelum (midpiece + tail). BART: total length, head, midpiece and tail. Wolfi: also flagellum
- Shall we use discriminant function analyses on the two most different sperm traits to predict the sperm morph? BART: Perhaps use all sperm traits, calculate PCA and plot PC1 vs PC2 with morph in different colours? Not sure this is necessary, but might be a nice way to illustrate the differences (if any)? WOLFGANG: As a note, DFA will always give stronger separation than PCA. PCA is likely to be unimpressive, and DFA is a bit HARKing, as it also separates groups (to some extent) when there is no true difference (e.g. randomly generated data).
- Do we need to control for the pedigree, if residuals of the models are randomly distributed, i.e. there is no pedigree structure in residuals? WOLFGANG:This is a matter of heritability. I would say that the pedigree is important to control for analyses of sperm morphology, but not velocity or CV.
- Which packages to use for modeling: [brms](https://mikheyev.github.io/brms-wam/), MCMCglmm, pedigreeMM
- When recording sperm velocity, shall we also make an additional record in a cellulose solution - a viscous environment, which "mimics" conditions within the female reproductive tract [@Schmoll2020]? BART: No, if it compromises number of sampled individuals. a viscous environment; TOM: If we proceed fast with basic measurements. If cellulose treatment slows us significantly, I would skip it at least for the first year

# TO ADD
- videos to accompany the protocols (vas deferens preparation, slide preparation, sperm collection and video recording)

# DECIDED
- We use 10sperm/male for the CV; after seeing  within male variation in sperm sizes and using Wolfgang's extrapolation we will decide whether more sperm per individual are needed and worth it. *We can randomly subset into 2x5 sperm per male, calculate 2x CV from 5, then get a repeatability of those CVs from 5, and we can then extrapolate the repeatability of our CV from 10. Surely this will be quite low.*
- Clemens will provide us with four of their faeders.

# References