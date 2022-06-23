require(car)
require(data.table)
require(ggplot)
require(ggpubr)
require(grid)
require(here)
require(lmodel2)
require(multcomp)
require(readxl)
require(smatr)

b = data.table(read_excel(here::here('Data/Ruff brains organs collected 2019 for Suzana.xlsx'), sheet = 1))#, range = "A1:G161"))
b[, morph := factor(morph, levels=c("R", "S", "C"))] 
b[morph == 'R', morph := 'Independent']
b[morph == 'S', morph := 'Satellite']
b[morph == 'C', morph := 'Faeder']
b[, soma := bodyweight - brainweight]

# allometry
  
  leveneTest(log10(brainweight)~morph,b) # ok

  ai=aov(log10(brainweight)~morph+log10(soma) +morph:log10(soma),b)
  Anova(ai, type="III")
  
  as=aov(log10(brainweight)~morph+log10(soma),b)
  Anova(as, type="III")


  mi = lm(log10(brainweight)~morph+log10(soma) + morph:log10(soma), b)
  summary(glht(mi))
  summary((mi))

  ms = lm(log10(brainweight)~morph+log10(soma), b)

  AIC(mi,ms)

# plot all
  summary(lm(log10(brainweight)~log10(soma), data = b))


  line_log = data.table(brainweight= exp(mean(log(b$brainweight))+c(-0.654,0,0.5)), soma = exp(mean(log(b$soma))+c(-0.65,0,0.5)), morph = 'Faeder')
  line_log10 = data.table(brainweight= 10^(mean(log10(b$brainweight))+c(-0.2,0,0.2)), soma = 10^(mean(log10(b$soma))+c(-0.2,0,0.2)), morph = 'Faeder')


  ggplot(b, aes(y = brainweight, x = soma, col = morph)) + geom_point() + 
      stat_smooth(method = 'lm') +
      stat_cor(aes(label = ..r.label..),  size = 3)+ #label.x = 1, 
      geom_line(data = line_log, aes(y = brainweight, x = soma), col = 'black', lty = 3)+
      geom_point(y = mean(log(b$brainweight)), x=mean(log(b$soma)), col = 'black')+
      scale_x_continuous(trans = 'log', 'Body - brain mass [g]')+
      scale_y_continuous(trans = 'log','Brain mass [g]')

  ggsave(file = 'Outputs/brain_body_allometry_ln.png', width = 6, height = 5)


  ggplot(b, aes(y = brainweight, x = soma, col = morph)) + geom_point() + 
      stat_smooth(method = 'lm') +
      stat_cor(aes(label = ..r.label..),  size = 3)+ #label.x = 1, 
      geom_line(data = line_log10, aes(y = brainweight, x = soma), col = 'black', lty = 3)+
      geom_point(y = mean(log10(b$brainweight)), x=mean(log10(b$soma)), col = 'black')+
      scale_x_continuous(trans = 'log10', 'Body - brain mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')

  ggsave(file = 'Outputs/brain_body_allometry_log10.png', width = 6, height = 5)

# plot F
  b = data.table(read_excel(here::here('Data/Ruff brains organs collected 2019 for Suzana.xlsx'), sheet = 1))#, range = "A1:G161"))
  b[, morph := factor(morph, levels=c("R", "S", "C"))] 
  b[morph == 'R', morph := 'Independent']
  b[morph == 'S', morph := 'Satellite']
  b[morph == 'C', morph := 'Faeder']
  b[, soma := bodyweight - brainweight]
  b=b[sex == 'F']
  summary(lm(log10(brainweight)~log10(soma), data = b))
  line_log10 = data.table(brainweight= 10^(mean(log10(b$brainweight))+c(-0.2,0,0.2)), soma = 10^(mean(log10(b$soma))+c(-0.2,0,0.2)), morph = 'Faeder')
  ggplot(b, aes(y = brainweight, x = soma, col = morph)) + geom_point() + 
      stat_smooth(method = 'lm') +
      stat_cor(aes(label = ..r.label..),  size = 3)+ #label.x = 1, 
      geom_line(data = line_log10, aes(y = brainweight, x = soma), col = 'black', lty = 3)+
      geom_point(y = mean(log10(b$brainweight)), x=mean(log10(b$soma)), col = 'black')+
      scale_x_continuous(trans = 'log10', 'Body - brain mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')+
      labs(tag = '♀')

  ggsave(file = 'Outputs/brain_body_allometry_log10_F.png', width = 6, height = 5)

# plot M
  b = data.table(read_excel(here::here('Data/Ruff brains organs collected 2019 for Suzana.xlsx'), sheet = 1))#, range = "A1:G161"))
  b[, morph := factor(morph, levels=c("R", "S", "C"))] 
  b[morph == 'R', morph := 'Independent']
  b[morph == 'S', morph := 'Satellite']
  b[morph == 'C', morph := 'Faeder']
  b[, soma := bodyweight - brainweight]
  b=b[sex == 'M']
  summary(lm(log10(brainweight)~log10(soma), data = b))
  line_log10 = data.table(brainweight= 10^(mean(log10(b$brainweight))+c(-0.2,0,0.2)), soma = 10^(mean(log10(b$soma))+c(-0.2,0,0.2)), morph = 'Faeder')
  # OSL
    ggplot(b, aes(y = brainweight, x = soma, col = morph)) + geom_point() + 
        stat_smooth(method = 'lm') +
        stat_cor(aes(label = ..r.label..),  size = 3)+ #label.x = 1, 
        geom_line(data = line_log10, aes(y = brainweight, x = soma), col = 'black', lty = 3)+
        geom_point(y = mean(log10(b$brainweight)), x=mean(log10(b$soma)), col = 'black')+
        scale_x_continuous(trans = 'log10', 'Body - brain mass [g]')+
        scale_y_continuous(trans = 'log10','Brain mass [g]')+
        labs(tag = '♂')
    ggsave(file = 'Outputs/brain_body_allometry_log10_M.png', width = 6, height = 5)

  # reduced axis regression
     # slope of allometry   
      summary(lm(log10(brainweight)~log10(soma), data = b[morph == 'Independent']))
      summary(lm(log10(brainweight)~log10(soma), data = b[morph == 'Satellite']))
      summary(lm(log10(brainweight)~log10(soma), data = b[morph == 'Faeder']))
      

    fit2 <- with(b[morph == 'Independent'], line.cis(log10(brainweight), log10(soma)))    
    fit2_i = data.frame(morph = 'Independent', intercept = fit2[1,1],slope=fit2[2,1], int_lwr = fit2[1,2], sl_lwr=fit2[2,2], int_upr = fit2[1,3], sl_upr = fit2[2,3])
 
    fit2 <- with(b[morph == 'Satellite'], line.cis(log10(brainweight), log10(soma)))  
    fit2_s = data.frame(morph = 'Satellite', intercept = fit2[1,1],slope=fit2[2,1], int_lwr = fit2[1,2], sl_lwr=fit2[2,2], int_upr = fit2[1,3], sl_upr = fit2[2,3]) 

    fit2 <- with(b[morph == 'Faeder'], line.cis(log10(brainweight), log10(soma)))   
    fit2_f = data.frame(morph = 'Faeder', intercept = fit2[1,1],slope=fit2[2,1], int_lwr = fit2[1,2], sl_lwr=fit2[2,2], int_upr = fit2[1,3], sl_upr = fit2[2,3])

    fit2 = rbind(fit2_i,fit2_s,fit2_f)

    ggplot(b, aes(y = brainweight, x = soma, col = morph)) + geom_point() + 
      geom_abline(data = fit2, aes(intercept=intercept, slope=slope, col=morph))+ 
      geom_abline(data = fit2, aes(intercept=int_lwr, slope=sl_lwr, col=morph), lty = 3)+ 
      geom_abline(data = fit2, aes(intercept=int_upr, slope=sl_upr, col=morph),lty = 3)+ 
      geom_line(data = line_log10, aes(y = brainweight, x = soma), col = 'black', lty = 3)+
      scale_x_continuous(trans = 'log10', 'Body - brain mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')+
      labs(tag = '♂')

   
    ggsave(file = 'Outputs/brain_body_allometry_log10_M_red
      uced.png', width = 6, height = 5)

# simulation 1
  # generate species data
    #d = data.table(body_mass=c(seq(5,200, by =2), seq(220, 500, by = 20), seq(550, 1000, by = 50), seq(2000, 8000, by = 1000)))
    d = data.table(body_mass=c(5,50,500,5000))
    d[ , body_mass_log10 := log10(body_mass)]
    d[ , brain_mass_log10 := -2 +0.57*body_mass_log10]
    d[ , brain_mass := 10^brain_mass_log10]
    d[ , species := 1:nrow(d)]
    d = d[rep(1:nrow(d),50),]
    d = d[order(species)]
    d[, pk :=1:nrow(d)]

    ggplot(d, aes(x = body_mass,  y =brain_mass))+
      geom_abline(intercept=-2, slope=0.57, col = 'red', lty = 3) +
      geom_point() +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')

  # add within species variation along the allometry   
    d[ , var_sp_body_log10 := rnorm(n = 1 , mean = 0, sd =abs(body_mass_log10*0.1)), by = pk]
    d[ , body_mass_log10_with_spVar := body_mass_log10+var_sp_body_log10]
    d[ , body_mass_with_spVar := 10^body_mass_log10_with_spVar]
    d[ , brain_mass_log10_spVar_based := -2 +0.57*body_mass_log10_with_spVar]
    d[ , brain_mass_spVar_based := 10^brain_mass_log10_spVar_based]

    ggplot(d, aes(x = body_mass_with_spVar,  y =brain_mass_spVar_based))+
      geom_abline(intercept=-2, slope=0.57, col = 'red', lty = 3) +
      geom_point() +
      #geom_smooth(method = 'lm', se = FALSE) +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')  

  # add random noise
    d[ , var_sp_brain_log10 := rnorm(n = 1 , mean = 0, sd =abs(brain_mass_log10*0.1)), by = pk]
    d[ , brain_mass_log10_spVar_based_noise := brain_mass_log10_spVar_based + var_sp_brain_log10]
    d[ , brain_mass_spVar_based_noise := 10^brain_mass_log10_spVar_based_noise]

    g1=
    ggplot(d, aes(x = body_mass_with_spVar,  y =brain_mass_spVar_based_noise, col = factor(species), fill = factor(species)))+
      geom_abline(intercept=-2, slope=0.57, col = 'red', lty = 3) +
      geom_point(pch = 21, alpha = 0.6) +
      geom_smooth(method = 'lm', se = FALSE, col = 'black') +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')+
      theme(legend.position = 'none')  

# simulation 2
  # generate species data
    #d = data.table(body_mass=c(seq(5,200, by =2), seq(220, 500, by = 20), seq(550, 1000, by = 50), seq(2000, 8000, by = 1000)))
    d = data.table(body_mass=c(5,50,500,5000))
    d[ , body_mass_log10 := log10(body_mass)]
    d[ , brain_mass_log10 := -2 +0.57*body_mass_log10]
    d[ , brain_mass := 10^brain_mass_log10]
    d[ , species := 1:nrow(d)]
    d = d[rep(1:nrow(d),50),]
    d = d[order(species)]
    d[, pk :=1:nrow(d)]

    ggplot(d, aes(x = body_mass,  y =brain_mass))+
      geom_abline(intercept=-2, slope=0.57, col = 'red', lty = 3) +
      geom_point() +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')

  # add within species variation along the allometry   
    d[ , body_mass_with_sp:= rnorm(n = 1 , mean = body_mass, sd =body_mass*0.1), by = pk]
    d[ , body_mass_with_sp_log10 := log10(body_mass_with_sp)]
    d[ , brain_mass_with_sp_log10 := -2 +0.57*body_mass_with_sp_log10]
    d[ , brain_mass_with_sp := 10^brain_mass_with_sp_log10]

    ggplot(d, aes(x = body_mass_with_sp,  y =brain_mass_with_sp, col = factor(species), fill = factor(species)))+
      geom_abline(intercept=-2, slope=0.57, col = 'red', lty = 3) +
      geom_point(pch = 21, alpha = 0.6) +
      #geom_smooth(method = 'lm', se = FALSE) +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')+
      theme(legend.position = 'none')    

  # add random noise
    d[ , brain_mass_with_sp_var := rnorm(n = 1 , mean = brain_mass_with_sp, sd =abs(brain_mass_with_sp*0.1)), by = pk]
    d[ , brain_mass_with_sp_var_log10 := log10(brain_mass_with_sp_var)]

    g2=
    ggplot(d, aes(x = body_mass_with_sp,  y =brain_mass_with_sp_var, col = factor(species), fill = factor(species)))+
      geom_abline(intercept=-2, slope=0.57, col = 'red', lty = 3) +
      geom_point(pch = 21, alpha = 0.6) +
      geom_smooth(method = 'lm', se = FALSE, col = 'black') +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')+
      theme(legend.position = 'none')  
# simulation 3
  # generate species data
    #d = data.table(body_mass=c(seq(5,200, by =2), seq(220, 500, by = 20), seq(550, 1000, by = 50), seq(2000, 8000, by = 1000)))
    d = data.table(body_mass=c(5,50,500,5000))
    d[ , body_mass_log10 := log10(body_mass)]
    d[ , brain_mass_log10 := -2 +0.57*body_mass_log10]
    d[ , brain_mass := 10^brain_mass_log10]
    d[ , species := 1:nrow(d)]
    d = d[rep(1:nrow(d),50),]
    d = d[order(species)]
    d[, pk :=1:nrow(d)]

    ggplot(d, aes(x = body_mass,  y =brain_mass))+
      geom_abline(intercept=-2, slope=0.57, col = 'red', lty = 3) +
      geom_point() +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')

  # add within species variation along the allometry   
    d[ , body_mass_with_sp:= rnorm(n = 1 , mean = body_mass, sd =abs(body_mass*0.2)), by = pk]
    d[ , body_mass_with_sp_log10 := log10(body_mass_with_sp)]
    d[ , brain_mass_with_sp_log10 := -2 +0.57*body_mass_with_sp_log10]
    d[ , brain_mass_with_sp := 10^brain_mass_with_sp_log10]

    ggplot(d, aes(x = body_mass_with_sp,  y =brain_mass_with_sp, col = factor(species), fill = factor(species)))+
      geom_abline(intercept=-2, slope=0.57, col = 'red', lty = 3) +
      geom_point(pch = 21, alpha = 0.6) +
      #geom_smooth(method = 'lm', se = FALSE) +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')+
      theme(legend.position = 'none')    

  # add random noise
    d[ , brain_mass_with_sp_var := rnorm(n = 1 , mean = brain_mass_with_sp, sd =abs(brain_mass_with_sp*0.2)), by = pk]
    d[ , brain_mass_with_sp_var_log10 := log10(brain_mass_with_sp_var)]

    g3=
    ggplot(d, aes(x = body_mass_with_sp,  y =brain_mass_with_sp_var, col = factor(species), fill = factor(species)))+
      geom_abline(intercept=-2, slope=0.57, col = 'red', lty = 3) +
      geom_point(pch = 21, alpha = 0.6) +
      geom_smooth(method = 'lm', se = FALSE, col = 'black') +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')+
      theme(legend.position = 'none')  
# simulation 2 also x noise - USE
  # generate species data
    #d = data.table(body_mass=c(seq(5,200, by =2), seq(220, 500, by = 20), seq(550, 1000, by = 50), seq(2000, 8000, by = 1000)))
    d = data.table(body_mass=c(5,50,500,5000))
    int_ = -1.081081 # get a precise intercept from Tsuboi et al 2018
    d[ , body_mass_log10 := log10(body_mass)]
    d[ , brain_mass_log10 := int_ +0.57*body_mass_log10] # !!! NO CLUE WHY -1.14 and not -2.3
    d[ , brain_mass := 10^brain_mass_log10]
    d[ , species := 1:nrow(d)]
    d = d[rep(1:nrow(d),1000),]
    d = d[order(species)]
    d[, pk :=1:nrow(d)]

    ggplot(d, aes(x = body_mass,  y =brain_mass))+
      geom_abline(intercept=int_, slope=0.57, col = 'red', lty = 3) +
      geom_point() +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10', 'Brain mass [g]')

  # add within species variation along the allometry   
    d[ , body_mass_with_sp:= rnorm(n = 1 , mean = body_mass, sd =body_mass*0.05), by = pk]
    d[ , body_mass_with_sp_log10 := log10(body_mass_with_sp)]
    d[ , brain_mass_with_sp_log10 := int_ +0.57*body_mass_with_sp_log10]
    d[ , brain_mass_with_sp := 10^brain_mass_with_sp_log10]

    ggplot(d, aes(x = body_mass_with_sp,  y =brain_mass_with_sp, col = factor(species), fill = factor(species)))+
      geom_abline(intercept=int_, slope=0.57, col = 'red', lty = 3) +
      geom_point(pch = 21, alpha = 0.6) +
      #geom_smooth(method = 'lm', se = FALSE) +
      annotate("text", x=1000, y=0.6, label= "Evolutionary allometry for Aves", col = 'red',  angle = 45) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')+
      theme(legend.position = 'none')    

  # add random noise along y axis
    d[ , brain_mass_with_sp_var := rnorm(n = 1 , mean = brain_mass_with_sp, sd =abs(brain_mass_with_sp*0.045)), by = pk]
    d[ , brain_mass_with_sp_var_log10 := log10(brain_mass_with_sp_var)]

  # add random noise along x axis
    d[ , body_mass_with_sp_var := rnorm(n = 1 , mean = body_mass_with_sp, sd =abs(body_mass_with_sp*0.08)), by = pk]
    d[ , body_mass_with_sp_var_log10 := log10(body_mass_with_sp_var)]  

    dc = data.table(species = unique(d$species))
    for(i in unique(d$species)){
      di = d[species == i]
      dc[species == i, OLS := coef(lm(log10(brain_mass_with_sp_var)~log10(body_mass_with_sp_var), data = di))[2]]
      dc[species == i, RMS :=  with(di, line.cis(log10(brain_mass_with_sp_var), log10(body_mass_with_sp_var)))[2,1]]  
      dc[species == i, RMS_int :=  with(di, line.cis(log10(brain_mass_with_sp_var), log10(body_mass_with_sp_var)))[1,1]]  
      dc[ species == i, min_body := min(di$body_mass_with_sp_var_log10)]  
      dc[ species == i, max_body := max(di$body_mass_with_sp_var_log10)]  
    }
    
    dc[ ,body_mass_with_sp_var := c(12, 120, 1200, 1800)] # coordinates for labels
    dc[ ,brain_mass_with_sp_var := c(0.2, 0.8, 3, 10)] # coordinates for labels

    dc[, label := paste0('OLS = ', round(OLS,2), '\nRMS = ', round(RMS,2))]

    dc[ ,min_brain := RMS_int+RMS*min_body] # need addin of the intercept
    dc[ ,max_brain := RMS_int+RMS*max_body]  

    g4=
    ggplot(d, aes(x = body_mass_with_sp_var,  y =brain_mass_with_sp_var, col = factor(species), fill = factor(species)))+
      geom_abline(intercept=int_, slope=0.57, col = 'red', lty = 3) +
      geom_point(pch = 21, alpha = 0.6) +
      geom_smooth(method = 'lm', se = FALSE, col = 'black') +
      geom_segment(data = dc, aes(x = 10^min_body, xend = 10^max_body, y = 10^min_brain, yend = 10^max_brain), col = 'darkgrey', lwd=1.5) +
      geom_text(data = dc, aes(x = body_mass_with_sp_var,  y =brain_mass_with_sp_var, col = factor(species), label = label), size = 3)+
      annotate("text", x=160, y=1.5, label= "Evolutionary allometry\nfor Aves = 0.57", col = 'red',  angle = 45, size = 3) + 
      scale_x_continuous(trans = 'log10', 'Body mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')+
      theme(legend.position = 'none')  
    g4
    d[, range(body_mass_with_sp_var), by = species] 
    d[, range(brain_mass_with_sp_var), by = species] 

    ggsave('Outputs/brain-allo-sim2mod.png', g4, width = 4, height =4)

# END