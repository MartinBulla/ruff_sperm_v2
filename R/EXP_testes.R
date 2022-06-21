require(car)
require(data.table)
require(ggplot)
require(ggpubr)
require(here)
require(multcomp)
require(readxl)

t = data.table(read_excel(here::here('Data/testes.xlsx'), sheet = 1))#, range = "A1:G161"))
t[, Morph := factor(Morph, levels=c("Res", "Sat", "Faed"))] 
t[Morph == 'Res', Morph := 'Independent']
t[Morph == 'Sat', Morph := 'Satellite']
t[Morph == 'Faed', Morph := 'Faeder']
t[, soma := Bodymass - Gonadmass]

# allometry
  
  leveneTest(log10(Gonadmass)~Morph,t) # ok

  ai=aov(log10(Gonadmass)~Morph+log10(soma) +Morph:log10(soma),t)
  Anova(ai, type="III")
  
  as=aov(log10(Gonadmass)~Morph+log10(soma),t)
  Anova(as, type="III")


  mi = lm(log10(Gonadmass)~Morph+log10(soma) + Morph:log10(soma), t)
  summary(glht(mi))
  summary((mi))

  ms = lm(log10(Gonadmass)~Morph+log10(soma), t)

  AIC(mi,ms)
# raw data
  #t = t[!Gonadmass<2] 
  leveneTest(Gonadmass~Morph,t) # ok

  ai=aov(Gonadmass~Morph+soma +Morph:soma,t)
  Anova(ai, type="III")
  
  as=aov(Gonadmass~Morph+soma,t)
  Anova(as, type="III")


  mi = lm(Gonadmass~Morph+soma + Morph:soma, t)
  summary(glht(mi))
  summary((mi))

  ms = lm(Gonadmass~Morph+soma, t)

  AIC(mi,ms)

  mil = lm(log(Gonadmass)~Morph+log(soma) + Morph:log(soma), t)
  msl = lm(log(Gonadmass)~Morph+log(soma), t)

  AIC(mil,msl)

# plot
  line_log = data.table(Gonadmass= exp(mean(log(t$Gonadmass))+c(-0.37,0,0.2)), soma = exp(mean(log(t$soma))+c(-0.37,0,0.2)), Morph = 'Faeder')
  line_log10 = data.table(Gonadmass= 10^(mean(log10(t$Gonadmass))+c(-0.15,0,0.1)), soma = 10^(mean(log10(t$soma))+c(-0.15,0,0.1)), morph = 'Faeder')

  ggplot(t, aes(y = Gonadmass, x = soma, col = Morph)) + geom_point() + 
      stat_smooth(method = 'lm') +
      stat_cor(aes(label = ..r.label..),  size = 3)+ #label.x = 1, 
      geom_line(data = line_log, aes(y = Gonadmass, x = soma), col = 'black', lty = 3)+
      geom_point(y = mean(log(t$Gonadmass)), x=mean(log(t$soma)), col = 'black')+
      scale_x_continuous(trans = 'log', 'Body - gonad mass [g]')+
      scale_y_continuous(trans = 'log','Gonad mass [g]')

  ggsave(file = 'Outputs/testes_body_allometry_ln.png', width = 6, height = 5)

  ggplot(t, aes(y = Gonadmass, x = soma, col = Morph)) + geom_point() + 
      stat_smooth(method = 'lm') +
      stat_cor(aes(label = ..r.label..),  size = 3)+ #label.x = 1, 
      geom_line(data = line_log10, aes(y = Gonadmass, x = soma), col = 'black', lty = 3)+
      geom_point(y = mean(log10(t$Gonadmass)), x=mean(log10(t$soma)), col = 'black')+
      scale_x_continuous(trans = 'log10', 'Body - gonad mass [g]')+
      scale_y_continuous(trans = 'log10','Gonad mass [g]')

  ggsave(file = 'Outputs/testes_body_allometry_log10.png', width = 6, height = 5)

# other  
  ggplot(t, aes(y = Gonadmass, x = soma, col = Morph)) + geom_point() + 
      stat_smooth(method = 'lm') +
      stat_cor(aes(label = ..r.label..),  label.x = 125, size = 3)+
      #geom_line(data = line_orig, aes(y = Gonadmass, x = soma), col = 'black', lty = 3)+
      #geom_point(x=mean(t$soma), y = mean(t$Gonadmass), col = 'black')+
      scale_x_continuous('Body - gonad mass [g]')+
      scale_y_continuous('Gonad mass [g]')    

  ggplot(t, aes(y = Gonadmass, x = SacDateDD, col = Morph)) + 
      stat_smooth(method = 'lm', se = FALSE) +
      geom_point() 

  ggplot(t, aes(y = Gonadmass, x = SacDateDD, col = Morph, shape = as.factor(SacDateYYYY))) + 
      stat_smooth(method = 'lm', se = FALSE) +
      geom_point() 
 
  ggplot(t, aes(y = soma, x = SacDateDD, col = Morph, shape = as.factor(SacDateYYYY))) + 
      stat_smooth(method = 'lm', se = FALSE) +
      geom_point() 