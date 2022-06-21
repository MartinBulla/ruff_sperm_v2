require(car)
require(data.table)
require(ggplot)
require(ggpubr)
require(here)
require(multcomp)
require(readxl)

b = data.table(read_excel(here::here('Data/Ruff brains organs collected 2019 for Suzana.xlsx'), sheet = 1))#, range = "A1:G161"))
b[, morph := factor(morph, levels=c("R", "S", "C"))] 
b[morph == 'R', morph := 'Independent']
b[morph == 'S', morph := 'Satellite']
b[morph == 'C', morph := 'Faeder']
b[, soma := bodyweight - brainweight]
b=b[sex == 'F']

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
  line_log10 = data.table(brainweight= 10^(mean(log10(b$brainweight))+c(-0.2,0,0.2)), soma = 10^(mean(log10(b$soma))+c(-0.2,0,0.2)), morph = 'Faeder')
  ggplot(b, aes(y = brainweight, x = soma, col = morph)) + geom_point() + 
      stat_smooth(method = 'lm') +
      stat_cor(aes(label = ..r.label..),  size = 3)+ #label.x = 1, 
      geom_line(data = line_log10, aes(y = brainweight, x = soma), col = 'black', lty = 3)+
      geom_point(y = mean(log10(b$brainweight)), x=mean(log10(b$soma)), col = 'black')+
      scale_x_continuous(trans = 'log10', 'Body - brain mass [g]')+
      scale_y_continuous(trans = 'log10','Brain mass [g]')+
      labs(tag = '♂')

  ggsave(file = 'Outputs/brain_body_allometry_log10_M.png', width = 6, height = 5)
