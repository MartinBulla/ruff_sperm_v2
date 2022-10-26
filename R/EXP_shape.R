# TOOLS
  require(here)
  source(here::here('R/tools.R'))

# DATA
  source(here::here('R/DAT_prepare.R'))       

# relative values
  bw[,h_rel:=Head/Total]
  bw[,n_rel_h:=Nucleus/Head]
  bw[,t_rel:=Tail/Total]
  bw[,h_rel_tail:=Head/Tail]
  bw[,h_rel_flag:=Head/Flagellum]
  bw[,t_rel_flag:=Tail/Flagellum]

  summary(bw)
  ggplot(bw, aes(x=h_rel_tail, y = h_rel_flag)) + geom_point()
  ggplot(bw, aes(x=h_rel_tail, y = Flagellum_rel)) + geom_point()
 
 # correlations
    cor((bw$Flagellum), (bw$Midpiece))
    cor(log(bw$Flagellum), log(bw$Midpiece))
    cor(log(bw$Tail), log(bw$Midpiece))
    cor(log(bw$Tail), log(bw$Head))
    cor((bw$Tail), (bw$Head))
 # END