
#' Clustering
  require('apcluster') 
  require('Rcpp')                                                                                                                       
  bw[, ID := .I]
  aw[, ID := .I]

  get_clusters <- function(apc) {
    x <- lapply(apc@clusters, function(x) data.table(ID = names(x) %>% as.numeric()))
    for (i in 1:length(x)) x[[i]][, clust_ID := i]
    rbindlist(x)
  }

#  scale variables
    aw[, Midpiece_z := scale(Midpiece)]
    aw[, Nucleus_z := scale(Nucleus)]
    aw[, Tail_z := scale(Tail)]
    aw[, VSL_z := scale(VSL)]
#  all data
    apc = apclusterK
    apc <- apcluster(corSimMat(r = 1), bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)
    apc = apclusterK(corSimMat(r = 1), bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)
    
    apc <- apcluster(corSimMat(r = 1), bw[, .(Head, Midpiece, Tail, Total)], q = 0 ,  details = TRUE)
    apc = apclusterK(corSimMat(r = 1), bw[, .(Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total)],  K = 3,  details = TRUE)


    apc <- apcluster(negDistMat(r=2) , bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)
    apc = apclusterK(negDistMat(r = 2), bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)
    apc <- apcluster(expSimMat(r=2) , bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)
    #apc = apclusterK(expSimMat(r = 2), bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)
    apc <- apcluster(linKernel() , bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)

    cc <- get_clusters(apc)

    O <- merge(cc, bw, by = "ID")

    #O[, .N, .(clust_ID, Morph)]

    #dev.new()
    ggplot(O, aes(x = clust_ID, y = Morph)) +
    geom_count()

    plot(apc, bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)])
    heatmap(apc)

#  averages   
    #apc <- apcluster(corSimMat(r = 1), aw[, .(Midpiece)], q = 0 ,  details = TRUE)
    #apc <- apclusterK(corSimMat(r = 1), aw[, .(Midpiece)], K = 3 ,  details = TRUE)

    apc <- apcluster(corSimMat(r = 1), aw[, .(Nucleus_z, Midpiece_z)], q = 0 ,  details = TRUE)
    apc <- apcluster(corSimMat(r = 1), aw[, .(Nucleus, Midpiece)], q = 0 ,  details = TRUE)
    apc <- apcluster(corSimMat(r = 1), aw[, .(Nucleus, Midpiece)], maxits = 2000, details = TRUE)
    
    apc <- apclusterK(corSimMat(r = 1), aw[, .(Nucleus_z, Midpiece_z)], K = 3 ,  bimaxit = 60, details = TRUE)

    apc <- apclusterK(corSimMat(r = 1), aw[, .(Nucleus, Midpiece)], K = 3 ,  bimaxit = 40, details = TRUE)

    apc <- apcluster(corSimMat(r = 1), aw[, .(Nucleus, Midpiece, Tail)], q = 0 ,  details = TRUE)
    apc <- apcluster(corSimMat(r = 1), aw[, .(Nucleus_z, Midpiece_z, Tail_z)], q = 0 ,  details = TRUE)
    apc <- apcluster(corSimMat(r = 1), aw[, .(Nucleus, Midpiece, Tail)],   details = TRUE)
    apc <- apclusterK(corSimMat(r = 1), aw[, .(Nucleus_z, Midpiece_z, Tail_z)], K = 3 ,  details = TRUE)
      (19+26+4)/92
    apc <- apclusterK(corSimMat(r = 1), aw[, .(Nucleus, Midpiece, Tail)], K = 3 ,  details = TRUE)

    apc <- apcluster(corSimMat(r = 1), aw[, .(Nucleus, Midpiece, Tail, VSL)], q = 0 ,  details = TRUE)
    apc <- apcluster(corSimMat(r = 1), aw[, .(Nucleus_z, Midpiece_z, Tail_z, VSL_z)], q = 0 ,  details = TRUE)
      (15+9+5)/92
    apc <- apclusterK(corSimMat(r = 1), aw[, .(Nucleus, Midpiece, Tail, VSL)], K = 3 ,  details = TRUE)
    apc <- apclusterK(corSimMat(r = 1), aw[, .(Nucleus_z, Midpiece_z, Tail_z, VSL_z)], K = 3 ,  details = TRUE)
       (24+14+6)/92

    apc <- apcluster(negDistMat(r = 1), aw[, .(Nucleus, Midpiece, Tail, VSL)], q = 0 ,  details = TRUE)
    apc <- apclusterK(negDistMat(r = 1), aw[, .(Nucleus_z, Midpiece_z, Tail_z, VSL_z)], K = 3 ,  details = TRUE)

    apc <- apcluster(negDistMat(r = 2), aw[, .(Nucleus_z, Midpiece_z, Tail_z, VSL_z)], q = 0 ,  details = TRUE)
    apc <- apclusterK(negDistMat(r = 2), aw[, .(Nucleus_z, Midpiece_z, Tail_z, VSL_z)], K = 3 ,  details = TRUE)

    apc <- apcluster(expSimMat(r = 2), aw[, .(Nucleus, Midpiece, Tail, VSL)], q = 0 ,  details = TRUE)
    apc <- apclusterK(expSimMat(r = 2), aw[, .(Nucleus_z, Midpiece_z, Tail_z, VSL_z)], K = 3 ,  details = TRUE)


    apc <- apcluster(expSimMat(r = 2), aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)
    apc = apclusterK(expSimMat(r = 2), aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)
    
    

    apc <- apcluster(corSimMat(r = 1), aw[, .(Head, Midpiece, Tail, Total)], q = 0 ,  details = TRUE)
    apc <- apclusterK(corSimMat(r = 1), aw[, .(Head, Midpiece, Tail, Total)], K = 3 ,  details = TRUE)

    apc <- apcluster(corSimMat(r = 1), aw[, .(VAP, VSL, VCL)], q = 0 ,  details = TRUE)
    apc <- apclusterK(corSimMat(r = 1), aw[, .(VAP, VSL, VCL)], K = 3 ,  details = TRUE)

    apc <- apcluster(negDistMat(r=2) , aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0.5 ,  details = TRUE)
     apc = apclusterK(negDistMat(r = 2), aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)
    
    apc <- apcluster(expSimMat(r=2) , aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)
     
     apc = apclusterK(expSimMat(r = 2), aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)

    cc <- get_clusters(apc)

    O <- merge(cc, aw, by = "ID")

    O[, .N, .(clust_ID, Morph)]

    #dev.new()
    ggplot(O, aes(x = clust_ID, y = Morph)) + geom_count()

    plot(apc, aw[, .(Nucleus, Midpiece, Tail)])
    plot(apc, aw[, .(Head, Midpiece, Tail, Total)])
    heatmap(apc)