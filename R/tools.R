 # =============================================================
 # ❗ Runs relative to the project's root directory and
 # loads packages, functions and pictures of the morhps
 # =============================================================

#' Loads packages and installs those that are not in the library
#' @param  vector of package names
#' @export

using<-function(...) {
    libs<-unlist(list(...))
    req<-unlist(lapply(libs,require,character.only=TRUE))#, quietly  = TRUE
    need<-libs[req==FALSE]
    if(length(need)>0){ 
        install.packages(need)
        lapply(need,require,character.only=TRUE)
    }
}

# load/install packages
    packages = c('data.table', 'ggplot2', 'ggthemes','glue','grid', 'gridExtra','htmlTable', 'lattice', 'lubridate', 'magrittr', 'maptools', 'plyr','png','raster','RColorBrewer','readxl','scales','scatterplot3d','stringr','zoo','viridis','writexl')
  sapply(packages, function(x) suppressPackageStartupMessages(using(x)) )

# constants
    set.seed = 5
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    fae = '#d4b691' # 'ffd6af'
    sat = 'white'
    ind = '#303030'
    colors = c(ind,sat,fae)
    fills = c(ind,sat,fae)
    cols = c('black','grey30','#bf925a') #grey30 instead of darkgrey
    fill_zf ='#a53708' 
    col_zf = '#f89f79'
 
# Set system time
   Sys.setenv(TZ="UTC")

# Customized ggplot theme
    theme_MB = theme(  
              axis.line = element_blank(),
              #axis.line = element_line(colour="grey70", size=0.25),
              axis.title = element_text(size=7.25, colour="grey10"),
              axis.title.y = element_text(vjust=1.5),
              axis.title.x = element_text(vjust=1),
              axis.text = element_text(size=6.25),#, vjust = 0.5, hjust=1),# margin=units(0.5,"mm")),
              axis.ticks.length=unit(0.5,"mm"),
              axis.ticks = element_line(colour = "grey70", size = 0.1),
              #axis.ticks.margin,
              
              strip.text.x = element_text(size = 6, color="grey20",  margin=margin(1,1,1,1,"mm")), #grey50
              strip.text.y = element_text(size = 6, color="grey20",  margin=margin(1,1,1,1,"mm")), #grey50
              strip.background = element_rect(fill="grey99",colour="grey70", size=0.25),
                #strip.background = element_blank(), 
                #strip.text = element_blank(),
              panel.spacing = unit(0, "mm"),
              panel.background=element_blank(),
              panel.border = element_rect(colour="grey40", size=0.1, fill = NA), #panel.border=element_blank(),
              #panel.grid = element_blank(),
              #panel.grid = element_line(colour = "grey92", size = 0.1), 
              #panel.grid.minor = element_line(size = rel(0.5)),
              
              legend.text=element_text(size=6),
              legend.title=element_text(size=6),
              legend.key = element_rect(colour = NA, fill = NA),
              legend.key.height= unit(0.5,"line"),
              legend.key.width = unit(0.25, "cm"),
              legend.margin = margin(0,0,0,0, unit="cm"),
              legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
              legend.background = element_blank()
              )

# Functions
  getime = function (x) {ifelse(is.na(x), as.numeric(NA), as.numeric(difftime(x, trunc(x,"day"), units = "hours")))}
      
  getDay = function (x) {as.Date(trunc(x, "day"))}
  
  # for adding single images to single panels in ggplot
    annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
        layer(data = data, stat = StatIdentity, position = PositionIdentity, 
            geom = ggplot2:::GeomCustomAnn,
            inherit.aes = TRUE, params = list(grob = grob, 
                                              xmin = xmin, xmax = xmax, 
                                              ymin = ymin, ymax = ymax))
      }
  # for standardized model outputs
    m_out = function(model = m, type = "mixed", 
        name = "define", dep = "define", fam = 'Gaussian',
        round_ = 3, nsim = 5000, aic = FALSE, save_sim = 'Data/model_sim/', back_tran = FALSE, perc_ = 1){
          # perc_ 1 = proportion or 100%
        bsim = sim(model, n.sim=nsim)  
        
        if(save_sim!=FALSE){save(bsim, file = paste0(save_sim, name,'.RData'))}
       
        if(type != "mixed"){
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 

          if(back_tran == TRUE & fam == "binomial"){
           v = perc_*plogis(v)
           ci = perc_*plogis(ci)
           }
          if(back_tran == TRUE & fam == "binomial_logExp"){
                v = perc_*(1-plogis(v))
                ci = perc_*(1-plogis(ci))
                ci = rbind(ci[2,],ci[1,])
               }

          if(back_tran == TRUE & fam == "Poisson"){
           v = exp(v)
           ci = exp(ci)
          }

         oi=data.frame(type='fixed',effect=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,])
          rownames(oi) = NULL
          oi$estimate_r=round(oi$estimate,round_)
          oi$lwr_r=round(oi$lwr,round_)
          oi$upr_r=round(oi$upr,round_)
          if(perc_ == 100){
           oi$estimate_r = paste0(oi$estimate_r,"%")
           oi$lwr_r = paste0(oi$lwr_r,"%")
           oi$upr_r = paste0(oi$upr_r,"%")
          }
         x=data.table(oi[c('type',"effect", "estimate_r","lwr_r",'upr_r')]) 
       
        }else{
         v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
         ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 

         if(back_tran == TRUE & fam == "binomial"){
          v = perc_*plogis(v)
          ci = perc_*plogis(ci)
         }
          if(back_tran == TRUE & fam == "binomial_logExp"){
                v = perc_*(1-plogis(v))
                ci = perc_*(1-plogis(ci))
                ci = rbind(ci[2,],ci[1,])
               }

          if(back_tran == TRUE & fam == "Poisson"){
            v = exp(v)
            ci = exp(ci)
         }

        oi=data.table(type='fixed',effect=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,])
            rownames(oi) = NULL
            oi[,estimate_r := round(estimate,round_)]
            oi[,lwr_r := round(lwr,round_)]
            oi[,upr_r :=round(upr,round_)]
            if(perc_ == 100){
             oi[,estimate_r := paste0(estimate_r,"%")]
             oi[,lwr_r := paste0(lwr_r,"%")]
             oi[,upr_r := paste0(upr_r,"%")]
            }
         oii=oi[,c('type',"effect", "estimate_r","lwr_r",'upr_r')] 
        
         l=data.frame(summary(model)$varcor)
         l = l[is.na(l$var2),]
         l$var1 = ifelse(is.na(l$var1),"",l$var1)
         l$pred = paste(l$grp,l$var1)

         q050={}
         q025={}
         q975={}
         pred={}
         
         # variance of random effects
         for (ran in names(bsim@ranef)) {
           ran_type = l$var1[l$grp == ran]
           for(i in ran_type){
            q050=c(q050,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.5)))
            q025=c(q025,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.025)))
            q975=c(q975,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.975)))
            pred= c(pred,paste(ran, i))
            }
           }
         # residual variance
         q050=c(q050,quantile(bsim@sigma^2, prob=c(0.5)))
         q025=c(q025,quantile(bsim@sigma^2, prob=c(0.025)))
         q975=c(q975,quantile(bsim@sigma^2, prob=c(0.975)))
         pred= c(pred,'Residual')

         ri=data.table(type='random %',effect=pred, estimate_r=round(100*q050/sum(q050)), lwr_r=round(100*q025/sum(q025)), upr_r=round(100*q975/sum(q975)))
           
         ri[lwr_r>upr_r, lwr_rt := upr_r]
         ri[lwr_r>upr_r, upr_rt := lwr_r]
         ri[!is.na(lwr_rt), lwr_r := lwr_rt]
         ri[!is.na(upr_rt), upr_r := upr_rt]
         ri$lwr_rt = ri$upr_rt = NULL

         ri[,estimate_r := paste0(estimate_r,'%')]
         ri[,lwr_r := paste0(lwr_r,'%')]
         ri[,upr_r := paste0(upr_r,'%')]
        
        x = data.table(rbind(oii,ri))
        }
        
        x[1, model := name]                                                                
        x[1, response := dep]                                                                
        x[1, error_structure := fam]      
        N = length(resid(model))                                                          
        x[1, N := N ]                                                                

        x=x[ , c('model', 'response', 'error_structure', 'N', 'type',"effect", "estimate_r","lwr_r",'upr_r')] 

        if (aic == TRUE){   
            x[1, AIC := AIC(update(model,REML = FALSE))] 
            }
        if (aic == "AICc"){
            aicc = AICc(model)
            x[1, AICc := aicc] 
        }
        if(type == "mixed" & nrow(x[type=='random %' & estimate_r =='0%'])==0){
          R2_mar = as.numeric(r2_nakagawa(model)$R2_marginal)
          R2_con = as.numeric(r2_nakagawa(model)$R2_conditional)
          x[1, R2_mar := R2_mar]
          x[1, R2_con := R2_con]
         }
        x[is.na(x)] = ""
        return(x)
      } 
# model ass
m_ass_s = function(name = 'define', title = 'define', binomial = FALSE, mo = m0, dat = d, fixed = NULL, categ = NULL, trans = NULL, spatial = TRUE, temporal = TRUE, PNG = TRUE, outdir = 'outdir'){
    # binomial - shall a plot visualizing response means per sequence of fitted data be visualized?
    # trans - vector containing transformation function used to transform each predictor
   if(PNG == TRUE){
    png(paste(outdir,name, ".png", sep=""), width=6,height=9,units="in",res=600)
     }else{dev.new(width=6,height=9)}
   
   n = length(fixed)+length(categ) + 4 + if(temporal==TRUE){1}else{0} + if(spatial==TRUE){1}else{0} 
   par(mfrow=c(ceiling(n/3),3))
   
   scatter.smooth(fitted(mo),resid(mo),col='grey');abline(h=0, lty=2, col ='red')
   scatter.smooth(fitted(mo),sqrt(abs(resid(mo))), col='grey')
   if(binomial == TRUE){
      plot(fitted(mo), jitter(mo$model[,1], amount=0.05), xlab="Fitted values", ylab="Original values", las=1, cex.lab=1, cex=0.8,  main=list(paste("Probability of", names(mo$model)[1]),cex=0.8) )
      abline(0,1, lty=3)
      t.breaks <- cut(fitted(mo), quantile(fitted(mo)))
      means <- tapply(mo$model[,1], t.breaks, mean)
      semean <- function(x) sd(x)/sqrt(length(x))
      means.se <- tapply(mo$model[,1], t.breaks, semean)
      points(quantile(fitted(mo),c(0.125,0.375,0.625,0.875)), means, pch=16, col="orange")
      segments(quantile(fitted(mo),c(0.125,0.375,0.625,0.875)), means-2*means.se, quantile(fitted(mo),c(0.125,0.375,0.625,0.875)), means+2*means.se,lwd=2, col="orange")
   }

   qqnorm(resid(mo), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='grey');qqline(resid(mo))
  
   # variables
     scatter={} 
     for (i in rownames(summary(mo)$coef)) {
          #i = "lat_abs"
        j=sub("\\).*", "", sub(".*\\(", "",i)) 
        scatter[length(scatter)+1]=j
      }
      x = data.frame(scatter=unique(scatter)[2:length(unique(scatter))],
                      log_ = grepl("log",rownames(summary(mo)$coef)[2:length(unique(scatter))]), stringsAsFactors = FALSE)
      if(length(fixed)>0){
      for (i in 1:length(fixed)){
          jj =fixed[i]
          variable=dat[, ..jj][[1]]
          if(trans[i]=='log'){
          scatter.smooth(resid(mo)~log(variable),xlab=paste('log(',jj,')',sep=''), col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
          }else if(trans[i]=='abs'){
          scatter.smooth(resid(mo)~abs(variable),xlab=paste('abs(',jj,')',sep=''), col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
          }else if(trans[i]=='sin'){
            scatter.smooth(resid(mo)~sin(variable),xlab=paste('sin(',jj,')',sep=''), col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
          }else if(trans[i]=='cos'){
            scatter.smooth(resid(mo)~cos(variable),xlab=paste('cos(',jj,')',sep=''), col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
          } else {  
          scatter.smooth(resid(mo)~variable,xlab=jj,col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
        }
       }
      }
      if(length(categ)>0){
        for(i in categ){
           variable=dat[, ..i][[1]]
            boxplot(resid(mo)~variable, medcol='grey', whiskcol='grey', staplecol='grey', boxcol='grey', outcol='grey');abline(h=0, lty=3, lwd=1, col = 'red')
           }
      }     
          
   if(temporal == TRUE){
        acf(resid(mo), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
        }
   if(spatial == TRUE){    
      spdata=data.frame(resid=resid(mo), x=dat$Longitude, y=dat$Latitude)
        spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
        #cex_=c(1,2,3,3.5,4)
        cex_=c(1,1.5,2,2.5,3)
        spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
      plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
        legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
      plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals (<0)', cex=0.8))
      plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals (>=0)', cex=0.8))
        }
   
   mtext(title, side = 3, line = -1, cex=0.7,outer = TRUE)
   
   if(PNG==TRUE){dev.off()}
  }  

# prepare images
  img_i=readPNG('Illustrations/independent.png')
  img_s=readPNG('Illustrations/satelite.png')
  img_f=readPNG('Illustrations/faeder_crop.png')
  img_fc=readPNG('Illustrations/faeder_crop.png')
  img_z=readPNG('Illustrations/ZebraFinch_AusMale_crop.png')
  gi <- rasterGrob(img_i, interpolate=TRUE)
  gs <- rasterGrob(img_s, interpolate=TRUE)
  gf <- rasterGrob(img_f, interpolate=TRUE)
  gfc <- rasterGrob(img_fc, interpolate=TRUE)      
  gz <- rasterGrob(img_z, interpolate=TRUE)      
# END