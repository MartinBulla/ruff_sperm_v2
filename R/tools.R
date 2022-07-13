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
    packages = c('data.table', 'ggplot2', 'ggthemes','glue','grid', 'gridExtra','htmlTable', 'lattice', 'lubridate', 'magrittr', 'maptools', 'plyr','raster','RColorBrewer','readxl','scales','scatterplot3d','stringr','zoo','viridis','writexl')
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
    cols = c('black','darkgrey','#bf925a')
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
      
# END