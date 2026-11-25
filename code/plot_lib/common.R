options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
suppressPackageStartupMessages({library(ggplot2);library(grid);library(patchwork);library(magick);library(ggrepel);library(circlize)})
blue<-'#0072B2';orange<-'#D55E00';pal<-c(stationary_1g=blue,clinorotation=orange);font<-'Arial'
pal_fill<-c(stationary_1g='#C8E1EF',clinorotation='#F4D7C2');pal_dark<-c(stationary_1g='#004F73',clinorotation='#8D3E00')
glab<-c(stationary_1g='1 g',clinorotation='s\u03bcg')
mm<-function(n)unit(n,'mm')
rd<-function(fig,n,...)read.csv(file.path(plot_data,aliases[[fig]],paste0(n,'.csv')),check.names=FALSE,...)
write_source<-function(...) invisible(NULL)
qfmt<-function(q)ifelse(q<.001,formatC(q,format='g',digits=2),ifelse(q<.01,sprintf('%.4f',q),sprintf('%.3f',q)))
theme_set(theme_classic(base_size=10,base_family=font)+theme(axis.text=element_text(colour='black',size=9),axis.title=element_text(size=10),plot.title=element_text(size=10,face='bold'),plot.subtitle=element_text(size=9),strip.background=element_blank(),strip.text=element_text(size=10,face='bold'),legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=9),legend.key.height=mm(3),legend.key.width=mm(6),plot.margin=margin(3,5,3,3),panel.grid=element_blank()))
txt<-function(s,x,y,size=10,bold=FALSE,col='black',just='left')grid.text(s,x=mm(x),y=unit(1,'npc')-mm(y),just=c(just,'top'),gp=gpar(fontfamily=font,fontsize=size,fontface=if(bold)'bold'else'plain',col=col,lineheight=1.1))
place<-function(x,y,w,h,fun){pushViewport(viewport(x=mm(x),y=unit(1,'npc')-mm(y),width=mm(w),height=mm(h),just=c('left','top')));fun();popViewport()}
drawp<-function(p)grid.draw(if(inherits(p,'patchwork'))patchworkGrob(p)else ggplotGrob(p))
header<-function(letter,title){txt(letter,0,0,12,TRUE);txt(title,6,0,10.5,TRUE)}
panel<-function(p,letter,title,w,h){header(letter,title);place(0,10,w,h-10,function()drawp(p))}
image_at<-function(path,x,y,w,h,crop=NULL){im<-image_read(path)[1];if(!is.null(crop))im<-image_crop(im,crop,repage=TRUE);place(x,y,w,h,function()grid.raster(as.raster(im),interpolate=TRUE))}
export<-function(...) invisible(NULL)

offset<-function(d,keys,base,spread=.07){d$x<-base;for(ii in split(seq_len(nrow(d)),interaction(d[keys],drop=TRUE))){d$x[ii]<-base[ii]+if(length(ii)>1)seq(-spread,spread,length.out=length(ii))else 0};d}
bracket<-function(p,x1,x2,y,label,drop){p+annotate('segment',x=x1,xend=x2,y=y,yend=y,linewidth=.35)+annotate('segment',x=c(x1,x2),xend=c(x1,x2),y=y,yend=y-drop,linewidth=.35)+annotate('text',x=(x1+x2)/2,y=y+drop*.5,label=label,vjust=0,size=3,family=font)}
