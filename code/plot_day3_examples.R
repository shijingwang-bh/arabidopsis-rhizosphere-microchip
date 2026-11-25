# Numeric redraw of Figure S7 from the complete archived contrast family and
# original provider intensities. The submitted raster appearance is archived
# separately because its original vector drawing source was not available.
meta<-read.csv(file.path(pkg,'data/lcms/corrected_sample_metadata_with_QC.csv'))
meta<-meta[meta$gravity!='QC'&meta$day==3,]
lines<-readLines(file.path(pkg,'data/lcms/metabolites_exp.xls'),warn=FALSE,encoding='bytes')
lines<-iconv(lines,from='GB18030',to='UTF-8',sub='byte')
raw<-read.delim(text=paste(lines,collapse='\n'),check.names=FALSE,quote='',comment.char='')
names(raw)[1]<-'feature_id'
stats<-rd('Figure_4','feature_contrasts_primary_named_QC30')
ids<-c('neg_5318','pos_12437')
labels<-c('9,10-Dihydroxy-12,13-\nepoxyoctadecanoic acid','9,12-Octadecadiynoic acid')
plots<-lapply(seq_along(ids),function(k){
 id<-ids[k];row<-match(id,raw$feature_id);stopifnot(!is.na(row))
 d<-meta;d$value<-log2(as.numeric(raw[row,meta$sample])+1)
 d$x<-ifelse(d$gravity=='stationary_1g',0,3)+ifelse(d$bacteria=='uninoculated',1,2)
 s<-do.call(rbind,lapply(split(d,d$x),function(a){
  m<-mean(a$value);se<-sd(a$value)/sqrt(nrow(a));ci<-qt(.975,nrow(a)-1)*se
  data.frame(x=a$x[1],mean=m,low=m-ci,high=m+ci)
 }))
 pair<-stats[stats$feature_id==id&stats$contrast%in%c('inoculation_stationary_day3','inoculation_clinorotation_day3'),]
 j<-stats[stats$feature_id==id&stats$contrast=='interaction_day3',]
 p<-ggplot(d,aes(x,value,colour=gravity))+
  geom_point(position=position_jitter(width=.085,height=0,seed=20261125),size=2)+
  geom_errorbar(data=s,aes(x=x,ymin=low,ymax=high),inherit.aes=FALSE,width=.15,colour='black',linewidth=.5)+
  geom_segment(data=s,aes(x=x-.18,xend=x+.18,y=mean,yend=mean),inherit.aes=FALSE,linewidth=.8)+
  scale_colour_manual(values=pal,guide='none')+
  scale_x_continuous(breaks=c(1,2,4,5),labels=c('1 g\nU','1 g\nI','s\u03bcg\nU','s\u03bcg\nI'),limits=c(.5,5.5))+
  scale_y_continuous(expand=expansion(mult=c(.05,.20)))+
  labs(x=NULL,y=expression(log[2](intensity+1)),title=labels[k],subtitle=sprintf('Gravity × inoculation: q = %.3f',j$q))+
  theme(plot.title=element_text(size=9.5,hjust=.5),plot.subtitle=element_text(size=9,hjust=.5))
 span<-diff(range(c(d$value,s$low,s$high)))
 for(g in c('stationary','clinorotation')){
  q<-pair$q[pair$contrast==paste0('inoculation_',g,'_day3')]
  a<-if(g=='stationary')1 else 4
  yy<-max(d$value[d$x%in%c(a,a+1)],s$high[s$x%in%c(a,a+1)])+.07*span
  p<-bracket(p,a,a+1,yy,paste0('q = ',qfmt(q)),span*.03)
 }
 p
})
save_panels('Figure_S7',190,116,function(){
 place(0,0,94,116,function()panel(plots[[1]],'a','',94,116))
 place(98,0,92,116,function()panel(plots[[2]],'b','',92,116))
})
