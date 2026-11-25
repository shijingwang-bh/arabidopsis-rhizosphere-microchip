options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows') invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
pkg<-normalizePath(args[1],winslash='/');run<-normalizePath(args[2],winslash='/')
suppressPackageStartupMessages({library(car);library(jsonlite)})
out<-file.path(run,'results/phenotypes');dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_input<-function(n)read.csv(file.path(pkg,'data/observations',n),check.names=FALSE)
f1<-read_input('phenotype_pre.csv');fl<-read_input('rfp.csv')
cc<-read_input('bacterial_counts.csv');f3<-read_input('growth_day7.csv')
save_table<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,fileEncoding='UTF-8')
glev <- c('stationary_1g','clinorotation')
m1 <- c('Lateral root number','Root length (mm)','Leaf area (cm\u00b2)','Leaf fresh weight (mg)')
m3 <- c('Root elongation (mm)','Root diameter (\u03bcm)','Leaf fresh weight (mg)','Leaf area (cm\u00b2)')
for(z in c('f1','fl','cc','f3')){x<-get(z);x$gravity<-factor(x$gravity,levels=glev);assign(z,x)}
f3$bacteria<-factor(f3$bacteria,levels=c('uninoculated','inoculated'))
stopifnot(all(is.finite(f1$value)),all(is.finite(f3$value)),all(is.finite(fl$value)),all(is.finite(cc$cells_1e8_per_ml)))
welch <- function(x,col){
 a<-x[x$gravity=='clinorotation',col];b<-x[x$gravity=='stationary_1g',col];t<-t.test(a,b)
 data.frame(n_clinorotation=length(a),n_stationary=length(b),mean_clinorotation=mean(a),mean_stationary=mean(b),mean_difference=mean(a)-mean(b),ci95_low=t$conf.int[1],ci95_high=t$conf.int[2],p_value=t$p.value)
}
f1stats<-do.call(rbind,lapply(m1,function(m)cbind(metric=m,welch(f1[f1$metric==m,],'value'))));f1stats$fdr<-p.adjust(f1stats$p_value,'BH')
cstats<-do.call(rbind,lapply(c(3,7),function(d)cbind(day=d,welch(cc[cc$day==d,],'cells_1e8_per_ml'))));cstats$fdr<-p.adjust(cstats$p_value,'BH')
astats<-list();hstats<-list();effects<-list()
for(m in m3){
 x<-f3[f3$metric==m,];fit<-lm(value~gravity*bacteria,x)
 a<-as.data.frame(car::Anova(fit,type=2));h<-as.data.frame(car::Anova(fit,type=2,white.adjust='hc3'))
 astats[[m]]<-data.frame(metric=m,term=rownames(a),df=a$Df,f_value=a[['F value']],p_value=a[['Pr(>F)']])
 hstats[[m]]<-data.frame(metric=m,term=rownames(h),df=h$Df,f_value=h[['F']],p_value=h[['Pr(>F)']])
 L<-c(0,0,1,.5);est<-sum(L*coef(fit));se<-sqrt(drop(t(L)%*%vcov(fit)%*%L));df<-df.residual(fit)
 effects[[m]]<-data.frame(metric=m,estimate=est,standard_error=se,df=df,ci95_low=est-qt(.975,df)*se,ci95_high=est+qt(.975,df)*se,p_value=2*pt(-abs(est/se),df))
}
a<-do.call(rbind,astats);h<-do.call(rbind,hstats);e<-do.call(rbind,effects);rownames(a)<-rownames(h)<-rownames(e)<-NULL
summarize <- function(x,keys,col='value'){
 ss<-split(x,interaction(x[keys],drop=TRUE));do.call(rbind,lapply(ss,function(v)cbind(v[1,keys,drop=FALSE],n=nrow(v),mean=mean(v[[col]]),se=sd(v[[col]])/sqrt(nrow(v)))))
}
s1<-summarize(f1,c('metric','gravity'));s3<-summarize(f3,c('metric','gravity','bacteria'));sf<-summarize(fl,c('gravity','region','day'));sc<-summarize(cc,c('gravity','day'),'cells_1e8_per_ml')
stopifnot(nrow(sf)==32,all(sf$n==3),nrow(fl)==96,nrow(f1)==81,nrow(f3)==90,nrow(cc)==17)

# November 25 primary rule: all FOUR growth contrasts form one BH family.
e$q_value<-p.adjust(e$p_value,'BH')
cc_model<-lm(cells_1e8_per_ml~gravity*factor(day),cc)
cc_type2<-as.data.frame(car::Anova(cc_model,type=2))
cc_type2<-data.frame(term=rownames(cc_type2),df=cc_type2$Df,F=cc_type2[['F value']],P=cc_type2[['Pr(>F)']])
save_table(cc_type2,'Bacterial_concentration_TypeII')
for(n in c('f1stats','cstats','a','h','e','s1','s3','sf','sc'))save_table(get(n),switch(n,
 f1stats='Figure1d_Welch_effects',cstats='Figure3c_daywise_Welch_effects',
 a='Growth_traits_TypeII',h='Growth_traits_HC3',e='Growth_traits_average_inoculation_effects',
 s1='Figure1d_group_summaries',s3='Growth_traits_group_summaries',
 sf='Figure3b_fluorescence_summaries',sc='Figure3c_cell_summaries'))
stopifnot(abs(e$q_value[e$metric=='Leaf fresh weight (mg)']-0.00226342911804885)<1e-12)
stopifnot(abs(e$q_value[e$metric=='Root diameter (μm)']-0.00408478329865372)<1e-12)
stopifnot(abs(cc_type2$P[cc_type2$term=='gravity:factor(day)']-0.160220503684902)<1e-10)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('PASS: phenotype and bacterial concentration models; four-trait joint BH correction\n')
