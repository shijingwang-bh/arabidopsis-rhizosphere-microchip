options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
args<-commandArgs(TRUE);pkg<-normalizePath(args[1],winslash='/');run<-normalizePath(args[2],winslash='/')
out<-file.path(run,'results/conditioned_medium')
suppressPackageStartupMessages({library(dplyr);library(jsonlite)})
raw<-read.csv(file.path(out,'data/growth_measurements_raw.csv'),check.names=FALSE)
raw$included_main<-as.logical(raw$included_main)
blank<-raw %>% filter(role!='assay') %>% group_by(collection_day,group,time_index,time_h) %>%
 summarise(blank_n=n(),blank_mean=mean(od_raw),blank_sd=sd(od_raw),.groups='drop') %>% rename(blank_group=group)
curves<-raw %>% filter(role=='assay') %>% mutate(blank_group=ifelse(group=='LB_control','LB_blank','MS_blank')) %>%
 left_join(blank,by=c('collection_day','blank_group','time_index','time_h')) %>%
 mutate(od_corrected=od_raw-blank_mean)
stopifnot(nrow(curves)==38*289,all(is.finite(curves$od_corrected)),sum(curves$included_main)==37*289)
cached<-read.csv(file.path(out,'data/workbook_cached_calculations.csv'))
bc<-cached %>% filter(kind=='blank_mean') %>% rename(blank_group=group) %>%
 left_join(blank,by=c('collection_day','blank_group','time_index'))
cc<-cached %>% filter(kind=='corrected') %>% select(-group) %>%
 left_join(curves,by=c('collection_day','well','time_index'))
errors<-c(blank_means=max(abs(bc$cached_value-bc$blank_mean)),corrected_day7=max(abs(cc$cached_value-cc$od_corrected)))
stopifnot(all(errors<1e-12))
summarize_curves<-function(x)x %>% group_by(collection_day,group,time_index,time_h) %>%
 summarise(n_wells=n(),mean=mean(od_corrected),sd=sd(od_corrected),se=sd/sqrt(n_wells),.groups='drop')
s<-summarize_curves(filter(curves,included_main))
sall<-summarize_curves(curves)
counts<-curves %>% group_by(collection_day,group) %>% summarise(source_wells=n_distinct(well),retained_wells=n_distinct(well[included_main]),time_points=n_distinct(time_index),.groups='drop')
metrics<-curves %>% arrange(time_h) %>% group_by(collection_day,group,well,included_main) %>%
 summarise(initial_OD=first(od_corrected),final_OD=last(od_corrected),max_OD=max(od_corrected),time_of_max_h=time_h[which.max(od_corrected)],.groups='drop')
for(n in c('raw','blank','curves','s','sall','counts','metrics')){
 nm<-c(raw='all_well_raw_OD',blank='background_timepoint_summary',curves='all_well_corrected_OD',s='main_curve_mean_SD_SE',sall='all_well_curve_mean_SD_SE',counts='well_counts',metrics='per_well_descriptive_summary')[[n]]
 write.csv(get(n),file.path(out,'tables',paste0(nm,'.csv')),row.names=FALSE,fileEncoding='UTF-8')
}
saveRDS(list(raw=raw,blank=blank,curves=curves,s=s,sall=sall,counts=counts),file.path(out,'data/growth_analysis.rds'))
write_json(list(cached_workbook_max_errors=as.list(errors),all_time_points=289,time_step_h=.25,time_range_h=c(0,72),
 source_assay_wells=38,retained_assay_wells=37,main_plot_wells=27,
 excluded_day7_E2='Uninoculated E2 well, retained in the archived individual-well display.',
 error_band='SE across retained assay wells; descriptive, not independent collection-batch inference.',
 no_smoothing=TRUE,no_time_point_deletion=TRUE,no_negative_clipping=TRUE),file.path(out,'qa/growth_analysis_checks.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'qa/sessionInfo_analysis.txt'))
print(counts);print(errors)
