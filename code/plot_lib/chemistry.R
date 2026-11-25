draw_Figure_4<-function(){
 d<-rd('Figure_4','PCA_biological_scores');pv<-rd('Figure_4','PCA_variance_explained');d$day<-factor(d$day,levels=c(3,7),labels=c('D3','D7'))
 pa<-ggplot(d,aes(PC1,PC2,colour=gravity,shape=bacteria,fill=day))+geom_point(size=2.7,stroke=.8)+scale_colour_manual(values=pal,labels=glab,breaks=names(pal))+scale_shape_manual(values=c(uninoculated=21,inoculated=24),labels=c(uninoculated='U',inoculated='I'),breaks=c('uninoculated','inoculated'))+scale_fill_manual(values=c(D3='white',D7='#929292'))+coord_equal()+labs(x=sprintf('PC1 (%.1f%%)',100*pv$variance_fraction[1]),y=sprintf('PC2 (%.1f%%)',100*pv$variance_fraction[2]))+guides(colour=guide_legend(order=1,nrow=1,override.aes=list(shape=21,fill='white')),shape=guide_legend(order=2,nrow=1,override.aes=list(colour='#555555',fill='white')),fill=guide_legend(order=3,nrow=1,override.aes=list(shape=21,colour='#555555')))+theme(legend.box='vertical',legend.text=element_text(size=8),legend.spacing.y=mm(0),legend.margin=margin(0,0,0,0))
 g<-rd('Figure_4','global_factorial_effects_family_sensitivity');g<-g[g$family=='primary_named_QC30',];g<-g[match(c('G','I','D','GI','GD','ID','GID'),g$term),];g$y<-7:1;gx<-max(g$partial_R2)*100;qx<-gx*1.50
 pb<-ggplot(g,aes(y=y))+geom_segment(aes(x=0,xend=100*partial_R2,yend=y),colour='#BFC9CD',linewidth=.7)+geom_point(aes(x=100*partial_R2),size=2.7,colour='#365C70')+geom_text(aes(x=100*partial_R2+gx*.06,label=sprintf('%.1f',100*partial_R2)),hjust=0,size=2.8)+geom_text(aes(x=qx,label=ifelse(q<.001,sprintf('%.5f',q),qfmt(q))),size=2.8)+annotate('text',x=qx,y=7.7,label='q (BH)',size=3,fontface='bold')+scale_y_continuous(breaks=g$y,labels=c('Gravity (G)','Inoculation (I)','Day (D)','G \u00d7 I','G \u00d7 D','I \u00d7 D','G \u00d7 I \u00d7 D'),limits=c(.5,8))+scale_x_continuous(limits=c(0,qx+gx*.17),breaks=c(0,5,10,15,20),expand=expansion(mult=0))+labs(x='Conditional partial R\u00b2 (%)',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=9))
 h<-rd('Figure_4','panel_d_four_context_effect_matrix_long');sel<-rd('Figure_4','panel_d_all_selected_rows_and_order');contexts<-c('inoculation_stationary_day3','inoculation_clinorotation_day3','inoculation_stationary_day7','inoculation_clinorotation_day7');h$x<-match(h$context,contexts);h$y<-580-h$row_order-5*(as.integer(sub('C','',h$cluster))-1);h$display_effect<-pmax(-4,pmin(4,h$log2_effect));mid<-aggregate(y~cluster,h,mean);cnt<-table(sel$cluster);mid$label<-paste0(mid$cluster,' (',cnt[mid$cluster],')')
 pc<-ggplot(h,aes(x,y,fill=display_effect))+geom_tile(width=.93,height=1)+scale_fill_gradient2(low='#3C6EA8',mid='#FAFAFA',high='#BF5139',limits=c(-4,4),name=expression('Inoculation effect'~(log[2])),breaks=c(-4,-2,0,2,4),labels=c('\u2264-4','-2','0','2','\u22654'))+scale_x_continuous(breaks=1:4,labels=c('D3\n1 g\n30 hits','D3\ns\u03bcg\n67 hits','D7\n1 g\n182 hits','D7\ns\u03bcg\n412 hits'),position='top',expand=expansion(add=.48))+scale_y_continuous(breaks=mid$y,labels=mid$label,expand=expansion(mult=0))+labs(x=NULL,y=NULL)+theme_void(base_family=font,base_size=10)+theme(axis.text.x=element_text(size=9,colour='black',margin=margin(b=5)),axis.text.y=element_text(size=10,colour='black',margin=margin(r=5)),legend.position='bottom',legend.title.position='top',legend.title=element_text(size=9),legend.text=element_text(size=9),plot.margin=margin(2,5,0,5))+guides(fill=guide_colourbar(direction='horizontal',barwidth=mm(75),barheight=mm(3),title.position='top'))
 write_source(h,'Figure_4','heatmap_plot_coordinates_and_display_values')
 export('Figure_4',190,230,function(){place(0,0,93,90,function()panel(pa,'a','Root-exudate profiles',93,90));place(97,0,93,90,function()panel(pb,'b','Factor effects',93,90));place(0,98,190,132,function()panel(pc,'c','Inoculation responses',190,132))})
}
draw_Figure_5<-function(){
 d<-rd('Figure_5','panel_a_interaction_time_scatter')
 stopifnot(nrow(d)==3126,sum(d$T_pass)==5,sum(d$J3_pass)==0,sum(d$J7_pass)==137)
 ids<-c('pos_5231','neg_3047')
 display_names<-c('Isorhamnetin glycoside','Sinapoyl glucose')
 d$direction<-ifelse(d$T>0,'Positive','Negative');lim<-ceiling(max(abs(c(d$J3,d$J7))))
 pa<-ggplot(d,aes(J3,J7))+
  geom_abline(slope=1,intercept=0,linetype=2,colour='#8598A5',linewidth=.4)+
  geom_hline(yintercept=0,colour='#DAE1E5',linewidth=.25)+geom_vline(xintercept=0,colour='#DAE1E5',linewidth=.25)+
  geom_point(data=d[!d$T_pass,],colour='#A9B4BC',alpha=.32,size=.75)+
  geom_point(data=d[d$T_pass,],aes(colour=direction),size=2.2)+
  scale_colour_manual(values=c(Positive=orange,Negative=blue),breaks=c('Positive','Negative'),labels=c('T > 0 (4)','T < 0 (1)'))+
  scale_x_continuous(limits=c(-lim,lim),breaks=c(-10,-5,0,5,10))+
  scale_y_continuous(limits=c(-lim,lim),breaks=c(-10,-5,0,5,10))+coord_equal()+
  labs(x=expression('D3 interaction'~J[3]~(log[2])),y=expression('D7 interaction'~J[7]~(log[2])),subtitle='Five responsive temporal interactions')+
  theme(plot.subtitle=element_text(hjust=.5,size=8.5),legend.text=element_text(size=8),legend.key.width=mm(4),plot.margin=margin(2,3,1,2))

 ann<-rd('Figure_5','feature_annotations_and_quality')
 screen<-d[,c('feature_id','J7','J7_q','J7_pass','T','T_q','T_pass')]
 ai<-match(screen$feature_id,ann$feature_id);stopifnot(!anyNA(ai),!anyDuplicated(screen$feature_id))
 for(k in c('putative_name','qc_cv','mass_error_ppm','any_group_half_at_or_below_one','strict_annotation_diagnostic','hmdb_superclass','kegg_id','hmdb_id'))screen[[k]]<-ann[[k]][ai]
 screen$minus_log10_J7_q<--log10(screen$J7_q)
 screen$screen_group<-ifelse(!screen$J7_pass,'Other',ifelse(screen$J7>0,'Positive','Negative'))
 screen$selected_example<-screen$feature_id%in%ids
 screen$annotation_subset<-screen$J7_pass & screen$strict_annotation_diagnostic & screen$hmdb_superclass=='Phenylpropanoids and polyketides'
 cutoff<-log2(1.5)
 stopifnot(identical(screen$J7_pass,screen$J7_q<.05 & abs(screen$J7)>cutoff))
 sel<-screen[match(ids,screen$feature_id),]
 sel$display_name<-display_names;sel$plot_label<-c('Isorhamnetin\nglycoside','Sinapoyl\nglucose')
 sel$label_x<-c(5.2,10);sel$label_y<-c(3.85,3.2)
 stopifnot(all(sel$J7_pass),all(sel$annotation_subset),!any(sel$T_pass))
 npos<-sum(screen$J7_pass & screen$J7>0);nneg<-sum(screen$J7_pass & screen$J7<0)
 xlim<-ceiling(max(abs(screen$J7)));ymax<-max(screen$minus_log10_J7_q)+.95
 pb<-ggplot(screen,aes(J7,minus_log10_J7_q))+
  geom_vline(xintercept=c(-cutoff,cutoff),linetype=2,colour='#8E9CA5',linewidth=.3)+
  geom_hline(yintercept=-log10(.05),linetype=2,colour='#8E9CA5',linewidth=.3)+
  geom_point(data=screen[!screen$J7_pass,],colour='#B2BCC3',alpha=.40,size=.8)+
  geom_point(data=screen[screen$J7_pass,],aes(colour=screen_group),size=1.4,alpha=.8)+
  geom_point(data=sel,shape=21,fill=orange,colour='black',size=2.65,stroke=.65)+
  geom_segment(data=sel,aes(x=J7,y=minus_log10_J7_q,xend=label_x,yend=label_y-.21),inherit.aes=FALSE,colour='#424A4F',linewidth=.3)+
  geom_text(data=sel,aes(x=label_x,y=label_y,label=plot_label),inherit.aes=FALSE,colour='black',size=2.8,family=font,lineheight=1.05)+
  scale_colour_manual(values=c(Positive=orange,Negative=blue),breaks=c('Positive','Negative'),
                      labels=c(sprintf('J7 > 0 (%d)',npos),sprintf('J7 < 0 (%d)',nneg)))+
  scale_x_continuous(limits=c(-xlim,xlim),breaks=c(-10,-5,0,5,10))+
  scale_y_continuous(limits=c(0,ymax),breaks=0:4,expand=expansion(mult=c(0,.01)))+
  labs(x=expression('Day-7 interaction'~J[7]~(log[2])),y=expression(-log[10](q[J[7]])),
       subtitle='137 responsive features; two examples labelled')+
  theme(plot.subtitle=element_text(hjust=.5,size=8.5),legend.text=element_text(size=8),legend.key.width=mm(4),plot.margin=margin(2,3,1,2))
 write_source(screen,'Figure_5','panel_b_day7_interaction_screen')
 write_source(sel,'Figure_5','panel_b_labelled_examples')
 write_source(screen[screen$annotation_subset,],'Figure_5','panel_b_annotation_subset_all_seven')

 cp<-rd('Figure_5','candidate_samples_all_days');cs<-rd('Figure_5','candidate_group_means_and_CI');ct<-rd('Figure_5','panel_d_candidate_all_contrasts')
 chem<-function(k){id<-ids[k];x<-cp[cp$day==7&cp$feature_id==id,];s<-cs[cs$day==7&cs$feature_id==id,];j<-ct[ct$feature_id==id&ct$contrast=='interaction_day7',];xp<-function(z)ifelse(z$gravity=='stationary_1g',0,2)+ifelse(z$bacteria=='uninoculated',1,2);x$x<-xp(x);s$x<-xp(s)
  ggplot(x,aes(x,log2_intensity))+
   geom_col(data=s,aes(x=x,y=mean_log2,fill=gravity),inherit.aes=FALSE,width=.55,colour=unname(pal_dark[s$gravity]),linewidth=.6)+
   geom_errorbar(data=s,aes(x=x,ymin=ci_low,ymax=ci_high),inherit.aes=FALSE,width=.13,colour='#6C7D86',linewidth=.45)+
   geom_point(aes(colour=gravity),position=position_jitter(width=.09,height=0,seed=20260912),size=1.9)+
   scale_colour_manual(values=pal,guide='none')+scale_fill_manual(values=pal_fill,guide='none')+
   scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.14)))+
   scale_x_continuous(breaks=1:4,labels=c('1 g\nU','1 g\nI','s\u03bcg\nU','s\u03bcg\nI'),limits=c(.5,4.5))+
   labs(x=NULL,y=expression(log[2](intensity+1)),title=display_names[k],
        subtitle=bquote(J[7]==.(sprintf('%.2f',j$log2_effect))*'; q = '*.(sprintf('%.4f',j$q))))+
   theme(plot.title=element_text(hjust=.5),plot.subtitle=element_text(hjust=.5,size=9))
 }
 write_source(cp[cp$day==7&cp$feature_id%in%ids,],'Figure_5','plotted_day7_candidate_observations')
 export('Figure_5',190,191,function(){
  place(0,0,93,103,function()panel(pa,'a','Interaction changes over time',93,103))
  place(97,0,93,103,function()panel(pb,'b','Day-7 interaction candidates',93,103))
  place(0,110,190,81,function(){header('c','Day-7 chemical responses');place(0,10,94,71,function()drawp(chem(1)));place(98,10,92,71,function()drawp(chem(2)))})
 })
 writeLines(c('a: temporal contrast T=J7-J3; five responsive features.',
              'b: separate day-7 J7 screen; all 3,126 tests retained; 137 responsive, 108 positive and 29 negative.',
              'Two existing annotated examples are highlighted within the day-7 screen, not identified as temporal hits.',
              'Seven day-7 hits meet the archived stricter annotation diagnostic and the phenylpropanoid/polyketide class; all seven are supplied.',
              'The two highlighted examples illustrate distinct inoculation response patterns; they are not a unique algorithmic top-two selection.',
              'c: original two group-profile plots, with unchanged observations, means and 95% intervals.'),file.path(pkg,'qa/Figure_5_selection_logic.txt'))
}

draw_Figure_S6<-function(){
 d<-rd('Figure_S6','panel_c_gravity_time_source_data');d$direction<-ifelse(d$time_pass,'Temporal change','Other');lim<-ceiling(max(abs(c(d$day3_log2_effect,d$day7_log2_effect))))
 p<-ggplot(d,aes(day3_log2_effect,day7_log2_effect))+geom_abline(slope=1,intercept=0,linetype=2,colour='#7B8F9B',linewidth=.4)+geom_hline(yintercept=0,colour='#DFE5E8',linewidth=.25)+geom_vline(xintercept=0,colour='#DFE5E8',linewidth=.25)+geom_point(data=d[!d$time_pass,],colour='#ADB8C0',alpha=.35,size=1)+geom_point(data=d[d$time_pass,],colour=blue,size=2.7)+coord_equal()+scale_x_continuous(limits=c(-lim,lim))+scale_y_continuous(limits=c(-lim,lim))+labs(x=expression('Day-3 gravity effect'~G[3]~(log[2])),y=expression('Day-7 gravity effect'~G[7]~(log[2])),subtitle='Day-specific responses: D3 = 0; D7 = 31\nBlue: three significant temporal changes in G')+theme(plot.subtitle=element_text(hjust=.5,size=9))
 export('Figure_S6',130,125,function()drawp(p))
}
