draw_Figure_S2<-function(){
 d<-rd('Figure_S2','RNA_PCA');p<-ggplot(d,aes(PC1,PC2,colour=gravity))+geom_point(size=3)+geom_text_repel(aes(label=sample),size=3.4,seed=20261104,min.segment.length=0,box.padding=.4,max.overlaps=Inf)+scale_colour_manual(values=pal,labels=glab,breaks=names(pal))+labs(x=sprintf('PC1 (%.2f%%)',d$PC1_percent[1]),y=sprintf('PC2 (%.2f%%)',d$PC2_percent[1]))+theme(legend.position='bottom')
 export('Figure_S2',135,110,function()drawp(p))
}
gene_heatmap<-function(fig,row_height=3.6){
 z<-rd(fig,'gene_statistics');m<-rd(fig,'gene_row_z');M<-as.matrix(m[,-1]);rownames(M)<-m$gene_id;M<-M[z$gene_id,];n<-nrow(z);top<-14;cw<-7;left<-62;yp<-function(i)top+(i-.5)*row_height
 zcol<-colorRamp2(c(-2,0,2),c('#3C6EA8','#FAFAFA','#BF5139'));fcol<-colorRamp2(c(-3,0,3),c('#3C6EA8','#FAFAFA','#BF5139'))
 txt('Gene',58,0,9,TRUE,just='right');txt(expression(log[2]*FC),137,0,9,TRUE,just='centre');txt('BH q',157,0,9,TRUE,just='centre')
 for(j in 1:9){x<-left+(j-.5)*cw;col<-if(j<=5)blue else orange;txt(colnames(M)[j],x,7,9,col=col,just='centre')}
 for(i in 1:n){
  yy<-yp(i);txt(z$display_symbol[i],58,yy-1.9,9.3,just='right')
  for(j in 1:9)grid.rect(x=mm(left+(j-.5)*cw),y=unit(1,'npc')-mm(yy),width=mm(cw),height=mm(row_height),gp=gpar(fill=zcol(max(-2,min(2,M[i,j]))),col=NA))
  grid.rect(x=mm(137),y=unit(1,'npc')-mm(yy),width=mm(8),height=mm(row_height),gp=gpar(fill=fcol(max(-3,min(3,z$logFC[i]))),col=NA));txt(qfmt(z$adj.P.Val[i]),157,yy-1.9,9.3,just='centre')
 }
 if(n>8){
  mods<-unique(z$module)
  for(mod in mods){ii<-which(z$module==mod);lab<-paste(strwrap(gsub('[\r\n]+',' ',mod),width=17),collapse='\n');txt(lab,1,mean(yp(ii))-4.2,9,TRUE)}
 }
 end<-top+n*row_height;txt('Row z score',left,end+4,9)
 for(k in 1:60)grid.rect(x=mm(left+(k-.5)*.9),y=unit(1,'npc')-mm(end+11),width=mm(.9),height=mm(2.8),gp=gpar(fill=zcol(-2+4*(k-1)/59),col=NA))
 txt('\u2264-2',left,end+14,8,just='centre');txt('0',left+27,end+14,8,just='centre');txt('\u22652',left+54,end+14,8,just='centre')
 txt(expression(log[2]*FC),130,end+4,9)
 for(k in 1:40)grid.rect(x=mm(130+(k-.5)*.7),y=unit(1,'npc')-mm(end+11),width=mm(.7),height=mm(2.8),gp=gpar(fill=fcol(-3+6*(k-1)/39),col=NA))
 txt('\u2264-3',130,end+14,8,just='centre');txt('0',144,end+14,8,just='centre');txt('\u22653',158,end+14,8,just='centre')
 write_source(data.frame(gene_id=z$gene_id,display_symbol=z$display_symbol,plot_order=seq_len(n),raw_log2FC=z$logFC,display_log2FC=pmax(-3,pmin(3,z$logFC)),q=z$adj.P.Val),fig,'heatmap_row_statistics_and_display_values')
 Md<-M;Md[]<-pmax(-2,pmin(2,Md));write_source(data.frame(gene_id=rownames(Md),Md,check.names=FALSE),fig,'heatmap_display_z_scores')
}
draw_Figure_S3<-function(){export('Figure_S3',180,80,function()gene_heatmap('Figure_S3',5.3))}
draw_Figure_S8<-function(){
 d<-rd('Figure_S8','selected_GO_statistics');d$y<-8:1;lab<-c('Plastid organization','Photosynthesis','Flavonoid biosynthesis','Cellular response to hypoxia','Detoxification','Phenol-containing\ncompound metabolism','Response to desiccation','Response to salt')
 p<-ggplot(d,aes(y=y))+geom_vline(xintercept=0,colour='#A2AFB7',linewidth=.35)+geom_segment(aes(x=0,xend=fgsea_NES,yend=y),colour='#B0BEC6',linewidth=.65)+geom_point(aes(x=fgsea_NES,colour=Direction,size=NGenes))+geom_text(aes(x=4.45,label=qfmt(FDR)),size=2.9)+geom_text(aes(x=6.65,label=qfmt(FRY_q)),size=2.9)+annotate('text',x=4.45,y=9,label='CAMERA q',size=2.9,fontface='bold')+annotate('text',x=6.65,y=9,label='FRY q',size=2.9,fontface='bold')+scale_x_continuous(limits=c(-3.3,7.7),breaks=c(-3,-2,-1,0,1,2,3),expand=expansion(mult=0))+scale_y_continuous(breaks=d$y,labels=lab,limits=c(.4,9.5),expand=expansion(mult=0))+scale_colour_manual(values=c(Up=orange,Down=blue),labels=c(Up='Positive',Down='Negative'))+scale_size_continuous(range=c(2.5,5),breaks=c(50,150,300),name='Genes')+labs(x='fgsea normalized enrichment score',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=9),legend.text=element_text(size=8),legend.title=element_text(size=8),legend.key.width=mm(5))
 export('Figure_S8',190,264,function(){place(0,0,190,98,function()panel(p,'a','Selected functional themes',190,98));place(0,105,190,159,function(){header('b','Representative genes from four functional themes');place(0,10,190,149,function()gene_heatmap('Figure_S8',3.6))})})
}
