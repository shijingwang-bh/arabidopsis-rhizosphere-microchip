options(stringsAsFactors=FALSE,width=150)
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
args<-commandArgs(TRUE);pkg<-normalizePath(args[1],winslash='/');run<-normalizePath(args[2],winslash='/')
suppressPackageStartupMessages({library(edgeR);library(limma);library(fgsea);library(BiocParallel)})
rd<-function(n,...)read.csv(file.path(pkg,'data/rna',n),check.names=FALSE,...)
map<-c(Figure_S2='rna_pca',Figure_S3='rna_degs',Figure_S8='rna_gene_sets')
wr<-function(x,n){
 if(startsWith(n,'shared/Salmon_analysis/'))p<-file.path(run,'results/rna',basename(n))
 else {parts<-strsplit(n,'/',fixed=TRUE)[[1]];p<-file.path(run,'plot_data',map[[parts[1]]],parts[2])}
 dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);write.csv(x,p,row.names=FALSE,na='')
}
ref<-'reference'
c0<-read.delim(file.path(pkg,'data/rna/cdna.count'),check.names=FALSE)
t0<-read.delim(file.path(pkg,'data/rna/cdna.tpm'),check.names=FALSE)
gene<-sub('\\.[0-9]+$','',c0$Name)
C<-rowsum(as.matrix(c0[,-1]),gene,reorder=FALSE);TPM<-rowsum(as.matrix(t0[,-1]),gene,reorder=FALSE)
meta<-rd('RNA_sample_metadata_pool_units.csv');meta<-meta[match(colnames(C),meta$sample),]
meta$gravity<-factor(meta$gravity,levels=c('stationary_1g','clinorotation'))
design<-model.matrix(~0+gravity,meta);colnames(design)<-levels(meta$gravity);contrast<-c(-1,1)
y<-DGEList(C,group=meta$gravity);keep<-filterByExpr(y,group=meta$gravity,min.count=5)
y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);v<-voom(y,design,plot=FALSE)
fit<-eBayes(contrasts.fit(lmFit(v,design),contrast),robust=TRUE)
tx<-topTable(fit,number=Inf,sort.by='none');tx$gene_id<-rownames(tx)
old<-rd(paste0(ref,'/gene_level_limma_voom.csv'));j<-match(tx$gene_id,old$gene_id)
stopifnot(nrow(tx)==20449,!anyNA(j),max(abs(tx$logFC-old$logFC[j]))<1e-8,max(abs(tx$adj.P.Val-old$adj.P.Val[j]))<1e-8,sum(tx$adj.P.Val<.05)==8)
tx$symbol<-old$symbol[j];tx$description<-old$description[j]
wr(tx,'shared/Salmon_analysis/all_gene_statistics.csv')
wr(data.frame(gene_id=rownames(C),C,check.names=FALSE),'shared/Salmon_analysis/gene_counts.csv')
wr(data.frame(gene_id=rownames(TPM),TPM,check.names=FALSE),'shared/Salmon_analysis/gene_TPM.csv')
wr(data.frame(gene_id=rownames(v$E),v$E,check.names=FALSE),'shared/Salmon_analysis/tested_gene_logCPM.csv')
pc<-prcomp(t(v$E),center=TRUE,scale.=FALSE);pv<-100*pc$sdev^2/sum(pc$sdev^2)
wr(data.frame(sample=meta$sample,gravity=meta$gravity,pc$x[,1:4],PC1_percent=pv[1],PC2_percent=pv[2]),'Figure_S2/RNA_PCA.csv')
wr(meta,'Figure_S2/sample_metadata.csv')
vq<-voomWithQualityWeights(y,design,plot=FALSE)
wr(data.frame(sample=meta$sample,gravity=meta$gravity,quality_weight=vq$targets$sample.weights,TMM_factor=y$samples$norm.factors),'shared/Salmon_analysis/sample_QC.csv')

# Frozen memberships avoid annotation database drift.
gm<-rd('GO_membership.csv');terms<-rd('GO_dictionary.csv')
paths<-split(gm$gene_id,gm$GO_ID)
paths<-paths[order(names(paths))];stopifnot(length(paths)==1705)
ix<-ids2indices(paths,rownames(v$E));ca<-camera(v,ix,design,contrast=contrast,inter.gene.cor=.01,sort=FALSE)
ce<-camera(v,ix,design,contrast=contrast,inter.gene.cor=NA,sort=FALSE);fr<-fry(v,ix,design,contrast=contrast,sort='none')
set.seed(20260911);stats<-sort(setNames(tx$t,tx$gene_id),decreasing=TRUE)
fg<-as.data.frame(fgseaMultilevel(pathways=paths,stats=stats,minSize=15,maxSize=500,eps=0,nPermSimple=10000,sampleSize=101,BPPARAM=SerialParam()))
fg$leadingEdge<-vapply(fg$leadingEdge,paste,collapse=';',character(1));wr(fg,'shared/Salmon_analysis/fgsea_all_GO.csv')
ca$GO_ID<-rownames(ca);ca$term<-terms$TERM[match(ca$GO_ID,terms$GOID)]
ca$estimated_rho_q<-ce$FDR[match(ca$GO_ID,rownames(ce))];ca$FRY_q<-fr$FDR[match(ca$GO_ID,rownames(fr))]
ca$fgsea_NES<-fg$NES[match(ca$GO_ID,fg$pathway)];ca$fgsea_q<-fg$padj[match(ca$GO_ID,fg$pathway)]
wr(ca,'shared/Salmon_analysis/all_GO_statistics.csv')
sel<-rd(paste0(ref,'/Figure2_functional_panel_source.csv'))
wr(ca[match(sel$GO_ID,ca$GO_ID),],'Figure_S8/selected_GO_statistics.csv')
for(spec in list(c('Figure_S3','Supplement_eight_DEG_heatmap_source.csv'),c('Figure_S8','Figure2c_selected_32_genes_and_statistics.csv'))){
 q<-rd(paste0(ref,'/',spec[2]));id<-q$gene_id
 z<-tx[match(id,tx$gene_id),];z$display_symbol<-q$display_symbol
 z$display_symbol[z$gene_id=='AT2G11810']<-'MGD3';z$display_symbol[z$gene_id=='AT4G23430']<-'TIC32';z$display_symbol[z$gene_id=='AT2G25625']<-'CV'
 z$module<-q$module;z$plot_order<-seq_along(id)
 wr(z,paste0(spec[1],'/gene_statistics.csv'))
 for(n in c('counts','TPM','logCPM','row_z')){
  M<-switch(n,counts=C[id,],TPM=TPM[id,],logCPM=v$E[id,],row_z=t(scale(t(v$E[id,]))))
  wr(data.frame(gene_id=id,M,check.names=FALSE),paste0(spec[1],'/gene_',n,'.csv'))
 }
}

stopifnot(sum(ca$estimated_rho_q<.05)==0L)
wr(data.frame(GO_ID=rownames(ce),ce,check.names=FALSE),'shared/Salmon_analysis/CAMERA_estimated_correlation.csv')
wr(data.frame(GO_ID=rownames(fr),fr,check.names=FALSE),'shared/Salmon_analysis/FRY_directional.csv')
capture.output(sessionInfo(),file=file.path(run,'results/rna/sessionInfo.txt'))
cat('PASS: 20,449 genes, eight DEGs, 1,705 GO sets, frozen annotations.\n')
