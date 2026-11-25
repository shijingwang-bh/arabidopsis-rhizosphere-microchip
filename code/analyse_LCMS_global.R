options(stringsAsFactors=FALSE,scipen=10)
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
args<-commandArgs(TRUE);pkg<-normalizePath(args[1],winslash='/');run<-normalizePath(args[2],winslash='/')
src<-file.path(pkg,'data/lcms');out<-file.path(run,'results/chemistry_global')
for(d in c('tables','data','qa'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(vegan);library(jsonlite)})
write_tab <- function(x, name) write.csv(x, file.path(out, "tables", paste0(name,".csv")), row.names=FALSE, na="", fileEncoding="UTF-8")
read_export <- function(n) {
  txt <- readLines(file.path(src,n), warn=FALSE, encoding="bytes")
  converted <- iconv(txt, from="GB18030", to="UTF-8", sub="byte")
  stopifnot(!anyNA(converted))
  bad <- grep("<[0-9a-f]{2}>", converted)
  if(length(bad)) write.csv(data.frame(line=bad,feature_id=sub("\t.*$","",converted[bad])),
    file.path(out,"qa",paste0(n,"_encoding_substitutions.csv")),row.names=FALSE)
  read.delim(text=paste(converted,collapse="\n"), check.names=FALSE, quote="", comment.char="", encoding="UTF-8")
}
message("Read and audit inputs")
named <- read_export("metabolites_exp.xls")
full <- read_export("All_metabolites_ori.xls")
ann <- read_export("metabolites_full_table.xls")
names(named)[1] <- names(full)[1] <- names(ann)[1] <- "feature_id"
meta_all <- read.csv(file.path(src,"corrected_sample_metadata_with_QC.csv"))
stopifnot(all(meta_all$bacteria[grepl("[37]N[0-9]+$",meta_all$sample)]=="uninoculated"), all(meta_all$bacteria[grepl("[37]P[0-9]+$",meta_all$sample)]=="inoculated"))
meta <- meta_all[meta_all$gravity != "QC",]; meta$day <- as.integer(meta$day)
qc_ids <- meta_all$sample[meta_all$gravity == "QC"]
stopifnot(nrow(meta)==41, length(qc_ids)==7, !anyDuplicated(named$feature_id), !anyDuplicated(full$feature_id), !anyDuplicated(meta$sample))
raw_all <- as.matrix(full[,meta$sample]); storage.mode(raw_all) <- "double"; rownames(raw_all) <- full$feature_id
raw_qc <- as.matrix(full[,qc_ids]); storage.mode(raw_qc) <- "double"; rownames(raw_qc) <- full$feature_id
stopifnot(all(is.finite(raw_all)),all(is.finite(raw_qc)),all(raw_all>=0),all(raw_qc>=0))
raw_named <- as.matrix(named[,meta$sample]); storage.mode(raw_named) <- "double"
delta_input <- max(abs(raw_named - raw_all[match(named$feature_id,full$feature_id),]))
stopifnot(delta_input < 1e-6)
# Preserve the exact named-export numbers used in the existing analysis;
# the full export differs only at floating-point text precision for these rows.
raw_all[match(named$feature_id,full$feature_id),] <- raw_named
raw_qc[match(named$feature_id,full$feature_id),] <- as.matrix(named[,qc_ids])
qc_cv <- apply(raw_qc,1,sd) / rowMeans(raw_qc)
groups <- sort(unique(meta$group))
low_group <- sapply(groups,function(g) rowMeans(raw_all[,meta$group==g,drop=FALSE]<=1)>=0.5)
feature <- data.frame(feature_id=full$feature_id, putative_name=full$name,
  named_family=full$feature_id %in% named$feature_id, qc_cv=qc_cv,
  qc_pass=is.finite(qc_cv) & qc_cv<=0.30,
  any_group_half_at_or_below_one=apply(low_group,1,any))
feature$primary <- feature$named_family & feature$qc_pass
annotation_map <- c(mass_error_ppm="Mass Error (ppm)",formula="Formula",mz="m/z",retention_min="Retention time (min)",
  kegg_id="KEGG_annotation",hmdb_id="HMDB_ID",hmdb_superclass="super_class(HMDB)",hmdb_class="HMDB_taxonomy",hmdb_subclass="sub_class(HMDB)")
for (nm in names(annotation_map)) feature[[nm]] <- ann[[annotation_map[nm]]][match(feature$feature_id,ann$feature_id)]
valid_id <- function(x) !is.na(x) & !x %in% c("--","-","","NA")
feature$strict_annotation_diagnostic <- feature$named_family & !is.na(feature$mass_error_ppm) & abs(feature$mass_error_ppm)<=5 & (valid_id(feature$kegg_id)|valid_id(feature$hmdb_id))
stopifnot(sum(feature$primary)==3126, sum(feature$qc_pass)==14031, sum(feature$named_family)==3414)
write_tab(feature,"feature_quality_and_family_membership")
write_tab(meta_all,"sample_metadata")
write_tab(as.data.frame(table(meta$group),responseName="n"),"biological_group_sizes")
L <- log2(raw_all+1)
coef_simple <- function(day, gravity=NULL, bacteria=NULL) {
  g <- meta[match(groups,meta$group),]
  v <- rep(0,length(groups)); names(v) <- groups
  eligible <- g$day==day
  if (!is.null(gravity)) {
    eligible <- eligible & g$gravity==gravity
    v[eligible] <- ifelse(g$bacteria[eligible]=="inoculated",1,-1)
  } else {
    eligible <- eligible & g$bacteria==bacteria
    v[eligible] <- ifelse(g$gravity[eligible]=="clinorotation",1,-1)
  }
  v
}
C <- list(
 gravity_uninoculated_day3=coef_simple(3,bacteria="uninoculated"),
 gravity_uninoculated_day7=coef_simple(7,bacteria="uninoculated"),
 inoculation_stationary_day3=coef_simple(3,gravity="stationary_1g"),
 inoculation_clinorotation_day3=coef_simple(3,gravity="clinorotation"),
 inoculation_stationary_day7=coef_simple(7,gravity="stationary_1g"),
 inoculation_clinorotation_day7=coef_simple(7,gravity="clinorotation"))
C$interaction_day3 <- C$inoculation_clinorotation_day3-C$inoculation_stationary_day3
C$interaction_day7 <- C$inoculation_clinorotation_day7-C$inoculation_stationary_day7
C$gravity_time_uninoculated <- C$gravity_uninoculated_day7-C$gravity_uninoculated_day3
C$inoculation_time_stationary <- C$inoculation_stationary_day7-C$inoculation_stationary_day3
C$inoculation_time_clinorotation <- C$inoculation_clinorotation_day7-C$inoculation_clinorotation_day3
C$gravity_inoculation_time <- C$interaction_day7-C$interaction_day3
coef_mat <- do.call(rbind,C)
write_tab(data.frame(contrast=rownames(coef_mat),coef_mat,check.names=FALSE),"contrast_definitions")
threshold <- log2(1.5)
contrast_stats <- function(z, md, family) {
  means <- sapply(groups,function(g) rowMeans(z[,md$group==g,drop=FALSE]))
  vars <- sapply(groups,function(g) apply(z[,md$group==g,drop=FALSE],1,var))
  ns <- sapply(groups,function(g) sum(md$group==g))
  result <- lapply(names(C),function(ct) {
    v <- C[[ct]]; est <- as.vector(means %*% v)
    terms <- sweep(vars,2,v^2/ns,"*"); se2 <- rowSums(terms)
    df <- se2^2 / rowSums(sweep(terms^2,2,ns-1,"/"))
    se <- sqrt(se2); t <- est/se
    p <- 2*pt(-abs(t),df); p[!is.finite(p)] <- 1
    q <- p.adjust(p,"BH"); ci <- qt(0.975,df)*se
    data.frame(family=family,contrast=ct,feature_id=rownames(z),log2_effect=est,se=se,df=df,t=t,p=p,q=q,
      ci95_low=est-ci,ci95_high=est+ci,pass=q<0.05 & abs(est)>threshold)
  })
  do.call(rbind,result)
}
summarize_contrasts <- function(x) do.call(rbind,lapply(split(x,x$contrast),function(d)
  data.frame(family=d$family[1],contrast=d$contrast[1],tested=nrow(d),q05=sum(d$q<0.05),
  displayed=sum(d$pass),positive=sum(d$pass & d$log2_effect>0),negative=sum(d$pass & d$log2_effect<0))))
families <- list(primary_named_QC30=which(feature$primary), original_named=which(feature$named_family), all_export_QC30=which(feature$qc_pass))
stats_list <- list(); global_list <- list()
for (fam in names(families)) {
  message("Feature contrasts: ",fam)
  stats_list[[fam]] <- contrast_stats(L[families[[fam]],,drop=FALSE],meta,fam)
  write_tab(stats_list[[fam]],paste0("feature_contrasts_",fam))
}
primary <- stats_list[[1]]
summary_contrasts <- do.call(rbind,lapply(stats_list,summarize_contrasts))
write_tab(summary_contrasts,"contrast_counts_family_sensitivity")
# Compare with the CURRENT corrected complete feature-test family.
old<-read.csv(file.path(pkg,'reference/plot_data/chemistry_global/feature_contrasts_primary_named_QC30.csv'))
key<-match(paste(old$contrast,old$feature_id),paste(primary$contrast,primary$feature_id))
stopifnot(!anyNA(key))
reproduction<-data.frame(metric=c('log2_effect','p','q'),maximum_absolute_difference=sapply(c('log2_effect','p','q'),function(n)max(abs(primary[[n]][key]-old[[n]]))))
stopifnot(all(reproduction$maximum_absolute_difference<1e-8))
write_tab(reproduction,'current_primary_family_reproduction')
# All seven equal-cell factorial contrasts are tested conditional on the other six.
design <- function(md) {
  G <- ifelse(md$gravity=="clinorotation",1,-1)
  I <- ifelse(md$bacteria=="inoculated",1,-1)
  D <- ifelse(md$day==7,1,-1)
  cbind(Intercept=1,G=G,I=I,D=D,GI=G*I,GD=G*D,ID=I*D,GID=G*I*D)
}
pareto <- function(z) {
  center <- rowMeans(z); s <- apply(z,1,sd); s[!is.finite(s)|s==0] <- 1
  list(y=t(sweep(sweep(z,1,center,"-"),1,sqrt(s),"/")), center=center, scale=sqrt(s))
}
projector <- function(x) {q <- qr.Q(qr(x)); tcrossprod(q)}
set.seed(20260912)
perms <- t(replicate(9999,sample.int(nrow(meta))))
saveRDS(perms,file.path(out,"data/permutation_indices_41.rds"))
global_model <- function(z,md,family,permmat,validate=FALSE) {
  yy <- pareto(z); y <- yy$y; X <- design(md); n <- nrow(X)
  stopifnot(qr(X)$rank==8)
  H <- projector(X); K <- tcrossprod(y); rss <- sum(diag(K))-sum(H*K); tss <- sum(diag(K)); rdf <- n-8
  # Eigen-coordinates retain the full sample Gram matrix, reducing computation only.
  ee <- eigen(K,symmetric=TRUE); keep <- ee$values>max(ee$values)*1e-12
  coordinates <- sweep(ee$vectors[,keep,drop=FALSE],2,sqrt(ee$values[keep]),"*")
  stopifnot(max(abs(tcrossprod(coordinates)-K))/max(abs(K))<1e-10)
  validations <- list()
  tab <- do.call(rbind,lapply(2:8,function(j) {
    fit <- rda(coordinates,X[,j,drop=FALSE],X[,-c(1,j),drop=FALSE])
    Hr <- projector(X[,-j,drop=FALSE]); Rr <- diag(n)-Hr
    qv <- as.vector(Rr %*% X[,j]); qv <- qv/sqrt(sum(qv^2)); Ht <- tcrossprod(qv)
    ss <- sum(Ht*K); ff <- ss/(rss/rdf)
    vegan_observed_F <- fit$CCA$tot.chi/(fit$CA$tot.chi/rdf)
    stopifnot(abs(vegan_observed_F-ff)<1e-8)
    # Freedman-Lane: Y* = H_reduced Y + P R_reduced Y, then refit the
    # original full and reduced models. The fitted component is annihilated
    # by both the conditional-term and full-residual projections.
    K0 <- Rr %*% K %*% Rr
    permuted_F <- vapply(seq_len(nrow(permmat)),function(b) {
      kp <- K0[permmat[b,],permmat[b,]]
      sum(Ht*kp)/((sum(diag(K0))-sum(H*kp))/rdf)
    },numeric(1))
    p_permutation <- (1+sum(permuted_F>=ff-1e-12))/(nrow(permmat)+1)
    if(validate) {
      e0 <- Rr %*% coordinates; fitted0 <- Hr %*% coordinates
      direct_refit_F <- vapply(seq_len(199),function(b) {
        ystar <- fitted0+e0[permmat[b,],,drop=FALSE]
        refit <- rda(ystar,X[,j,drop=FALSE],X[,-c(1,j),drop=FALSE])
        refit$CCA$tot.chi/(refit$CA$tot.chi/rdf)
      },numeric(1))
      delta <- max(abs(direct_refit_F-permuted_F[1:199]))
      stopifnot(delta<1e-7)
      validations[[colnames(X)[j]]] <<- data.frame(term=colnames(X)[j],F_difference=abs(vegan_observed_F-ff),permuted_F_max_difference=delta,explicit_vegan_refits=199)
    }
    data.frame(family=family,term=colnames(X)[j],n=n,features=nrow(z),df=1,residual_df=rdf,
      conditional_SS=ss,residual_SS=rss,total_SS=tss,partial_R2=ss/(ss+rss),conditional_SS_over_total=ss/tss,
      F=ff,p=p_permutation,full_model_R2=1-rss/tss)
  }))
  tab$q <- p.adjust(tab$p,"BH")
  if(validate) write_tab(do.call(rbind,validations),"global_model_projection_validation")
  tab
}
for (fam in names(families)) {
  message("Full-factorial global tests: ",fam)
  global_list[[fam]] <- global_model(L[families[[fam]],,drop=FALSE],meta,fam,perms,validate=fam==names(families)[1])
}
global <- do.call(rbind,global_list)
write_tab(global,"global_factorial_effects_family_sensitivity")
write_tab(data.frame(sample=meta$sample,design(meta),check.names=FALSE),"global_factorial_design_matrix")

message("PCA and multivariate dispersion")
primary_ids <- feature$feature_id[feature$primary]
scaled <- pareto(L[primary_ids,,drop=FALSE])
pca <- prcomp(scaled$y,center=FALSE,scale.=FALSE)
pca_fraction <- pca$sdev^2/sum(pca$sdev^2)
pca_scores <- cbind(meta,pca$x[,1:6])
qcz <- t(sweep(sweep(log2(raw_qc[primary_ids,,drop=FALSE]+1),1,scaled$center,"-"),1,scaled$scale,"/"))
qc_scores <- cbind(meta_all[meta_all$gravity=="QC",],qcz %*% pca$rotation[,1:6])
write_tab(pca_scores,"PCA_biological_scores")
write_tab(qc_scores,"PCA_QC_projected_scores")
write_tab(data.frame(PC=paste0("PC",seq_along(pca_fraction)),variance_fraction=pca_fraction),"PCA_variance_explained")
write_tab(data.frame(feature_id=primary_ids,center=scaled$center,pareto_divisor=scaled$scale,pca$rotation[,1:6]),"PCA_feature_loadings_and_scaling")
dispersion <- list(); disp_points <- list()
for (nm in c("group","gravity","bacteria","day")) {
  f <- factor(meta[[nm]])
  fit <- betadisper(dist(scaled$y),f,type="median",bias.adjust=TRUE)
  test <- permutest(fit,permutations=perms)
  dispersion[[nm]] <- data.frame(grouping=nm,F=test$tab$F[1],p=test$tab$`Pr(>F)`[1],groups=nlevels(f),bias_adjusted=TRUE)
  disp_points[[nm]] <- data.frame(sample=meta$sample,grouping=nm,group=as.character(f),distance_to_spatial_median=fit$distances)
}
dispersion <- do.call(rbind,dispersion); dispersion$q_four_diagnostics <- p.adjust(dispersion$p,"BH")
write_tab(dispersion,"multivariate_dispersion_diagnostics")
write_tab(do.call(rbind,disp_points),"multivariate_dispersion_source_data")

message("Six neutral leave-one-out analyses for the extra-entry cell")
loo_stats <- list(); loo_global <- list()
omissions <- meta$sample[meta$group=="clinorotation_inoculated_day7"]
set.seed(20260913); perms40 <- t(replicate(9999,sample.int(40)))
saveRDS(perms40,file.path(out,"data/permutation_indices_40.rds"))
for (omit in omissions) {
  use <- meta$sample!=omit; md <- meta[use,]; zz <- L[primary_ids,use,drop=FALSE]
  st <- contrast_stats(zz,md,paste0("omit_",omit))
  sm <- summarize_contrasts(st)
  sm$effect_correlation_with_primary <- vapply(sm$contrast,function(ct) {
    a <- st[st$contrast==ct,]; b <- primary[primary$contrast==ct,]; cor(a$log2_effect,b$log2_effect[match(a$feature_id,b$feature_id)])
  },numeric(1))
  sm$primary_hits_retained <- vapply(sm$contrast,function(ct) {
    a <- st[st$contrast==ct,]; b <- primary[primary$contrast==ct,]; sum(a$pass & b$pass[match(a$feature_id,b$feature_id)])
  },numeric(1))
  loo_stats[[omit]] <- sm
  loo_global[[omit]] <- global_model(zz,md,paste0("omit_",omit),perms40)
  # Save every feature-level result, so sensitivity is auditable beyond counts.
  write_tab(st,paste0("leave_one_out_",omit))
}
write_tab(do.call(rbind,loo_stats),"leave_one_out_contrast_summary")
write_tab(do.call(rbind,loo_global),"leave_one_out_global_effects")

message("Build panel source data and transparent heatmap selection")
get_ct <- function(ct) primary[primary$contrast==ct,]
a <- get_ct("gravity_uninoculated_day3"); b <- get_ct("gravity_uninoculated_day7"); t <- get_ct("gravity_time_uninoculated")
stopifnot(identical(a$feature_id,b$feature_id),identical(a$feature_id,t$feature_id))
scatter <- data.frame(feature_id=a$feature_id,day3_log2_effect=a$log2_effect,day3_q=a$q,day3_pass=a$pass,
  day7_log2_effect=b$log2_effect,day7_q=b$q,day7_pass=b$pass,delta_day7_minus_day3=t$log2_effect,
  time_p=t$p,time_q=t$q,time_pass=t$pass,time_ci95_low=t$ci95_low,time_ci95_high=t$ci95_high)
scatter <- cbind(scatter,feature[match(scatter$feature_id,feature$feature_id),c("putative_name","qc_cv","any_group_half_at_or_below_one","mass_error_ppm","kegg_id","hmdb_id")])
write_tab(scatter,"panel_c_gravity_time_source_data")
contexts <- c("inoculation_stationary_day3","inoculation_clinorotation_day3","inoculation_stationary_day7","inoculation_clinorotation_day7")
effect <- sapply(contexts,function(ct) get_ct(ct)$log2_effect); rownames(effect) <- primary_ids
qs <- sapply(contexts,function(ct) get_ct(ct)$q); rownames(qs) <- primary_ids
passes <- qs<0.05 & abs(effect)>threshold
included <- rowSums(passes)>0
heat <- effect[included,,drop=FALSE]; heatq <- qs[included,,drop=FALSE]
hc <- hclust(dist(heat),method="ward.D2")
cluster_raw <- cutree(hc,k=6)
ord <- hc$order
cluster_order <- unique(cluster_raw[ord]); cluster <- match(cluster_raw,cluster_order)
names(cluster) <- rownames(heat)
row_ids <- rownames(heat)[ord]
representatives <- unlist(lapply(1:6,function(cl) {
  ids <- names(cluster)[cluster==cl]; m <- heat[ids,,drop=FALSE]
  d <- rowSums(sweep(m,2,colMeans(m),"-")^2)
  ids[order(d,ids)][seq_len(min(4,length(ids)))]
}))
centroid_distance <- setNames(numeric(nrow(heat)),rownames(heat))
for(cl in 1:6) {ids<-names(cluster)[cluster==cl]; centroid_distance[ids]<-sqrt(rowSums(sweep(heat[ids,,drop=FALSE],2,colMeans(heat[ids,,drop=FALSE]),"-")^2))}
selection <- data.frame(row_order=seq_along(row_ids),feature_id=row_ids,cluster=paste0("C",cluster[row_ids]),
  representative=row_ids %in% representatives,distance_to_cluster_centroid=centroid_distance[row_ids],heat[row_ids,,drop=FALSE])
selection <- cbind(selection,feature[match(row_ids,feature$feature_id),c("putative_name","qc_cv","any_group_half_at_or_below_one","mass_error_ppm","hmdb_superclass")])
write_tab(selection,"panel_d_all_selected_rows_and_order")
long <- do.call(rbind,lapply(seq_along(contexts),function(j) data.frame(feature_id=row_ids,row_order=seq_along(row_ids),
  cluster=paste0("C",cluster[row_ids]),context=contexts[j],log2_effect=heat[row_ids,j],q=heatq[row_ids,j],
  context_pass=heatq[row_ids,j]<0.05 & abs(heat[row_ids,j])>threshold)))
write_tab(long,"panel_d_four_context_effect_matrix_long")
write_tab(selection[match(representatives,selection$feature_id),],"supplement_representative_features")
cluster_summary <- do.call(rbind,lapply(1:6,function(cl) {
  ids<-names(cluster)[cluster==cl]
  data.frame(cluster=paste0("C",cl),features=length(ids),t(colMeans(heat[ids,,drop=FALSE])),check.names=FALSE)
}))
write_tab(cluster_summary,"panel_d_cluster_descriptive_means")
diagnostics <- list(named_raw_match_max_abs=delta_input,named_features=3414,primary_features=length(primary_ids),all_features=nrow(feature),
 all_QC30_features=sum(feature$qc_pass),biological_samples=nrow(meta),QC_injections=length(qc_ids),
 heatmap_union=nrow(heat),heatmap_color_max_abs=max(abs(heat)),scatter_max_abs=max(abs(c(scatter$day3_log2_effect,scatter$day7_log2_effect))),
 primary_floor_flagged=sum(feature$primary & feature$any_group_half_at_or_below_one),
 selected_floor_flagged=sum(selection$any_group_half_at_or_below_one),
 primary_three_way_FDR05=sum(get_ct("gravity_inoculation_time")$q<0.05),
 PCA_PC1=pca_fraction[1],PCA_PC2=pca_fraction[2],permutations=9999,seed=20260912)
write_json(diagnostics,file.path(out,"qa/numerical_audit.json"),auto_unbox=TRUE,pretty=TRUE,digits=14)
capture.output(sessionInfo(),file=file.path(out,"qa/sessionInfo_analysis.txt"))
saveRDS(list(meta=meta,feature=feature,primary=primary,summary=summary_contrasts,global=global,
 pca_scores=pca_scores,qc_scores=qc_scores,pca_fraction=pca_fraction,dispersion=dispersion,
 dispersion_points=do.call(rbind,disp_points),scatter=scatter,heat=heat,heatq=heatq,cluster=cluster,
 order=row_ids,representatives=representatives,selection=selection,cluster_summary=cluster_summary,
 loo_summary=do.call(rbind,loo_stats),loo_global=do.call(rbind,loo_global),contexts=contexts,diagnostics=diagnostics),file.path(out,"data/figure4_analysis.rds"))
message("COMPLETE")
print(diagnostics)
print(summary_contrasts[summary_contrasts$family=="primary_named_QC30",])
print(global[global$family=="primary_named_QC30",c("term","partial_R2","F","p","q","full_model_R2")])
print(dispersion)
