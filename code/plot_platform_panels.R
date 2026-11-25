args<-commandArgs(TRUE);stopifnot(length(args)==2L)
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
pkg<-normalizePath(args[1],winslash='/');run<-normalizePath(args[2],winslash='/')
dir.create(file.path(run,'parts'),recursive=TRUE,showWarnings=FALSE)
data<-function(name){
 if(name=='Table_S7_loading_simulation.csv')p<-file.path(pkg,'data/observations/loading_simulation.csv')
 else if(name=='gene_statistics.csv')p<-file.path(run,'plot_data/rna_degs/main_gene_statistics.csv')
 else p<-file.path(run,'plot_data/pre_inoculation',name)
 read.csv(p,check.names=FALSE,fileEncoding='UTF-8-BOM')
}
save_pdf<-function(name,width,height,plot_fun){
 cairo_pdf(file.path(run,'parts',name),width/25.4,height/25.4,family='sans',bg='white')
 on.exit(dev.off());plot_fun()
}
observed <- data("Figure1d_observations.csv")
summaries <- data("Figure1d_group_summaries.csv")
tests <- data("Figure1d_Welch_effects.csv")
rates <- data("Table_S7_loading_simulation.csv")
genes <- data("gene_statistics.csv")
groups <- c("stationary_1g", "clinorotation")
fill <- c("#C8E1EF", "#F4D7C2")
border <- c("#004F73", "#8D3E00")
point <- c("#0072B2", "#D55E00")

bar_panel <- function(metric, axis_title, labels, ticks, y_max, q_digits=3L) {
  d <- observed[observed$metric == metric, ]
  s <- summaries[summaries$metric == metric, ]
  s <- s[match(groups, s$gravity), ]
  q <- tests$fdr[tests$metric == metric]
  stopifnot(nrow(d)>0L, length(q)==1L, all(table(factor(d$gravity, levels=groups)) == s$n))
  par(mar=c(2.8, 3.35, 0.35, 0.3), mgp=c(2.05, 0.62, 0),
      tcl=-0.25, xaxs="i", yaxs="i", ps=10)
  plot.new()
  plot.window(c(0.48, 2.52), c(0, y_max))
  rect(c(1,2)-0.29, 0, c(1,2)+0.29, s$mean, col=fill, border=border, lwd=1)
  for (i in 1:2) {
    y <- d$value[d$gravity == groups[i]]
    x <- numeric(length(y))
    for (ix in split(seq_along(y), y)) {
      x[ix] <- i + if (length(ix) > 1L) seq(-0.15, 0.15, length.out=length(ix)) else 0
    }
    points(x, y, pch=16, col=point[i], cex=0.78)
  }
  arrows(1:2, s$mean-s$se, 1:2, s$mean+s$se, angle=90, length=.04, code=3, lwd=.8)
  axis(1, 1:2, labels=labels, cex.axis=.81)
  axis(2, ticks, las=1, cex.axis=.84)
  title(ylab=axis_title, cex.lab=.91)
  y_bracket <- max(d$value, s$mean+s$se) + 0.09*y_max
  segments(1, y_bracket, 2, y_bracket)
  segments(1:2, y_bracket, 1:2, y_bracket-.03*y_max)
  text(1.5, y_bracket + .045*y_max, sprintf(paste0("q = %.",q_digits,"f"),q), cex=.82)
}
save_pdf("cultivation.pdf", 119, 56, function() {
  par(mfrow=c(1,3))
  bar_panel("Root length (mm)", "Root length (mm)", c("Control","Chip"), seq(0,60,20), 66)
  bar_panel("Leaf area (cm²)", expression(Leaf~area~(cm^2)), c("Control","Chip"), seq(0,.20,.05), .253)
  bar_panel("Leaf fresh weight (mg)", "Leaf fresh weight (mg)", c("Control","Chip"), seq(0,6,2), 6.7)
})
save_pdf("loading.pdf", 66, 56, function() {
  par(mar=c(3.15,3.75,.35,.55), mgp=c(2.3,.65,0), tcl=-.25, ps=10)
  x <- rates$Inlet_rate_uL_per_min
  y <- rates$Simulated_outlet_arrival_time_s
  stopifnot(identical(x, c(.1,1,3,5,10,50,100,200)), y[x==50] == 22.68)
  plot(x,y,log="xy", xlim=c(.075,275),ylim=c(4,21000), axes=FALSE,
       xlab=expression(Inlet~rate~(mu*L~min^{-1})), ylab="Outlet-arrival time (s)",
       cex.lab=.92, type="n")
  axis(1,c(.1,1,10,100),labels=c("0.1","1","10","100"),cex.axis=.84)
  axis(2,c(10,100,1000,10000),labels=c("10","100","1,000","10,000"),las=1,cex.axis=.80)
  points(x,y,pch=16,col="#405B6B",cex=.8)
  points(x[x==50],y[x==50],pch=21,col="#004F73",bg="#C8E1EF",lwd=1.2,cex=1.4)
  text(40,69,"50 µL/min",adj=1,cex=.76,col="#004F73")
  text(40,40,"22.68 s",adj=1,cex=.76,col="#004F73")
  box(bty="l")
})
save_pdf("lateral_roots.pdf", 64, 69, function() {
  bar_panel("Lateral root number","Lateral roots per plant",c("1 g", "sµg"),seq(0,30,10),30,q_digits=4L)
})
save_pdf("deg_effects.pdf", 116, 69, function() {
  g <- genes[order(genes$logFC),]
  stopifnot(nrow(g)==8L, all(g$adj.P.Val < .05), sum(g$logFC<0)==7L)
  y <- nrow(g):1
  par(mar=c(2.8,6.0,1.35,.55), mgp=c(2.05,.6,0), tcl=-.25, ps=10)
  plot.new()
  plot.window(c(-4.4,8.2),c(.45,8.85),xaxs="i",yaxs="i")
  abline(v=0,col="#AAB2B8",lty=2)
  col <- ifelse(g$logFC<0,point[1],point[2])
  segments(0,y,g$logFC,y,col=col,lwd=1.4)
  points(g$logFC,y,pch=16,col=col,cex=.93)
  axis(2,at=y,labels=g$display_symbol,las=1,tick=FALSE,line=-.15,cex.axis=.83)
  axis(1,at=c(-4,-2,0,2,4,6),cex.axis=.84)
  title(xlab=expression(log[2]~fold~change~(s*mu*g~"vs"~1~g)),cex.lab=.92)
  text(7.4,8.73,"BH q",cex=.8,font=2,xpd=NA)
  text(7.4,y,sprintf("%.3f",g$adj.P.Val),cex=.83)
})
