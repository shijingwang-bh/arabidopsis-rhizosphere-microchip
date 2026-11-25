args<-commandArgs(TRUE);stopifnot(length(args)==2L)
pkg<-normalizePath(args[1],winslash='/');run<-normalizePath(args[2],winslash='/')
out<-file.path(run,'figures');dir.create(out,recursive=TRUE,showWarnings=FALSE)
plot_data<-file.path(run,'plot_data')
aliases<-jsonlite::fromJSON(file.path(pkg,'metadata/plot_data_aliases.json'))
source(file.path(pkg,'code/plot_lib/common.R'))
source(file.path(pkg,'code/plot_lib/platform.R'))
source(file.path(pkg,'code/plot_lib/chemistry.R'))
source(file.path(pkg,'code/plot_lib/RNA.R'))
write_source <- function(...) invisible(NULL)
save_panels <- function(name, width, height, draw) {
  cairo_pdf(file.path(out, paste0(name, ".pdf")), width / 25.4,
            height / 25.4, family = font, bg = "white")
  on.exit(dev.off())
  grid.newpage()
  draw()
}
capture_panels <- function(draw) {
  code <- as.list(body(draw))[-1]
  code <- Filter(function(x) {
    !is.call(x) || !as.character(x[[1]])[1] %in% c("writeLines", "cat")
  }, code)
  body(draw) <- as.call(c(list(quote(`{`)), code, list(quote(environment()))))
  draw()
}
export <- function(fig, w, h, fun) invisible(NULL)

section_label <- function(label, y) {
  txt(label, 0, y, 10.5, TRUE, col = "#344B5C")
  grid.lines(x = mm(c(0, 190)), y = unit(1, "npc") - mm(y + 5),
             gp = gpar(col = "#CBD5DB", lwd = 0.7))
}


live<-capture_panels(draw_Figure_2)
growth<-capture_panels(draw_Figure_3)
chem_global<-capture_panels(draw_Figure_4)
chem_interaction<-capture_panels(draw_Figure_5)
rna_pca<-capture_panels(draw_Figure_S2)$p
rna_go<-capture_panels(draw_Figure_S8)$p
live_body<-as.list(body(draw_Figure_2))
original_export<-live_body[[length(live_body)]]
microscopy<-as.list(body(eval(original_export[[5]],live)))[[2]]
save_panels("Figure_2_bottom", 190, 74, function() {
  place(0, 0, 74, 74, function() panel(rna_pca, "c", "Leaf RNA sample PCA", 74, 74))
  go <- rna_go + guides(colour = "none") + theme(legend.position = "bottom")
  place(78, 0, 112, 74, function() panel(go, "d", "Selected functional themes", 112, 74))
})

save_panels("Figure_3", 190, 241, function() {
  section_label("Live-root co-culture in the chip", 0)
  place(0, 8, 190, 99, function() {
    eval(microscopy, live)
    place(107, 0, 83, 99, function() panel(live$pf, "b", "Root-zone fluorescence", 83, 99))
  })
  place(0, 110, 190, 51, function() {
    counts <- live$pc + labs(y = expression(atop("Concentration", (10^8~cells~mL^{-1}))))
    panel(counts, "c", "Bacterial concentration in chip-channel liquid", 190, 51)
  })
  section_label("Growth in medium recovered from uninoculated chips", 168)
  place(0, 177, 94, 55, function() panel(growth$curve(3), "d", "Day-3 collection", 94, 55))
  place(98, 177, 92, 55, function() panel(growth$curve(7), "e", "Day-7 collection", 92, 55))
  for (i in seq_along(growth$cols)) {
    x <- c(9, 71, 138)[i]
    grid.lines(x = mm(c(x, x + 7)), y = unit(1, "npc") - mm(235),
               gp = gpar(col = growth$cols[i], lwd = 1.7))
    txt(growth$labs[i], x + 9, 233, 8.5)
  }
})

save_panels("Figure_4", 190, 188, function() {
  place(0, 0, 93, 81, function() panel(chem_global$pa, "a", "Root-exudate profiles", 93, 81))
  place(97, 0, 93, 81, function() panel(chem_global$pb, "b", "Factor effects", 93, 81))
  heat <- chem_global$pc + theme(legend.text = element_text(size = 8),
                                 axis.text.x = element_text(size = 8)) +
    guides(fill = guide_colourbar(direction = "horizontal", barwidth = mm(57),
                                  barheight = mm(3), title.position = "top"))
  place(0, 86, 103, 98, function() panel(heat, "c", "Inoculation responses", 103, 98))
  place(107, 86, 83, 98, function() panel(chem_interaction$pa, "d", "Interaction changes over time", 83, 98))
})

save_panels("Figure_5", 190, 226, function() {
  section_label("Day-7 root-exudate chemistry in independent culture dishes", 0)
  place(0, 9, 190, 66, function() {
    panel(chem_interaction$pb, "a", "Day-7 interaction candidates", 190, 66)
  })
  place(0, 79, 94, 66, function() {
    place(0, 1, 94, 64, function() drawp(chem_interaction$chem(1)))
    txt("b", 0, 0, 12, TRUE)
  })
  place(98, 79, 92, 66, function() {
    place(0, 1, 92, 64, function() drawp(chem_interaction$chem(2)))
    txt("c", 0, 0, 12, TRUE)
  })
  section_label("Plant growth after 7 days of co-culture in the chip", 149)
  place(0, 158, 94, 61, function() panel(growth$plots[[1]], "d", "Mature-root diameter", 94, 61))
  place(98, 158, 92, 61, function() panel(growth$plots[[2]], "e", "Leaf fresh weight", 92, 61))
  for (i in seq_along(pal)) {
    x <- c(69, 111)[i]
    grid.rect(x = mm(x), y = unit(1, "npc") - mm(222), width = mm(5), height = mm(3),
              gp = gpar(fill = pal_fill[i], col = pal_dark[i]))
    txt(glab[i], x + 5, 220, 8.5)
  }
})

# Current supplementary numbering. PCA is already Figure 2c and has no
# duplicate standalone supplementary output.
export<-function(fig,w,h,fun){
 mapping<-c(Figure_S3='Figure_S2',Figure_S8='Figure_S3',Figure_S5a='Figure_S5a',
            Figure_S5b='Figure_S5b',Figure_S6='Figure_S6')
 save_panels(mapping[[fig]],w,h,fun)
}
draw_Figure_S3();draw_Figure_S8();draw_Figure_S5a();draw_Figure_S5b();draw_Figure_S6()
source(file.path(pkg,'code/plot_day3_examples.R'))
capture.output(sessionInfo(),file=file.path(run,'plotting_sessionInfo.txt'))
cat('PASS: current main and supplementary statistical plots\n')
