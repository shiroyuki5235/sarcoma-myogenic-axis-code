# Recreate all submitted figures from S1 Data in R.
# Set the working directory to the extracted S1_Code folder.
# Run Module 00 and 01, then any figure module. No expression models are refitted.

#Module 00. Paths, packages and figure settings ----
options(stringsAsFactors = FALSE)
CODE_DIR <- normalizePath(Sys.getenv("SARCOMA_CODE_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)
S1_DATA_FILE <- file.path(dirname(CODE_DIR), "S1_Data.xlsx")
OUTPUT_DIR <- file.path(CODE_DIR, "Results", "Figures")
FIGURE_DPI <- 600L
if (!file.exists(S1_DATA_FILE)) stop("Set CODE_DIR to S1_Code and place S1_Data.xlsx beside that folder.")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl", repos = "https://cloud.r-project.org")
stopifnot(requireNamespace("readxl", quietly = TRUE), capabilities("tiff"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
blue <- "#285F7C"; orange <- "#B55E35"
group_colors <- c(blue, orange)
short_labels <- c("Dedifferentiated liposarcoma" = "DDLPS", "Undifferentiated/pleomorphic sarcoma" = "UPS",
                  "Rhabdomyosarcoma" = "RMS", "Leiomyosarcoma" = "LMS", "Synovial sarcoma" = "SS", "BCOR-CCNB3 sarcoma" = "BCOR-CCNB3")
palette_colors <- c("#2878A5", "#D76A24", "#469B61", "#BC4C51", "#8964A7", "#876655", "#BE72A1", "#777777", "#A49A30", "#17A2A4")

#Module 01. Load the reported data and check sample IDs ----
geo_records <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "GEO_records", skip = 4, n_max = 1176))
clustering_results <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "Clustering", skip = 4, n_max = 30))
signature_genes <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "Signature_genes", skip = 4, n_max = 34))
c1_reference_scores <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "C1_reference", skip = 4, n_max = 592))
pca_coordinates <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "PCA", skip = 4, n_max = 1258))
tcga_scores <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "TCGA_scores", skip = 4, n_max = 206))
auc_results <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "AUC_comparison", skip = 4, n_max = 3))
representative_metadata <- pca_coordinates[pca_coordinates$stage == "before_ComBat", ]
stopifnot(nrow(geo_records) == 1176L, sum(geo_records$Sample_Decision == "KEEP") == 1115L,
          nrow(representative_metadata) == 629L, sum(representative_metadata$original_group == "C2") == 37L,
          !anyDuplicated(representative_metadata$sampleid), nrow(tcga_scores) == 206L,
          !anyDuplicated(tcga_scores$PATIENT_ID), !anyDuplicated(tcga_scores$SAMPLE_ID),
          sum(signature_genes$in_revised) == 30L, nrow(c1_reference_scores) == 592L,
          !anyNA(pca_coordinates[, c("PC1", "PC2")]), !anyNA(tcga_scores$myogenic_score))
pca_before_ids <- pca_coordinates$sampleid[pca_coordinates$stage == "before_ComBat"]
pca_after_ids <- pca_coordinates$sampleid[pca_coordinates$stage == "after_ComBat"]
stopifnot(identical(pca_before_ids, pca_after_ids))
print(data.frame(dataset = c("retained_GEO", "representatives", "C1_reference", "TCGA"),
                 n = c(1115, nrow(representative_metadata), nrow(c1_reference_scores), nrow(tcga_scores))))

#Module 02. Fig 1: Study workflow ----
for (figure_format in c("pdf", "tiff")) {
  figure_file <- file.path(OUTPUT_DIR, paste0("Fig1.", ifelse(figure_format == "tiff", "tif", "pdf")))
  if (figure_format == "pdf") cairo_pdf(figure_file, width = 7.2, height = 5.8, family = "serif")
  if (figure_format == "tiff") tiff(figure_file, width = 7.2, height = 5.8, units = "in", res = FIGURE_DPI, compression = "lzw", type = "cairo")
  par(mar = c(0.8, 0.8, 0.8, 0.8), family = "serif", bg = "white")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  workflow_boxes <- data.frame(x = c(.07, .07, .07, .53, .07, .53), y = c(.81, .62, .38, .38, .08, .08),
                              w = c(.86, .86, .40, .40, .40, .40), h = c(.12, .13, .18, .18, .21, .21))
  workflow_text <- c("Nine predefined GEO series\n1,176 deposited records; 1,115 retained records",
    "Explicit cross-series links define 629 connected-record groups\nOne representative per group; patient-level independence not established",
    "Structure sensitivity\n20,287 genes; matched spaces\nSelected k changes after reduction",
    "Fixed C1/C2 comparison groups\nC1: 592; C2: 37 (36 LMS)\n30-gene myogenic signature",
    "GEO clinical annotations\nSource-specific recorded metastasis\nNo verified event-time analysis",
    "Independent TCGA-SARC projection\n28 available signature genes\n206 patients; primary cohort n = 174\nHistology discrimination and OS/DFS")
  for (box_index in seq_len(nrow(workflow_boxes))) {
    b <- workflow_boxes[box_index, ]; rect(b$x, b$y, b$x+b$w, b$y+b$h, col = "#F0F5F9", border = blue)
    text(b$x+b$w/2, b$y+b$h/2, workflow_text[box_index], cex = .80)
  }
  arrows(c(.50,.36,.64,.67,.73), c(.81,.62,.62,.38,.38), c(.50,.27,.73,.27,.73),
         c(.755,.565,.565,.295,.295), length = .07, col = blue)
  dev.off()
}

#Module 03. Fig 2: Matched-space clustering ----
clustering_conditions <- c("all_records__full_gene", "all_records__top_variable_PC",
 "dedup_fixed_correction__full_gene", "dedup_fixed_correction__top_variable_PC",
 "dedup_refit_correction__full_gene", "dedup_refit_correction__top_variable_PC")
condition_labels <- c("All records / full genes", "All records / PCs", "Representatives, fixed / full genes",
                      "Representatives, fixed / PCs", "Representatives, refit / full genes", "Representatives, refit / PCs")
stopifnot(setequal(unique(clustering_results$condition), clustering_conditions))
for (figure_format in c("pdf", "tiff")) {
  figure_file <- file.path(OUTPUT_DIR, paste0("Fig2.", ifelse(figure_format == "tiff", "tif", "pdf")))
  if (figure_format == "pdf") cairo_pdf(figure_file, width = 7.5, height = 4.6, family = "serif")
  if (figure_format == "tiff") tiff(figure_file, width = 7.5, height = 4.6, units = "in", res = FIGURE_DPI, compression = "lzw", type = "cairo")
  layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE), heights=c(4,1))
  par(mar = c(4,4.3,2.8,1), family = "serif", bg = "white", cex = .85, las = 1, bty = "l")
  for (panel in 1:2) {
    variable <- c("silhouette", "mean_conditional_ARI")[panel]
    plot(NA, xlim=c(2,6), ylim=range(clustering_results[[variable]])*c(.9,1.05), xaxt="n",
         xlab="Number of clusters (k)", ylab=c("Mean silhouette","Mean adjusted Rand index")[panel],
         main=c("A  Matched-space silhouette","B  Conditional subsampling")[panel])
    axis(1, at=2:6)
    for (condition_index in seq_along(clustering_conditions)) {
      d <- clustering_results[clustering_results$condition == clustering_conditions[condition_index], ]; d <- d[order(d$k), ]
      lines(d$k, d[[variable]], type="o", pch=16, cex=.55, col=palette_colors[condition_index], lwd=1.2)
    }
  }
  par(mar=c(0,0,0,0)); plot.new()
  legend("center", legend=condition_labels, col=palette_colors[1:6], pch=16, lty=1, ncol=2, bty="n", cex=.90)
  dev.off()
}

#Module 04. Fig 3: Fixed-group composition ----
for (figure_format in c("pdf", "tiff")) {
  figure_file <- file.path(OUTPUT_DIR, paste0("Fig3.", ifelse(figure_format == "tiff", "tif", "pdf")))
  if (figure_format == "pdf") cairo_pdf(figure_file, width=7.5, height=4.6, family="serif")
  if (figure_format == "tiff") tiff(figure_file, width=7.5, height=4.6, units="in", res=FIGURE_DPI, compression="lzw", type="cairo")
  par(mfrow=c(1,2), mar=c(4,6.2,2.8,.8), family="serif", bg="white", cex=.83, las=1, bty="l")
  for (panel in 1:2) {
    key <- c("histology", "gse_id")[panel]
    categories <- representative_metadata[[key]]; categories[is.na(categories)] <- "Unknown"
    counts <- table(categories, factor(representative_metadata$original_group, levels=c("C1","C2")))
    counts <- counts[order(rowSums(counts)), , drop=FALSE]
    labels <- rownames(counts); labels[labels %in% names(short_labels)] <- short_labels[labels[labels %in% names(short_labels)]]
    barplot(t(counts), beside=TRUE, horiz=TRUE, col=group_colors, border=NA, names.arg=labels,
            xlab="Representative records", main=c("A  Histology", "B  Representative source")[panel], cex.names=.90)
    if(panel==2) legend("bottomright", legend=c("C1","C2"), fill=group_colors, bty="n", cex=.9)
  }
  dev.off()
}

#Module 05. Fig 4: Within-source and within-LMS contrasts ----
fixed_signature <- signature_genes[signature_genes$in_revised, ]
marker_genes <- c("MYH11", "LMOD1", "CNN1", "MYLK", "MYOCD")
marker_contrasts <- fixed_signature[match(marker_genes, fixed_signature$GeneSymbol), ]
stopifnot(!anyNA(marker_contrasts$GeneSymbol), nrow(fixed_signature)==30L)
for (figure_format in c("pdf", "tiff")) {
  figure_file <- file.path(OUTPUT_DIR, paste0("Fig4.", ifelse(figure_format == "tiff", "tif", "pdf")))
  if (figure_format == "pdf") cairo_pdf(figure_file, width=7.5, height=4.5, family="serif")
  if (figure_format == "tiff") tiff(figure_file, width=7.5, height=4.5, units="in", res=FIGURE_DPI, compression="lzw", type="cairo")
  layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE), heights=c(5,1))
  par(mar=c(4.4,4.8,2.8,.8), family="serif", bg="white", cex=.85, las=1, bty="l")
  plot(marker_contrasts$all_histologies_logFC, 5:1+.06, pch=16, col=blue, xlim=c(0,6), ylim=c(.5,5.5), yaxt="n",
       ylab="", xlab="log2 fold change (C2 versus C1)", main="A  Within-source marker contrast")
  points(marker_contrasts$within_LMS_logFC, 5:1-.06, pch=16, col=orange)
  axis(2, at=5:1, labels=marker_genes, tick=FALSE)
  plot(fixed_signature$all_histologies_logFC, fixed_signature$within_LMS_logFC, pch=16, col=blue, cex=.75,
       xlab="log2 fold change: all histologies", ylab="log2 fold change: LMS only", main="B  All 30 fixed signature genes")
  par(mar=c(0,0,0,0));plot.new();legend("center",legend=c("All histologies (275 / 37)","LMS only (53 / 36)"),
    col=group_colors,pch=16,ncol=2,bty="n")
  dev.off()
}

#Module 06. Fig 5: Descriptive C1 reference panels ----
reference_panels <- c("Smooth_muscle_reference", "Fibroblast_ECM_reference", "T_cell_reference", "Cytotoxic_reference", "Myeloid_reference", "Endothelial_reference")
c1_reference_scores$histology[is.na(c1_reference_scores$histology)] <- "Unknown"
c1_histology_n <- table(c1_reference_scores$histology)
heatmap_histologies <- sort(names(c1_histology_n)[c1_histology_n>=5L])
heatmap_medians <- matrix(NA_real_, nrow=length(heatmap_histologies), ncol=6)
for (histology_index in seq_along(heatmap_histologies)) for (panel_index in seq_along(reference_panels)) {
  heatmap_medians[histology_index,panel_index] <- median(c1_reference_scores[c1_reference_scores$histology==heatmap_histologies[histology_index],reference_panels[panel_index]],na.rm=TRUE)
}
heatmap_labels <- heatmap_histologies
heatmap_labels[heatmap_labels %in% names(short_labels)] <- short_labels[heatmap_labels[heatmap_labels %in% names(short_labels)]]
heatmap_labels <- paste0(heatmap_labels," (n=",as.integer(c1_histology_n[heatmap_histologies]),")")
heatmap_colors <- colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(101)
for (figure_format in c("pdf", "tiff")) {
  figure_file <- file.path(OUTPUT_DIR, paste0("Fig5.", ifelse(figure_format == "tiff", "tif", "pdf")))
  if (figure_format == "pdf") cairo_pdf(figure_file, width=7.5, height=4.3, family="serif")
  if (figure_format == "tiff") tiff(figure_file, width=7.5, height=4.3, units="in", res=FIGURE_DPI, compression="lzw", type="cairo")
  layout(matrix(c(1,2),nrow=1),widths=c(6,1))
  par(mar=c(4.5,9.4,1,0.3),family="serif",bg="white",cex=.82,las=1)
  image(1:6,seq_along(heatmap_histologies),t(heatmap_medians[nrow(heatmap_medians):1,,drop=FALSE]),
        col=heatmap_colors,zlim=c(-1,1),axes=FALSE,xlab="",ylab="",useRaster=TRUE)
  axis(1,at=1:6,labels=c("Smooth\nmuscle","Fibroblast /\nECM","T cell","Cytotoxic","Myeloid","Endothelial"),tick=FALSE,cex.axis=.85)
  axis(2,at=seq_along(heatmap_histologies),labels=rev(heatmap_labels),tick=FALSE,cex.axis=.9)
  for (i in seq_len(nrow(heatmap_medians))) for (j in 1:6) text(j,nrow(heatmap_medians)+1-i,sprintf("%.2f",heatmap_medians[i,j]),cex=.85)
  par(mar=c(4.5,1,1,5.2))
  image(c(0,1),seq(-1,1,length.out=101),matrix(rep(seq(-1,1,length.out=101),each=2),nrow=2),
        col=heatmap_colors,axes=FALSE,xlab="",ylab="")
  axis(4,at=c(-1,-.5,0,.5,1),las=1,cex.axis=.85)
  mtext("Median panel score (gene-wise z scores)",side=4,line=3,cex=.75,las=0)
  dev.off()
}

#Module 07. Fig 6: Matched before/after PCA colored by source ----
pca_stages <- c("before_ComBat", "after_ComBat")
source_levels <- sort(unique(pca_coordinates$gse_id))
for (figure_format in c("pdf", "tiff")) {
  figure_file <- file.path(OUTPUT_DIR, paste0("Fig6.", ifelse(figure_format == "tiff", "tif", "pdf")))
  if (figure_format == "pdf") cairo_pdf(figure_file, width=7.5, height=4.6, family="serif")
  if (figure_format == "tiff") tiff(figure_file, width=7.5, height=4.6, units="in", res=FIGURE_DPI, compression="lzw", type="cairo")
  layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE),heights=c(5,1))
  par(mar=c(4.3,4.3,2.8,1),family="serif",bg="white",cex=.85,las=1,bty="l")
  for (panel in 1:2) {
    pca_subset <- pca_coordinates[pca_coordinates$stage==pca_stages[panel], ]
    plot(pca_subset$PC1,pca_subset$PC2,pch=16,cex=.45,col=adjustcolor(palette_colors[match(pca_subset$gse_id,source_levels)],alpha.f=.6),
      xlab=sprintf("PC1 (%.1f%%)",100*pca_subset$PC1_variance[1]),ylab=sprintf("PC2 (%.1f%%)",100*pca_subset$PC2_variance[1]),
      main=c("A  Before ComBat","B  After ComBat")[panel])
  }
  par(mar=c(0,0,0,0));plot.new();legend("center",legend=source_levels,col=palette_colors[seq_along(source_levels)],pch=16,ncol=4,bty="n",cex=.9)
  dev.off()
}

#Module 08. S1 Fig: Matched PCA colored by platform and histology ----
for (figure_format in c("pdf", "tiff")) {
  figure_file <- file.path(OUTPUT_DIR, paste0("S1_Fig.", ifelse(figure_format == "tiff", "tif", "pdf")))
  if (figure_format == "pdf") cairo_pdf(figure_file, width=7.5, height=8, family="serif")
  if (figure_format == "tiff") tiff(figure_file, width=7.5, height=8, units="in", res=FIGURE_DPI, compression="lzw", type="cairo")
  layout(matrix(c(1,2,3,3,4,5,6,6),nrow=4,byrow=TRUE),heights=c(4,0.7,4,1.4))
  par(family="serif",bg="white",cex=.85,las=1,bty="l")
  for (key in c("platform_id","histology")) {
    pca_groups <- pca_coordinates[[key]];pca_groups[is.na(pca_groups)] <- "Unknown"
    pca_levels <- sort(unique(pca_groups))
    stopifnot(length(pca_levels) <= length(palette_colors))
    for(panel in 1:2) {
      par(mar=c(4.3,4.3,2.8,1)); selection <- pca_coordinates$stage==pca_stages[panel]
      pca_subset <- pca_coordinates[selection, ]
      plot(pca_subset$PC1,pca_subset$PC2,pch=16,cex=.45,col=adjustcolor(palette_colors[match(pca_groups[selection],pca_levels)],alpha.f=.6),
        xlab=sprintf("PC1 (%.1f%%)",100*pca_subset$PC1_variance[1]),ylab=sprintf("PC2 (%.1f%%)",100*pca_subset$PC2_variance[1]),
        main=paste(if(key=="platform_id") c("A","B")[panel] else c("C","D")[panel],
                   c(" Before ComBat"," After ComBat")[panel]))
    }
    pca_labels <- pca_levels;pca_labels[pca_labels %in% names(short_labels)] <- short_labels[pca_labels[pca_labels %in% names(short_labels)]]
    par(mar=c(0,0,0,0));plot.new();legend("center",legend=pca_labels,col=palette_colors[seq_along(pca_levels)],pch=16,ncol=3,bty="n",cex=.9)
  }
  dev.off()
}

#Module 09. Fig 7: Fixed TCGA score and histology discrimination ----
histology_order <- c("STLMS","ULMS","DDLPS","UPS","MFS","SS","MPNST")
stopifnot(setequal(unique(tcga_scores$histology),histology_order),identical(as.integer(auc_results$n),c(174L,201L,206L)))
tcga_scores$histology_plot <- factor(tcga_scores$histology, levels=rev(histology_order))
tcga_boxplot <- boxplot(myogenic_score~histology_plot, data=tcga_scores, plot=FALSE)
# Linear-interpolated quartiles match the original figure; whiskers stay within 1.5 IQR.
for (histology_index in seq_along(tcga_boxplot$names)) {
  values <- tcga_scores$myogenic_score[tcga_scores$histology==tcga_boxplot$names[histology_index]]
  quartiles <- as.numeric(quantile(values, probs=c(.25,.5,.75), type=7))
  box_iqr <- quartiles[3]-quartiles[1]
  inside_whiskers <- values[values>=quartiles[1]-1.5*box_iqr & values<=quartiles[3]+1.5*box_iqr]
  tcga_boxplot$stats[,histology_index] <- c(min(inside_whiskers),quartiles,max(inside_whiskers))
}
for (figure_format in c("pdf", "tiff")) {
  figure_file <- file.path(OUTPUT_DIR, paste0("Fig7.", ifelse(figure_format == "tiff", "tif", "pdf")))
  if (figure_format == "pdf") cairo_pdf(figure_file, width=7.5, height=4.3, family="serif")
  if (figure_format == "tiff") tiff(figure_file, width=7.5, height=4.3, units="in", res=FIGURE_DPI, compression="lzw", type="cairo")
  par(mfrow=c(1,2),mar=c(4.5,5,2.8,1),family="serif",bg="white",cex=.85,las=1,bty="l")
  bxp(tcga_boxplot,horizontal=TRUE,outline=FALSE,boxfill="#D6E5EC",boxcol=blue,whiskcol=blue,staplecol=blue,
          ylab="",xlab="Mean gene-wise Z score\n(206-patient scaling)",main="A  Fixed myogenic score",medcol=orange)
  par(mar=c(4.5,7,2.8,1))
  plot(auc_results$AUC,3:1,pch=16,col=blue,xlim=c(.85,1),ylim=c(.5,3.5),yaxt="n",ylab="",
       xlab="AUC with 95% bootstrap interval",main="B  LMS discrimination")
  segments(auc_results$CI_low,3:1,auc_results$CI_high,3:1,col=blue)
  segments(auc_results$CI_low,3:1-.05,auc_results$CI_low,3:1+.05,col=blue)
  segments(auc_results$CI_high,3:1-.05,auc_results$CI_high,3:1+.05,col=blue)
  axis(2,at=3:1,labels=c("Primary\n174 / 53 LMS","+ uterine LMS\n201 / 80 LMS","+ MPNST\n206 / 80 LMS"),tick=FALSE,cex.axis=.85)
  dev.off()
}

#Module 10. Save session information ----
writeLines(capture.output(sessionInfo()),file.path(OUTPUT_DIR,"R_sessionInfo.txt"))
message("Figs 1–7 and S1 Fig saved: ",OUTPUT_DIR)
