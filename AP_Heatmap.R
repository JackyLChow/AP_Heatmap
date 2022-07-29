# "All purpose" tool for initial heatmap generation
## something quick, rote, but also tidy and visually appealing

library(readr)
library(ComplexHeatmap)
library(matrixStats)

library(RColorBrewer)
library(pals)
library(circlize)

setwd("~/Documents/BFX_proj/AP_Heatmap")

################################################################################
#
# load data
#
################################################################################

n_counts <- read_csv("~/Documents/BFX_proj/AP_Heatmap/_input/Chow_Genetics_normcounts.csv")
c_mtx <- as.matrix(n_counts[, 2:ncol(n_counts)])
rownames(c_mtx) <- n_counts$gene

meta <- read_csv("~/Documents/BFX_proj/AP_Heatmap/_input/meta.csv")

################################################################################
#
# prep data for heatmap
#
################################################################################

### counts ---
h_mtx <- c_mtx[order(rowVars(c_mtx), decreasing = T)[1:50], ]
h_mtx <- t(scale(t(h_mtx)))
#h_mtx <- c_mtx[c("STE12", "DIG1", "RAS2", "RTG3", "SPT8", "MUC1", "RIM101"), ]

### heatmap annotation ---
# heatmap colors
col_htmp <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow2"))

# annotation colors
col_fg <- c("norm" = "grey75", "low" = "grey30", "high" = "black")
col_st <- brewer.pal(length(unique(meta$strain)), "Set1"); names(col_st) <- unique(meta$strain)

anno_col <- HeatmapAnnotation(strain = meta$strain,
                              fg = meta$filamentous,
                              col = list(strain = col_st,
                                         fg = col_fg))

################################################################################
#
# bake the heatmap
#
################################################################################

set.seed(415); hm_ <- draw(Heatmap(h_mtx,
                                   col = col_htmp,
                                   name = "z-score",
                                   height = unit(2.5, "mm") * nrow(h_mtx), width = unit(2.5, "mm") * ncol(h_mtx),
                                   show_row_dend = T, show_column_dend = T,
                                   row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
                                   top_annotation = anno_col,
                                   row_km = 2,
                                   column_km = 2))

png("_output/hm1.png", height = nrow(h_mtx) * 20, width = ncol(h_mtx) * 20 + 100, res = 100)
hm_
dev.off()

r_o <- row_order(hm_)
for(cl_ in names(r_o)){
  if(is.list(r_o)){
    r_o[[cl_]] <- rownames(h_mtx)[r_o[[cl_]]]
  } else {
    r_o <- rownames(h_mtx)[r_o]
  }
}
saveRDS(r_o, "_output/r_o.rds")

c_o <- column_order(hm_)
for(cl_ in names(c_o)){
  if(is.list(c_o)){
    c_o[[cl_]] <- colnames(h_mtx)[c_o[[cl_]]]
  } else {
    c_o <- colnames(h_mtx)[c_o]
  }
}
saveRDS(c_o, "_output/c_o.rds")


