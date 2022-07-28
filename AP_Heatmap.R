# "All purpose" tool for initial heatmap generation
## something quick, rote, but also tidy and visually appealing

library(readr)
library(ComplexHeatmap)

################################################################################
#
# load data
#
################################################################################

n_counts <- read_csv("Documents/BFX_proj/AP_Heatmap/_input/Chow_Genetics_normcounts.csv")
c_mtx <- as.matrix(n_counts[, 2:ncol(n_counts)])
rownames(c_mtx) <- n_counts$gene

meta <- read_csv("Documents/BFX_proj/AP_Heatmap/_input/meta.csv")

################################################################################
#
# prep data for heatmap
#
################################################################################

c_mtx <- t(scale(t(c_mtx)))

set.seed(415); Heatmap(c_mtx[1:10, ])

