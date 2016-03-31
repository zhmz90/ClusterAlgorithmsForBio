library(pheatmap)
# 
# test = matrix(rnorm(200), 20, 10)
# test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
# test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
# test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
# colnames(test) = paste("Test", 1:10, sep = "")
# rownames(test) = paste("Gene", 1:20, sep = "")
# 
# pheatmap(test)
# pheatmap(test, kmeans_k = 2)
# pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
# pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# pheatmap(test, cluster_row = FALSE)
# pheatmap(test, legend = FALSE)

matrix <- read_csv("../data/matrix.csv",col_names = FALSE)
matrix <- t(matrix)
#matrix <- as.data.frame(matrix)
sample_names <- read_csv("../data/sample_names.csv",col_names =FALSE)
feature_names <- read_csv("../data/gene_mutations.csv",col_names = FALSE)
rownames(matrix) <- feature_names
colnames(matrix) <- sample_names

pheatmap(matrix,cluster_rows = TRUE)



