gc()
library(ComplexHeatmap)
data_sub_2=matrix(data = c(0.74,0.67,0.35,0.71,0.85,
                           0.84,0.67,0.24,0.64,0.11,
                           0.23,0.18,0.13,0.87,0.64,
                           0.20,0.88,0.46,0.38,0.34,
                           0.78,0.12,0.4,0.72,0.85),
                  nrow = 5, ncol = 5)
colnames(data_sub_2) = c("Mouse", "CMS1", "CMS2", "CMS3", "CMS4")
ht = Heatmap(as.matrix(data_sub_2), cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        col=c("#87CEFA", "white", "#CC2121"),
        #left_annotation = ha2,
        row_names_gp = gpar(fontsize = 14),
        column_names_rot = 45,
        show_heatmap_legend = F,
        column_names_centered = T)

pdf("/home/shimw/project/conservation/showcase.pdf", width = 2.1, height = 2.5)
draw(ht)
dev.off()

pdf("/home/shimw/project/conservation/showcase.pdf", width = 2.3, height = 2.55)
draw(ht, padding= unit(c(2, 9, 2, 2), "mm"))
dev.off()

