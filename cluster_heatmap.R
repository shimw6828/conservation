################################################
##大热图
##对所有基因，分cluster画热图
##2021/12/05


##############################################
##读入数据
human_data <- read.csv("/home/shimw/project/conservation/filter_human_data.csv",row.names = 1, check.names = F)
human_info <- read.csv("/home/panxl/CRC2/filter_human_info.csv",row.names = 1, check.names = F)
mouse_data <- read.csv("/home/shimw/project/conservation/filter_mouse_data.csv",row.names = 1, check.names = F)
mouse_info <- read.csv("/home/panxl/CRC2/filter_mouse_info.csv",row.names = 1, check.names = F)
cluster_type = read.csv("/home/panxl/CRC2/all_cluater",check.names = F)
cluster_type = cluster_type[cluster_type$cluster!="Cluster 5: Not-significant in both species",]
human_data <- human_data[cluster_type$Gene.stable.ID,]%>%t()%>%
  scale()%>%t()%>%as.data.frame()
mouse_data <- mouse_data[cluster_type$Mouse.gene.stable.ID,]%>%t()%>%
  scale()%>%t()%>%as.data.frame()
row.names(mouse_data) <- cluster_type$Gene.stable.ID


human_info$species = "Human"
mouse_info$species = "Mouse"
names(mouse_info)[2] = "sample"
all_info = rbind(human_info,mouse_info)
all_data <- cbind(human_data,mouse_data)



library(ComplexHeatmap)
library(circlize)
colors=list(Group = c('Normal' = '#2fa1dd', 'Tumor' = '#f87669'),
            species=c('Human' = '#2fa1dd', 'Mouse' = '#f87669'))

# species = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
#                      labels = c("Human","Mouse"),
#                      labels_gp = gpar(col = "white", fontsize = 12))


library(GetoptLong)  # for the function qq()
group_block_anno = function(group, empty_anno, gp = gpar(), 
                            label = NULL, label_gp = gpar()) {
  
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gp)
  if(!is.null(label)) {
    grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp)
  }
}


top_annotation = HeatmapAnnotation(
                                   df = all_info[,c("Species","Group")],
                                   col = colors,
                                   which = 'col'
)


gene_cluster = factor(cluster_type$cluster,
                      levels = c("Cluster 1: Up-regulated in both species",
                                 "Cluster 2: Up-regulated in mouse reverse in human",
                                 "Cluster 3: Down-regulated in both species",
                                 "Cluster 4: Up-regulated in human reverse in mouse"
                      ))
row_annotation = rowAnnotation(cluster = gene_cluster)

# colors=list(Group = c('Normal' = '#2fa1dd', 'Tumor' = '#f87669'))
# top_annotation = HeatmapAnnotation(empty = anno_empty(border = FALSE, height = unit(5, "mm")),
#                                    Group = all_info[,c("Group")],
#                                    col = colors
# )



species = all_info$species
species = factor(species,levels = c("Human","Mouse"))
Group = all_info$Group
Group = factor(Group,levels = c("Normal","Tumor"))

kk = data.frame("species"=species,"Group"=Group)


m = Heatmap(all_data,
            top_annotation = top_annotation,
            left_annotation = row_annotation,
            name = "Expression",
            column_split = kk,
            cluster_column_slices = FALSE,
            row_split = gene_cluster,
            cluster_row_slices = FALSE,
            border = F,
            column_gap  = unit(c(0, 1,0), "mm"),
            show_row_names = F,
            show_column_names = F,
            cluster_rows = F,
            column_title = NULL,
            row_title = NULL)

draw(m,  annotation_legend_side="right",merge_legend=TRUE)








#########################################################
##调整
human_data <- read.csv("/home/shimw/project/conservation/filter_human_data.csv",row.names = 1, check.names = F)
human_info <- read.csv("/home/panxl/CRC2/filter_human_info.csv",row.names = 1, check.names = F)
mouse_data <- read.csv("/home/shimw/project/conservation/filter_mouse_data.csv",row.names = 1, check.names = F)
mouse_info <- read.csv("/home/panxl/CRC2/filter_mouse_info.csv",row.names = 1, check.names = F)
cluster_type = read.csv("/home/panxl/CRC2/all_cluater",check.names = F)
cluster_type = cluster_type[cluster_type$cluster!="Cluster 5: Not-significant in both species",]
cms_type = read.csv("/home/panxl/CRC2/human_cmstypes.csv",check.names = F)

human_info$Species = "Human"
mouse_info$Species = "Mouse"
names(mouse_info)[2] = "sample"

all_info = rbind(human_info,mouse_info)
all_info = left_join(all_info,cms_type,by = c("sample"="barcode"))
all_info$type[is.na(all_info$type)] = "Unknown"

all_info = arrange(all_info,type,Species)

human_data <- human_data[cluster_type$Gene.stable.ID,]%>%t()%>%
  scale()%>%t()%>%as.data.frame()
mouse_data <- mouse_data[cluster_type$Mouse.gene.stable.ID,]%>%t()%>%
  scale()%>%t()%>%as.data.frame()
row.names(mouse_data) <- cluster_type$Gene.stable.ID




all_data <- cbind(human_data,mouse_data)
all_data = all_data[,all_info$sample]





library(ComplexHeatmap)
library(circlize)
library(GetoptLong)  # for the function qq()
group_block_anno = function(group, empty_anno, gp = gpar(), 
                            label = NULL, label_gp = gpar()) {
  
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gp)
  if(!is.null(label)) {
    grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = gpar(col = "white",fontfamily="sans",fontface = "bold"))
  }
}
#8f7eab
top_annotation = HeatmapAnnotation(empty = anno_empty(border = FALSE, height = unit(6, "mm")),
                                   Group = all_info[,c("Group")],
                                   Type = all_info[,c("type")],
                                   col = list(Group = c('Normal' = '#b9e0dc', 'Tumor' = '#f18772'),
                                              Type = c('Unknown' = '#bbbcb9', 'CMS1' = '#b97573',
                                                        'CMS2' = '#686954', 'CMS3' = '#5b5e66', 'CMS4' = '#816147')),
                                   simple_anno_size = unit(3,"mm"),
                                   annotation_legend_param = list(
                                     Group = list(title = "Group",at = c("Normal","Tumor"),labels = c("Normal","Tumor")),
                                     Type = list(title = "Type",at = c("CMS1","CMS2","CMS3","CMS4","Unknown"),
                                                 labels = c("CMS1","CMS2","CMS3","CMS4","Unknown"))
                                   ),
                                   annotation_name_gp = gpar(fontsize = 8,fontfamily="sans"),
                                   # gp = gpar(lwd = 0.5),
                                   gap = unit(c(2,0), "mm")
)

species = all_info$Species
species = factor(species,levels = c("Human","Mouse"))
Group = all_info$Group
Group = factor(Group,levels = c("Normal","Tumor"))

kk = data.frame("species"=species,"Group"=Group)
gene_cluster = factor(cluster_type$cluster,
                      levels = c("Cluster 1: Up-regulated in both species",
                                 "Cluster 2: Up-regulated in mouse reverse in human",
                                 "Cluster 3: Down-regulated in both species",
                                 "Cluster 4: Up-regulated in human reverse in mouse"
                      ))
row_annotation = rowAnnotation(cluster = anno_block(gp = gpar(fill = c("#CB9C7A","#8696a7",
                                                                       "#CDB97D","#7b8b6f"),
                                                              col = 0),
                                                          labels = c("Cluster1", "Cluster2", "Cluster3",
                                                                     "Cluster4"), 
                                                          labels_gp = gpar(fontsize = 10),
                                                    width = unit(6, "mm")
                                                    ))





ht = Heatmap(all_data,
        top_annotation = top_annotation,
        left_annotation = row_annotation,
        name = "Expression",
        column_split = kk,
        cluster_column_slices = FALSE,
        row_split = gene_cluster,
        cluster_row_slices = FALSE,
        border = F,
        column_gap  = unit(c(0, 1,0), "mm"),
        show_row_names = F,
        show_column_names = F,
        cluster_rows = F,
        show_column_dend = F,
        column_title = NULL,
        row_title = NULL,
        cluster_column = F
        )



cairo_pdf(filename = "/home/shimw/project/conservation/cluster_heatmap.pdf",height=6.5,width = 5.5)
draw(ht,  annotation_legend_side="right",merge_legend=TRUE)
group_block_anno(1:2, "empty", gp = gpar(fill = "#a44a42",lwd = 0.5),label = "Human")
group_block_anno(3:4, "empty", gp = gpar(fill = "#d1af5b",lwd = 0.5),label = "Mouse")
dev.off()


