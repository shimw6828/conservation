########################################### 1.批次矫正热图 ########################################################
######################1.1 human
human_data <- read.csv("/home/shimw/project/conservation/filter_human_data.csv",row.names = 1, check.names = F)
# human_info <- read.csv("/home/shimw/project/conservation/filter_human_info.csv",row.names = 1, check.names = F)
# human_info$Group=gsub("native","Normal",human_info$Group)
# human_info$Group=gsub("tumor","Tumor",human_info$Group)
# human_info$DataSet=gsub("ourdata","GSE156451",human_info$DataSet)
# write.csv(human_info,file = "/home/panxl/CRC2/filter_human_info.csv")
human_info <- read.csv("/home/panxl/CRC2/filter_human_info.csv",row.names = 1, check.names = F)
hh=cor(human_data, method = "spearman")
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer") 
display.brewer.all()
display.brewer.pal(n =11, name = 'Set3')#'Paired' "RdBu"
color_list1 = brewer.pal(n = 11, name = "Set3")
ha1 = HeatmapAnnotation(DataSet=factor(human_info[,c("DataSet")]),
                        Group = factor(human_info[,c("Group")], levels = c("Normal", "Tumor")),
                        col = list(Group = c("Normal"=color_list1[5],"Tumor"=color_list1[4] ), 
                                   DataSet=c("GSE156451"=color_list1[6],"TCGA"=color_list1[3])))#c("CY"=color_list1[4], "SG"="#65c294","DC"="#e0861a")))#,"CY"= color_list1[1], "SG" = color_list1[2],"DC"= color_list1[6]

ht1 = Heatmap(hh, cluster_rows = T,
           cluster_columns = T,
           show_row_dend =F,
           show_column_dend =F,
           name = "Correlation",
           show_column_names = F,
           show_row_names = F,
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
           col= c("#87CEFA", "white", "#CC2121"),
           top_annotation = ha1)
           #left_annotation = ha2)#"navy","white","#ae545e"
cairo_pdf(filename = "/home/shimw/project/conservation/human_cluster_heatmap.pdf",height=6,width = 5.5)
draw(ht1,  annotation_legend_side="right",merge_legend=TRUE)
dev.off()
######################1.2 mouse
mouse_data <- read.csv("/home/shimw/project/conservation/filter_mouse_data.csv",row.names = 1, check.names = F)
mm=cor(mouse_data, method = "spearman")
# mouse_info <- read.csv("/home/shimw/project/conservation/filter_mouse_info.csv",row.names = 1)
# mouse_info$DataSet=gsub("ourdata" ,"GSE178145",mouse_info$DataSet)
# write.csv(mouse_info,file = "/home/panxl/CRC2/filter_mouse_info.csv")
#write.csv(Supplementary_Table,file = "/home/panxl/CRC2/Supplementary_Table.csv")
mouse_info <- read.csv("/home/panxl/CRC2/filter_mouse_info.csv",row.names = 1, check.names = F)
# annotation_mcol = data.frame(
# Group = factor(mouse_info[,c("Group")]),DataSet=factor(mouse_info[,c("DataSet")]),row.names=rownames(mouse_info))
# mm=cor(mouse_data, method = "spearman")
# library(pheatmap)
# detach("package:ComplexHeatmap", unload=TRUE)
# pheatmap::pheatmap(mm, annotation_col = annotation_mcol, show_rownames = F,
#                    show_colnames = F, clustering_method="ward.D")
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer") 
display.brewer.all()
display.brewer.pal(n =8, name = 'Dark2')#'Paired' "RdBu"
color_list2 = brewer.pal(n = 8, name = "Dark2")
display.brewer.pal(n =8, name = 'Dark2')#'Paired' "RdBu"
color_list2 = brewer.pal(n = 8, name = "Dark2")
display.brewer.pal(n =11, name = 'Set3')#'Paired' "RdBu"
color_list1 = brewer.pal(n = 11, name = "Set3")

ha2 = HeatmapAnnotation(DataSet=factor(mouse_info[,c("DataSet")]),
                        Group = factor(mouse_info[,c("Group")]),
                        col = list(Group = c("Normal"=color_list1[5],"Tumor"=color_list1[4] ), 
                                   DataSet=c("GSE102648"=color_list1[1],"GSE111615"=color_list1[3],
                                             "GSE116374"=color_list1[8],"GSE124029"=color_list1[9],
                                             "GSE132814"=color_list1[10], "GSE136029" =color_list1[11],
                                             "GSE142154"="#f8ebd8","GSE152032"=color_list2[1], 
                                             "GSE155777"=color_list2[2], "GSE158276"=color_list2[3], 
                                             "GSE164960"=color_list2[4], "GSE166277"=color_list2[5],  
                                             "GSE57533"=color_list2[6],  "GSE60837"=color_list2[7],
                                             "GSE97013"=color_list2[8],  "GSE98496"="#a27e7e","GSE178145"="#7b8b6f")))#c("CY"=color_list1[4], "SG"="#65c294","DC"="#e0861a")))#,"CY"= color_list1[1], "SG" = color_list1[2],"DC"= color_list1[6]
#colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
Heatmap(mm, cluster_rows = F,
           cluster_columns = F,
           show_column_names = F,
           show_row_names = F,
           col=c("#87CEFA", "white", "#CC2121"),
           top_annotation = ha2)
           #left_annotation = ha2)#"navy","white","#ae545e"
ht2 = Heatmap(mm, cluster_rows = T,
             cluster_columns = T,
             show_row_dend =F,
             show_column_dend =F,
             name = "Correlation",
             show_column_names = F,
             show_row_names = F,
             clustering_method_rows = "ward.D2",
             clustering_method_columns = "ward.D2",
             col= c("#87CEFA", "white", "#CC2121"),
             top_annotation = ha2)
cairo_pdf(filename = "/home/shimw/project/conservation/mouse_cluster_heatmap.pdf",height=4.5,width = 5.9)
###有点短，即使draw也不能放在一起，所以直接不draw更好看
ht2
dev.off()



########################################  2.人和小鼠差异表达 #########################################
######################2.1 human
library(limma)
human_data <- read.csv("/home/shimw/project/conservation/filter_human_data.csv",row.names = 1, check.names = F)
human_info <- read.csv("/home/panxl/CRC2/filter_human_info.csv",row.names = 1, check.names = F)
fit=lmFit(human_data,design=model.matrix(~human_info$Group))  
fit=eBayes(fit) 
options(digits = 4) 
topTable(fit,coef=2,adjust='BH') 
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
#write.csv(deg,file = "/home/panxl/CRC2/human_DEG.csv")
######################2.2 mouse
mouse_data <- read.csv("/home/shimw/project/conservation/filter_mouse_data.csv",row.names = 1, check.names = F)
mouse_info <- read.csv("/home/panxl/CRC2/filter_mouse_info.csv",row.names = 1, check.names = F)
fit=lmFit(mouse_data,design=model.matrix(~mouse_info$Group))  
fit=eBayes(fit) 
options(digits = 4) 
topTable(fit,coef=2,adjust='BH') 
deg=topTable(fit,coef=2,adjust='BH',number = Inf)#22413
#write.csv(deg,file = "/home/panxl/CRC2/mouse_DEG.csv")
########################################  3.DEG 火山图#########################################
rm(list=ls())
dev.off()
#画vocano plot看差异基因
all_mouse_DEG=read.csv("/home/panxl/CRC2/mouse_DEG.csv",row.names = 1)#22413
ensembl_ID_change=read.csv("/home/panxl/mouse_gene_9.txt.txt")#55416
a=ensembl_ID_change[ensembl_ID_change$Gene.stable.ID%in%rownames(all_mouse_DEG),]#22194
b=all_mouse_DEG[rownames(all_mouse_DEG)%in%ensembl_ID_change$Gene.stable.ID,]#
b$Gene.stable.ID=rownames(b)
all_mouse_DEG_changeID=left_join(a,b,by="Gene.stable.ID")#21047
all_mouse_DEG=all_mouse_DEG_changeID
#标记基因上调或者下调
Data =all_mouse_DEG[,c(1,2,3,6,7)]
Data$Change = ifelse(Data$P.Value < 0.05 & abs(Data$logFC) >= 1, 
                     ifelse((Data$logFC)> 1 ,'up-regulated','down-regulated'),'not-significant')
table(Data$Change )
# down-regulated not-significant    up-regulated 
# 2043           18246            1905 
Data_gene=Data[Data$P.Value<0.05&abs(Data$logFC)>6,]#4
#作图
library(tidyverse)
library(ggrepel)  #标签用
library(ggplot2) 
#添加几何对象  geom_point散点图，将direction映射给点颜色  aes映射颜色
p=ggplot(data = Data,aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = Change))+#,size = abs(logFC))) + 
  scale_color_manual(values=c("#228B22","#A9A9A9","#DC143C"))+
  ggrepel::geom_label_repel( 
    aes(label =Data_gene$Gene.name),
    data = Data_gene,
    color="black",
    size = 6,
    box.padding=unit(0.5, "lines"), point.padding=unit(0.5, "lines"), 
    segment.color = "#A9A9A9", segment.size = 1,
    arrow = arrow(length=unit(0.01, "npc")),force = 1, max.iter = 3e3, 
  )+
  ylab('-log10 (P_value)')+#修改y轴名称
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-1,1),lty=5,lwd=0.5) +#添加横线|FoldChange|>1
  geom_hline(yintercept = -log10(0.05),lty=5,lwd=0.5)+ 
  theme_classic(base_line_size = 1)+# 坐标轴的粗细
  theme(axis.title.x = element_text(size =30, 
                                    color = "black"),
        #face = "italic"),
        axis.title.y = element_text(size = 30,
                                    color = "black",
                                    #face = "italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 25),
        #face = "bold"),
        axis.text.x = element_text(size = 30,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   #face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 30,  
                                   color = "black",
                                   #face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p
#ggsave(p, file='/home/panxl/CRC2/volcano_deg.pdf', width=10, height=10) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
######################################### 4.1Mouse kegg pathway条形图（20）################################################################
all_mouse_DEG=read.csv("/home/panxl/CRC2/mouse_DEG.csv",row.names = 1)#22520
ensembl_ID_change=read.csv("/home/panxl/mouse_gene_9.txt.txt")#55416
a=ensembl_ID_change[ensembl_ID_change$Gene.stable.ID%in%rownames(all_mouse_DEG),]#22297
b=all_mouse_DEG[rownames(all_mouse_DEG)%in%ensembl_ID_change$Gene.stable.ID,]#
b$Gene.stable.ID=rownames(b)
all_mouse_DEG_changeID=left_join(a,b,by="Gene.stable.ID")#22297
all_mouse_DEG=all_mouse_DEG_changeID
#标记基因上调或者下调
res =all_mouse_DEG[,c(1,2,3,6,7)]
res_up=all_mouse_DEG[all_mouse_DEG$logFC>=1 &all_mouse_DEG$P.Value<0.05,]# 1905   #2188
res_down=all_mouse_DEG[all_mouse_DEG$logFC<=-1 &all_mouse_DEG$P.Value<0.05,]#2043   #2079
#BiocManager::install("org.Mm.eg.db") #这个包里存有人的注释文件
# 载入包dian
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Mm.eg.db)
#作图前处理——提gene ID --> 转换为ENTREZID
##up gene
DEG.gene_ENSEMBL_up = as.character(res_up$Gene.stable.ID) #获得基因 gene ID
DEG.entrez_id_up = mapIds(x = org.Mm.eg.db,
                          keys = DEG.gene_ENSEMBL_up,
                          keytype = "ENSEMBL",
                          column = "ENTREZID")#
DEG.entrez_id_up = na.omit(DEG.entrez_id_up)#去除NA 1459
kegg_up <- enrichKEGG(gene       = DEG.entrez_id_up,
                      keyType       = "kegg",
                      organism      = "mmu",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)
kegg_up_result=kegg_up@result
#write.csv(kegg_up_result,file = "/home/panxl/CRC2/mouse_DEG_kegg_up.txt")
library(ggplot2)
###up kegg pathway
kegg <- read.csv("/home/panxl/CRC2/mouse_DEG_kegg_up.txt",stringsAsFactors=F)
#对富集结果按照p.adjust进行从小到大排序，保证最显著的通路在前
kegg <- kegg[order(kegg$pvalue),]
kegg <- kegg[1:20,]
top20 <- data.frame(kegg$Description,kegg$Count ,kegg$pvalue,kegg$p.adjust)
colnames(top20) <- c("Description","count","pvalue","padj")
top20=top20%>%mutate(Description = fct_reorder(Description, -pvalue))
p1=ggplot(data=top20,aes(x=Description,y=-log10(pvalue),fill=padj))+#fill=padj fill颜色填充，使用连续值padj 
  geom_bar(stat="identity") + coord_flip()+#coord_flip()颠倒坐标轴
  labs(x="",y="-log10 (P_value)",title="KEGG Up Pathway")+
  scale_fill_gradient(low="#DC143C",high="#DC143C")+
  theme_classic()+
  theme(axis.title.x = element_text(size =22, 
                                    color = "black"),
        #face = "italic"),
        axis.title.y = element_text(size = 25,
                                    color = "black",
                                    #face = "italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        title = element_text(size = 22,
                             color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 17 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   #face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 19,  
                                   color = "black",
                                   vjust = 0.4))
p1
#ggsave(p1, file='/home/panxl/CRC2/kegg_up.pdf', width=12, height=8)
#信号通路图
##down
DEG.gene_ENSEMBL_down = as.character(res_down$Gene.stable.ID) #获得基因 gene ID
DEG.entrez_id_down = mapIds(x = org.Mm.eg.db,
                            keys = DEG.gene_ENSEMBL_down,
                            keytype = "ENSEMBL",
                            column = "ENTREZID")#
DEG.entrez_id_down = na.omit(DEG.entrez_id_down)#去除NA 817
kegg_down <- enrichKEGG(gene       = DEG.entrez_id_down,
                        keyType       = "kegg",
                        organism      = "mmu",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)


kegg_down_result=kegg_down@result
#write.csv(kegg_down_result,file = "/home/panxl/CRC2/mouse_DEG_kegg_down.txt")
library(ggplot2)
#install.packages('patchwork')
library(patchwork)
kegg <- read.csv("/home/panxl/CRC2/mouse_DEG_kegg_down.txt",stringsAsFactors=F)
kegg <- kegg[order(kegg$pvalue),]
top20 <- kegg[1:20,]
top20 <- data.frame(top20$Description,top20$Count ,top20$pvalue,top20$p.adjust)
colnames(top20) <- c("Description","count","pvalue","padj")
top20=top20%>%mutate(Description = fct_reorder(Description, -pvalue))
p2=ggplot(data=top20,aes(x=Description,y=-log10(pvalue),fill=padj))+#fill=padj fill颜色填充，使用连续值padj 
  geom_bar(stat="identity") + coord_flip()+#coord_flip()颠倒坐标轴
  labs(x="",y="-log10 (P_value)",title="KEGG Down Pathway")+
  scale_fill_gradient(low="#228B22",high="#228B22")+
  theme_classic()+
  theme(axis.title.x = element_text(size =22, 
                                    color = "black"),
        #face = "italic"),
        axis.title.y = element_text(size = 27,
                                    color = "black",
                                    #face = "italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        title = element_text(size = 22,
                             color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 17 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   #face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 19,  
                                   color = "black",
                                   vjust = 0.4))
#ggsave(p2, file='/home/panxl/CRC2/kegg_down.pdf', width=12, height=8)

######################4.2 Mouse kegg pathway挑出来的通路热图（10）####################
####mouse kegg up
mouse_kegg_up<- read.csv("/home/panxl/CRC2/mouse_DEG_kegg_up.txt",stringsAsFactors=F)
mouse_kegg_up$Description=gsub("Cytokine-cytokine receptor interaction","Cytokine interaction",mouse_kegg_up$Description)
mouse_kegg_up$Description=gsub("IL-17 signaling pathway","IL-17 signaling",mouse_kegg_up$Description)
mouse_kegg_up$Description=gsub("Cell adhesion molecules","Cell adhesion",mouse_kegg_up$Description)
mouse_kegg_up$Description=gsub("TNF signaling pathway","TNF signaling",mouse_kegg_up$Description)
mouse_kegg_up$Description=gsub("Wnt signaling pathway","Wnt signaling",mouse_kegg_up$Description)
mouse_kegg_up$Description=gsub("Calcium signaling pathway","Calcium signaling",mouse_kegg_up$Description)
mouse_kegg_up$Description=gsub("cAMP signaling pathway","cAMP signaling",mouse_kegg_up$Description)
mouse_kegg_up$Description=gsub("Bile secretion","Bile secretion",mouse_kegg_up$Description)
mouse_kegg_up$Description=gsub("Drug metabolism - cytochrome P450","Cytochrome P450",mouse_kegg_up$Description)
mouse_kegg_up$Description=gsub("Chemical carcinogenesis - DNA adducts","Chemical carcinogenesis",mouse_kegg_up$Description)

mouse_kegg_up20 <- mouse_kegg_up[order(mouse_kegg_up$p.adjust),][c(1:20),]
mouse_kegg_up20$Description
mouse_kegg_up5=mouse_kegg_up20[c(1,8,9,10,20),]
mouse_kegg_up5$Description
mouse_kegg_up5$type="Mouse Up"
####mouse kegg down
mouse_kegg_down <- read.csv("/home/panxl/CRC2/mouse_DEG_kegg_down.txt",stringsAsFactors=F)
mouse_kegg_down$Description=gsub("Cytokine-cytokine receptor interaction","Cytokine interaction",mouse_kegg_down$Description)
mouse_kegg_down$Description=gsub("IL-17 signaling pathway","IL-17 signaling",mouse_kegg_down$Description)
mouse_kegg_down$Description=gsub("Cell adhesion molecules","Cell adhesion",mouse_kegg_down$Description)
mouse_kegg_down$Description=gsub("TNF signaling pathway","TNF signaling",mouse_kegg_down$Description)
mouse_kegg_down$Description=gsub("Wnt signaling pathway","Wnt signaling",mouse_kegg_down$Description)
mouse_kegg_down$Description=gsub("Calcium signaling pathway","Calcium signaling",mouse_kegg_down$Description)
mouse_kegg_down$Description=gsub("cAMP signaling pathway","cAMP signaling",mouse_kegg_down$Description)
mouse_kegg_down$Description=gsub("Bile secretion","Bile secretion",mouse_kegg_down$Description)
mouse_kegg_down$Description=gsub("Drug metabolism - cytochrome P450","Cytochrome P450",mouse_kegg_down$Description)
mouse_kegg_down$Description=gsub("Chemical carcinogenesis - DNA adducts","Chemical carcinogenesis",mouse_kegg_down$Description)
mouse_kegg_down20 <- mouse_kegg_down[order(mouse_kegg_down$p.adjust),][c(1:20),]

mouse_kegg_down5=mouse_kegg_down20[c(2,3,5,6,15),]
mouse_kegg_down5$Description
mouse_kegg_down5$type="Mouse Down"
###mouse pathway
mouse_pathway=rbind(mouse_kegg_up5,mouse_kegg_down5)
up_pathway=mouse_kegg_up[mouse_kegg_up$Description%in%mouse_pathway$Description,][,c(3,7)]
up_pathway$Description
colnames(up_pathway)=c("Description", "Mouse Up")
down_pathway=mouse_kegg_down[mouse_kegg_down$Description%in%mouse_pathway$Description,][,c(3,7)]
colnames(down_pathway)=c("Description", "Mouse Down")
pathway=left_join(up_pathway,down_pathway,by = "Description")
rownames(pathway)=pathway$Description
pathway=pathway[,-1]
ha1 = HeatmapAnnotation(Group=factor(c("Mouse Up", "Mouse Down"),
                                     levels = c("Mouse Up", "Mouse Down")),
                        col = list(Group = c("Mouse Up"="#FF9999","Mouse Down"="#96A48B")),show_legend = F)
ht_list=Heatmap(-log10(pathway), cluster_rows = F,
                cluster_columns = F,
                show_column_names = T,
                heatmap_legend_param = list(
                  title = "-log10(padj)", at = c(0, 5,30),
                  labels = c("low","", "high")
                ),
                show_row_names = T,
                col=colorRamp2(c(0,5, 30), c( "#85C1E9", "white","#EF8A62" )),
                top_annotation = ha1,
                column_names_gp  =gpar(
                  col = "black",
                  fontsize = 10
                ),
                column_names_rot =45,
                row_names_gp = gpar(fontsize = 13),
                show_heatmap_legend = T)
draw(ht_list, heatmap_legend_side = "left")

# lgd = Legend(col_fun = colorRamp2(c(0,15, 30), c( "#67A9D1", "white","#EF8A62" )), title = "-log10(padj)", at = c(0, 30), 
#              labels = c("low", "high"))#85C1E9
# draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
########################################### 5.Human and mouse fpkm#######################################################################################
#################### 5.1 mouse fpkm#######################
# mouse_info <- read.csv("/home/panxl/CRC2/filter_mouse_info.csv",row.names = 1, check.names = F)
# summary_readcounts1 <- read.table("/home/panxl/CRC/all_mouse_data/summary_count.txt",sep = ",",header=TRUE,stringsAsFactors = F)
# rownames(summary_readcounts1)=summary_readcounts1[,1]
# summary_readcounts=summary_readcounts1[,colnames(summary_readcounts1)%in%mouse_info$SRR_id]
#write.csv(summary_readcounts,file = "/home/panxl/CRC2/mouse78_count.csv")
summary_readcounts=read.csv("/home/panxl/CRC2/mouse78_count.csv",row.names = 1)
#summary_readcounts=summary_readcounts[grepl("ENS",rownames(summary_readcounts)),]
library(GenomicFeatures)
## Calculate gene length for 
txdb <- makeTxDbFromGFF("//home/panxl/reference/mouse_reference/gencode.vM20.annotation.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(GenomicRanges::reduce(x)))})
names(exons_gene_lens) = gsub("\\..*","",names(exons_gene_lens))
exons_genelens=unlist(exons_gene_lens)
total_reads=apply(summary_readcounts,2,sum)
exons_gene_len=exons_genelens[match(rownames(summary_readcounts),names(exons_genelens))]
FPKM_summary=(10^9*t(t(summary_readcounts)/total_reads))/(exons_gene_len)# n26 t38 n3 t6
tumor=mouse_info[mouse_info$Group=="Tumor",]$SRR_id
FPKM_summary_tumor=FPKM_summary[,colnames(FPKM_summary)%in%tumor]#44
normal=mouse_info[mouse_info$Group=="Normal",]$SRR_id
FPKM_summary_normal=FPKM_summary[,colnames(FPKM_summary)%in%normal]#34
FPKM_summary1=cbind(FPKM_summary_normal,FPKM_summary_tumor)#34 44
#write.table(FPKM_summary1,"/home/panxl/CRC2/mouse_FPKM_summary78.txt", col.names=NA,sep="\t",quote=F) #normal29  tumor 44
#################### 5.2   human fpkm###################
human_info <- read.csv("/home/panxl/CRC2/filter_human_info.csv",row.names = 1, check.names = F)
summary_readcounts1 <- read.csv("/home/panxl/CRC/RNA_seq/wuda_Gene_exp.csv",row.names = 1)
summary_readcounts1$gene=rownames(summary_readcounts1)
summary_readcounts2 <- read.csv("/home/shimw/project/conservation/COAD_hg19_count.csv",check.names = F,row.names = 1)
rownames(summary_readcounts2)=gsub("\\..*","",rownames(summary_readcounts2))
summary_readcounts2$gene=rownames(summary_readcounts2)

summary_readcounts3=left_join(summary_readcounts1,summary_readcounts2,c("gene"="gene"))
summary_readcounts4=na.omit(summary_readcounts3)#473
rownames(summary_readcounts4)=summary_readcounts4$gene
summary_readcounts8=summary_readcounts4[,colnames(summary_readcounts4)%in%human_info$sample]
#write.csv(summary_readcounts8,file = "/home/panxl/CRC2/human439_count.csv")
summary_readcounts=read.csv("/home/panxl/CRC2/human439_count.csv",row.names = 1)
library(GenomicFeatures)
## Calculate gene length for human-92.
txdb <- makeTxDbFromGFF("/home/panxl/reference/human_reference/gencode.v37lift37.annotation.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(GenomicRanges::reduce(x)))})
names(exons_gene_lens) = gsub("\\..*","",names(exons_gene_lens))
exons_genelens=unlist(exons_gene_lens)  
total_reads=apply(summary_readcounts,2,sum)
exons_gene_len=exons_genelens[match(rownames(summary_readcounts),names(exons_genelens))]
FPKM_summary=(10^9*t(t(summary_readcounts)/total_reads))/(exons_gene_len)
colnames(FPKM_summary)=gsub("[.]","-",colnames(FPKM_summary))
tumor=human_info[human_info$Group=="Tumor",]$sample
FPKM_summary_tumor=FPKM_summary[,colnames(FPKM_summary)%in%tumor]#335
normal=human_info[human_info$Group=="Normal",]$sample
FPKM_summary_normal=FPKM_summary[,colnames(FPKM_summary)%in%normal]#104
FPKM_summary1=cbind(FPKM_summary_normal,FPKM_summary_tumor)#104 335
write.table(FPKM_summary1,"/home/panxl/CRC2/human_FPKM_summary.txt", col.names=NA,sep="\t",quote=F) #normal29  tumor 44
a=read.csv("/home/panxl/CRC2/human_FPKM_summary.txt",sep = "\t",check.names = F,row.names = 1)




#################### 5.3  CMS fpkm#######################
library(GenomicFeatures)
## Calculate gene length for human-92.
txdb <- makeTxDbFromGFF("/home/panxl/reference/human_reference/gencode.v37lift37.annotation.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(GenomicRanges::reduce(x)))})
names(exons_gene_lens) = gsub("\\..*","",names(exons_gene_lens))
exons_genelens=unlist(exons_gene_lens)
save(exons_genelens,file = "/home/panxl/CRC2/exons_genelens.RData")

df=c("CMS1","CMS2","CMS3","CMS4")
CMS_FPKM= function(x){
  CMS=summary_readcounts[,colnames(summary_readcounts)%in%cmstypes[cmstypes$type==x,]$barcode]
  total_reads=apply(CMS,2,sum)
  exons_gene_len=exons_genelens[match(rownames(CMS),names(exons_genelens))]
  FPKM_summary=(10^9*t(t(CMS)/total_reads))/(exons_gene_len)
  write.csv(FPKM_summary,paste("/home/panxl/CRC2/","fpkm_",x,".csv",sep=""))
}
CMS_FPKM("CMS1")
CMS_FPKM("CMS2")
CMS_FPKM("CMS3")
CMS_FPKM("CMS4")

########################################### 6.Wnt pathway heatmap###############################################################################
#mmu04310	mmu04310	Wnt signaling pathway
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(KEGGREST)
load("/home/panxl/CRC2/hom_kegg.RData")
Wnt_pathway=hom_kegg[hom_kegg$PATH=="04310",]
library(tidyverse)
library(pheatmap)

all_mouse_fpkm=read.csv("/home/panxl/CRC_Conservation/Mouse/FPKM_summary73.txt",sep="\t",row.names = 1)
Wnt_gene=all_mouse_fpkm[rownames(all_mouse_fpkm)%in%Wnt_pathway$Mouse.gene.stable.ID,]
Wnt_gene_df=rownames(Wnt_gene)
Wnt_gene_SYMBOL= mapIds(x = org.Mm.eg.db,
                        keys = Wnt_gene_df,
                        keytype = "ENSEMBL",
                        column ="SYMBOL")#133
Wnt_gene=cbind(Wnt_gene,Wnt_gene_SYMBOL)
rownames(Wnt_gene)=Wnt_gene$Wnt_gene_SYMBOL
Wnt_gene=Wnt_gene[,-c(74,75)]
annotation_col= data.frame(Group = factor(c(rep("Normal", 29), rep("Tumor", 44)), levels = c("Normal", "Tumor")))
row.names(annotation_col) <- colnames(Wnt_gene)
library(pheatmap)
 detach("package:ComplexHeatmap", unload=TRUE)
p=pheatmap(Wnt_gene, cluster_rows = F, show_rownames = T,
           cluster_cols = T,
           scale = "row",
           annotation_col = annotation_col,
           color =colorRampPalette(c("#87CEFA","white", "#CC2121"))(100))
p
kegg <- read.csv("/home/panxl/CRC_Conservation/Mouse/mouse_DEG_kegg_up.txt",stringsAsFactors=F)
Wnt_pathway=kegg[kegg$X=="mmu04310",]

Wnt_gene_entrez=data.frame("22418","22422","22402","12443","24117","228770","14283",
                           "22421","17393","93960","16842","56811","68010","93735","14160","27373","93897","77583",
                           "21414","170758","329252","22415","20671","93840","22420","22409","207742","12006","407821")
Wnt_gene_entrez=as.character(Wnt_gene_entrez)
Wnt_gene_ENSEMBL= mapIds(x = org.Mm.eg.db,
                         keys = Wnt_gene_entrez,
                         keytype = "ENTREZID",
                         column ="ENSEMBL")#27
Wnt_gene_ENSEMBL = na.omit(Wnt_gene_ENSEMBL)#去除NA 27
library(tidyverse)
library(pheatmap)

all_mouse_fpkm=read.csv("/home/panxl/CRC_Conservation/Mouse/FPKM_summary73.txt",sep="\t",row.names = 1)
Wnt_gene=all_mouse_fpkm[rownames(all_mouse_fpkm)%in%Wnt_gene_ENSEMBL,]
Wnt_gene_df=rownames(Wnt_gene)
Wnt_gene_SYMBOL= mapIds(x = org.Mm.eg.db,
                        keys = Wnt_gene_df,
                        keytype = "ENSEMBL",
                        column ="SYMBOL")#27
Wnt_gene=cbind(Wnt_gene,Wnt_gene_SYMBOL)
rownames(Wnt_gene)=Wnt_gene$Wnt_gene_SYMBOL
Wnt_gene=Wnt_gene[,-74]
annotation_col= data.frame(Group = factor(c(rep("Normal", 29), rep("Tumor", 44)), levels = c("Normal", "Tumor")))
row.names(annotation_col) <- colnames(Wnt_gene)
p=pheatmap(Wnt_gene, cluster_rows = F, show_rownames = T,
           cluster_cols = F,
           scale = "row",
           annotation_col = annotation_col,
           color =colorRampPalette(c("#87CEFA","white", "#CC2121"))(100))

######################################## 7.1.Fpkm correlation of normal and tumor #######################################################################################
#human   tumor467, native,41
#mouse   norma: the former29    tumor: the latter 44
#human_fpkm_normal 41
human_fpkm_summary=read.csv("/home/panxl/CRC2/human_FPKM_summary.txt",sep = "\t",row.names = 1)
human_fpkm_normal=human_fpkm_summary[,c(1:104)]
human_fpkm_normal$human_normal_rowmean=rowMeans(human_fpkm_normal)
human_fpkm_normal$human_normal_geneid=rownames(human_fpkm_normal)
human_fpkm_normal=human_fpkm_normal[,c(105,106)]
#human_fpkm_tumor 467
human_fpkm_tumor=human_fpkm_summary[,c(105:439)]
human_fpkm_tumor$human_tumor_rowmean=rowMeans(human_fpkm_tumor)
human_fpkm_tumor$human_tumor_geneid=rownames(human_fpkm_tumor)
human_fpkm_tumor=human_fpkm_tumor[,c(336,337)]
#mouse_fpkm_normal 29
mouse_fpkm_summary=read.csv("/home/panxl/CRC2/mouse_FPKM_summary78.txt",sep="\t",row.names=1)
mouse_fpkm_normal=mouse_fpkm_summary[,c(1:34)]
mouse_fpkm_normal$mouse_normal_rowmean=rowMeans(mouse_fpkm_normal)
mouse_fpkm_normal$mouse_normal_geneid=rownames(mouse_fpkm_normal)
mouse_fpkm_normal=mouse_fpkm_normal[,c(35,36)]
#mouse_fpkm_tumor44
mouse_fpkm_tumor=mouse_fpkm_summary[,c(35:78)]
mouse_fpkm_tumor$mouse_tumor_rowmean=rowMeans(mouse_fpkm_tumor)
mouse_fpkm_tumor$mouse_tumor_geneid=rownames(mouse_fpkm_tumor)
mouse_fpkm_tumor=mouse_fpkm_tumor[,c(45,46)]
#mouse human homologues
HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]#16550

#normal
##mouse
mouse_normal_HOM_fpkm =mouse_fpkm_normal[mouse_fpkm_normal$mouse_normal_geneid%in%HOM_one2one$Mouse.gene.stable.ID,]#16248
HOM_mouse_normal=HOM_one2one[HOM_one2one$Mouse.gene.stable.ID%in%mouse_fpkm_normal$mouse_normal_geneid,]#16248

##H0M_data
human_normal=human_fpkm_normal[human_fpkm_normal$human_normal_geneid%in%HOM_mouse_normal$Gene.stable.ID,]#16110
HOM_normal=HOM_one2one[HOM_one2one$Gene.stable.ID%in%human_normal$human_normal_geneid,]#16141
mouse_normal=mouse_fpkm_normal[mouse_fpkm_normal$mouse_normal_geneid%in%HOM_normal$Mouse.gene.stable.ID,]#16110

#tumor

##mouse
mouse_tumor_HOM_fpkm =mouse_fpkm_tumor[mouse_fpkm_tumor$mouse_tumor_geneid%in%HOM_one2one$Mouse.gene.stable.ID,]#16248
HOM_mouse_tumor=HOM_one2one[HOM_one2one$Mouse.gene.stable.ID%in%mouse_fpkm_tumor$mouse_tumor_geneid,]#16248

##H0M_data
human_tumor=human_fpkm_tumor[human_fpkm_tumor$human_tumor_geneid%in%HOM_mouse_tumor$Gene.stable.ID,]#16110
HOM_tumor=HOM_one2one[HOM_one2one$Gene.stable.ID%in%human_tumor$human_tumor_geneid,]#16110
mouse_tumor=mouse_fpkm_tumor[mouse_fpkm_tumor$mouse_tumor_geneid%in%HOM_tumor$Mouse.gene.stable.ID,]# 16110 16141

#fpkm_data
##normal
colnames(human_normal)=c("human_normal_rowmean","Gene.stable.ID")
colnames(mouse_normal)=c("mouse_normal_rowmean","Mouse.gene.stable.ID")
fpkm_normal1=merge(human_normal,HOM_normal,by="Gene.stable.ID")
fpkm_normal=merge(fpkm_normal1,mouse_normal,by="Mouse.gene.stable.ID")
##tumor
colnames(human_tumor)=c("human_tumor_rowmean","Gene.stable.ID")
colnames(mouse_tumor)=c("mouse_tumor_rowmean","Mouse.gene.stable.ID")
fpkm_tumor1=merge(human_tumor,HOM_tumor,by="Gene.stable.ID")
fpkm_tumor=merge(fpkm_tumor1,mouse_tumor,by="Mouse.gene.stable.ID")


##画图
library(ggplot2)
library("ggsci")
library(dplyr)  
p1=ggplot(fpkm_normal, aes(x=log10(human_normal_rowmean+1),y=log10(mouse_normal_rowmean+1))) +
  geom_point(size=1) +
  #geom_abline(slope=1,intercept=0,color="red",size=1)+
  scale_color_uchicago() +
  annotate("text", x = 1.5, y = 3.9, label = "r=0.70, p-value <2e-16",
           color="#350E20FF",size = 9 )+
  #geom_smooth(method=lm,color="red")+
  labs(x="Human log10(FPKM+1)",y="Mouse log10(FPKM+1)",title="Normal FPKM Correlation")+
  theme_classic(base_line_size = 1) +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24, color="black"),
        title = element_text(size = 30,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"))
# axis.ticks.length.y = unit(2,"mm")
#)#+geom_text(data=ann,aes(label=Gene.name))
#correlation
cor(fpkm_normal$human_normal_rowmean, fpkm_normal$mouse_normal_rowmean) #0.7  2e-16      以前的值0.6989655
cor.test(fpkm_normal$human_normal_rowmean, fpkm_normal$mouse_normal_rowmean,method = "kendall", alternative = "greater")#normal0.6473  p-value < 2.2e-16
p2=ggplot(fpkm_tumor, aes(x=log10(human_tumor_rowmean+1),y=log10(mouse_tumor_rowmean+1))) +
  geom_point(size=1) +
  annotate("text",  x = 1.5, y = 3.9, label = "r=0.57, p-value <2e-16",
           color="#350E20FF",size = 9 )+
 # geom_smooth(method=lm,color="red")+
  scale_color_uchicago() +
  labs(x="Human log10(FPKM+1)",y="Mouse log10(FPKM+1)",title="Tumor FPKM Correlation")+
  theme_classic(base_line_size = 1) +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24, color="black"),
        title = element_text(size = 30,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        # axis.ticks.length.y = unit(2,"mm")
  )#+geom_text(data=ann,aes(label=Gene.name))
#correlation
cor(fpkm_tumor$human_tumor_rowmean, fpkm_tumor$mouse_tumor_rowmean) # 0.565   2.2e-16            以前的值 0.5971579
cor.test(fpkm_tumor$human_tumor_rowmean, fpkm_tumor$mouse_tumor_rowmean,method = "kendall", alternative = "greater")#  tumor0.6491  p-value < 2.2e-16
library(patchwork)
p1+p2
ggsave(p1, file='/home/panxl/CRC2/fpkm_normal.pdf', width=7, height=7) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
ggsave(p2, file='/home/panxl/CRC2/fpkm_tumor.pdf', width=7, height=7) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
######################################## 7.2.CMS Fpkm correlation  #######################################################################################
#mouse human homologues
HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]#16550

#CMS1 59
df=c("CMS1","CMS2","CMS3","CMS4")
CMS_meanfpkm=lapply(df, function(x){
  CMS1_fpkm_summary=read.csv(paste("/home/panxl/CRC2/fpkm_",x,".csv",sep=""),row.names = 1)
  CMS1_fpkm_summary$CMS1_normal_rowmean=rowMeans(CMS1_fpkm_summary)
  CMS1_fpkm_summary$CMS1_tumor_geneid=rownames(CMS1_fpkm_summary)
  CMS1_fpkm_summaryr=CMS1_fpkm_summary[,c(length(CMS1_fpkm_summary)-1,length(CMS1_fpkm_summary))]
  CMS1_HOM_fpkm =na.omit(left_join(CMS1_fpkm_summaryr,HOM_one2one,c("CMS1_tumor_geneid"="Gene.stable.ID")))[,-4]
data=data.frame( fpkm=CMS1_HOM_fpkm$CMS1_normal_rowmean, 
                 geneid= CMS1_HOM_fpkm$CMS1_tumor_geneid, 
                 Mouse.gene.stable.ID=CMS1_HOM_fpkm$Mouse.gene.stable.ID)
  return(data )
  })%>%dplyr::bind_cols()
CMS_meanfpkm1=CMS_meanfpkm[,c(1,4,7,10,11,12)]
colnames(CMS_meanfpkm1)=c("CMS1","CMS2","CMS3","CMS4","Gene.stable.ID","Mouse.gene.stable.ID")
save(CMS_meanfpkm1,file = "/home/panxl/CRC2/CMS_meanfpkm.RData")
#fpkm_data
load("/home/panxl/CRC2/CMS_meanfpkm.RData")
#mouse_fpkm_tumor44
mouse_fpkm_summary=read.csv("/home/panxl/CRC2/mouse_FPKM_summary78.txt",sep="\t",row.names=1)
mouse_fpkm_tumor=mouse_fpkm_summary[,c(35:78)]
# mouse_info=read.csv("/home/panxl/CRC2/filter_mouse_info.csv")
# mouse_tumor_info=mouse_info[mouse_info$Group=="Tumor",]$SRR_id
# intersect(mouse_tumor_info,colnames(mouse_fpkm_tumor))
mouse_fpkm_tumor$mouse_tumor_rowmean=rowMeans(mouse_fpkm_tumor)
mouse_fpkm_tumor$mouse_tumor_geneid=rownames(mouse_fpkm_tumor)
mouse_fpkm_tumor=mouse_fpkm_tumor[,c(45,46)]
mouse_tumor_HOM_fpkm =na.omit(left_join(mouse_fpkm_tumor,HOM_one2one,c("mouse_tumor_geneid"="Mouse.gene.stable.ID")))[,-4]
data=left_join(CMS_meanfpkm1,mouse_tumor_HOM_fpkm,by = "Gene.stable.ID")

##画图
library(ggplot2)
library("ggsci")
library(dplyr)  
p2=ggplot(data, aes(x=log10(CMS1+1),y=log10(mouse_tumor_rowmean+1))) +
  geom_point(size=1) +
  annotate("text", x = 1.5, y = 3.9, label = "r=0.565, p-value <2e-16",
           color="#350E20FF",size = 9 )+
  geom_smooth(method=lm,color="red")+
  scale_color_uchicago() +
  labs(x="Human log10(FPKM+1)",y="Mouse log10(FPKM+1)",title="Tumor FPKM Correlation")+
  theme_classic(base_line_size = 1) +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24, color="black"),
        title = element_text(size = 30,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        # axis.ticks.length.y = unit(2,"mm")
  )#+geom_text(data=ann,aes(label=Gene.name))
#correlation
data=na.omit(data)
cor(data$CMS4,data$mouse_tumor_rowmean) # 0.563   0.55 0.622  0.557        以前的值 0.5971579
cor.test(fpkm_tumor$human_tumor_rowmean, fpkm_tumor$mouse_tumor_rowmean,method = "kendall", alternative = "greater")#  tumor0.6491  p-value < 2.2e-16
library(patchwork)
p1+p2
ggsave(p1, file='/home/panxl/CRC2/fpkm_normal.pdf', width=7, height=7) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
ggsave(p2, file='/home/panxl/CRC2/fpkm_tumor.pdf', width=7, height=7) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径



####################################### 8.1 CRC 同源基因logfc在cluster分布################################################
################# 8.1  cluster  human HOM logfc protein coding 
##human 未分型HOM logfc
HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]#16550
#mouse
library(tidyverse)
all_mouse_DEG=read.csv("/home/panxl/CRC2/mouse_DEG.csv")#22413
mouse_logFC1 =all_mouse_DEG[all_mouse_DEG$X%in%HOM_one2one$Mouse.gene.stable.ID & all_mouse_DEG$adj.P.Val<=0.05,][,c(1,2)]#10980
a=HOM_one2one[HOM_one2one$Mouse.gene.stable.ID%in%mouse_logFC1$X,]#10980
human_data= read.csv("/home/panxl/CRC2/human_DEG.csv")#27542
human_logFC1  =human_data[human_data $X%in%HOM_one2one$Gene.stable.ID &human_data$adj.P.Val<=0.05,][,c(1,2)]#12526
b=HOM_one2one[HOM_one2one$Gene.stable.ID%in%human_logFC1$X,]#12526
#H0M_logFC1_data
c=unique(rbind(a,b))#14233
hh=na.omit(left_join(human_data,c,c("X"="Gene.stable.ID" )))#14027
mm=na.omit(left_join(all_mouse_DEG,c,c("X"="Mouse.gene.stable.ID")))#13376
hom=na.omit(left_join(hh,mm,c("X"="Gene.stable.ID")))#13170
H0M_logFC1=hom[,c(1,2,length(hh)-1,length(hh)+2)]#12871
colnames(H0M_logFC1)=c("Gene.stable.ID","Human_logFC1","Mouse.gene.stable.ID","Mouse_logFC1")
H0M_logFC1$cluster=ifelse(abs(H0M_logFC1$Mouse_logFC1)<=1&abs(H0M_logFC1$Human_logFC1)<=1,"Cluster 5: Not-significant in both species",case_when(
  H0M_logFC1$Mouse_logFC1>=0&H0M_logFC1$Human_logFC1>=0 ~ "Cluster 1: Up-regulated in both species",
  H0M_logFC1$Mouse_logFC1>=0&H0M_logFC1$Human_logFC1<=0 ~ "Cluster 2: Up-regulated in mouse reverse in human",
  H0M_logFC1$Mouse_logFC1<=0&H0M_logFC1$Human_logFC1<=0 ~ "Cluster 3: Down-regulated in both species",
  H0M_logFC1$Mouse_logFC1<=0&H0M_logFC1$Human_logFC1>=0 ~ "Cluster 4: Up-regulated in human reverse in mouse"))
# library(biomaRt)
# library(ggsci)
# library(ggrepel)
# library(reprex)
# mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "ensembl.org")
# mart <- useDataset("hsapiens_gene_ensembl", mart)
# annotLookup <- getBM(mart = mart,
#                      attributes = c("external_gene_name","ensembl_gene_id", "gene_biotype"),
#                      filters = "ensembl_gene_id", values = H0M_logFC1$Gene.stable.ID,
#                      uniqueRows=TRUE)
# 
# names(annotLookup)[2] = "Gene.stable.ID"
# H0M_logFC1 = dplyr::left_join(H0M_logFC1, annotLookup)
#write.csv(H0M_logFC1,file = "/home/panxl/CRC2/expression_conservation_and_divergence_p<0.05.csv",row.names = F)

all_mouse_DEG=read.csv("/home/panxl/CRC2/mouse_DEG.csv")#22413
a=HOM_one2one[HOM_one2one$Mouse.gene.stable.ID%in%all_mouse_DEG$X,]#13947
human_data= read.csv("/home/panxl/CRC2/human_DEG.csv")
b=HOM_one2one[HOM_one2one$Gene.stable.ID%in%human_data$X,]#14836
#H0M_logFC1_data
c=unique(rbind(a,b))#15123
hh=na.omit(left_join(human_data,c,c("X"="Gene.stable.ID" )))#14836
mm=na.omit(left_join(all_mouse_DEG,c,c("X"="Mouse.gene.stable.ID")))#13947
hom=na.omit(left_join(hh,mm,c("X"="Gene.stable.ID")))#13660
cluster5=hom[!(hom$X%in%H0M_logFC1$Gene.stable.ID),c(1,2,length(hh)-1,length(hh)+2)]#490
colnames(cluster5)=c("Gene.stable.ID","Human_logFC1","Mouse.gene.stable.ID","Mouse_logFC1")
cluster5$cluster="Cluster 5: Not-significant in both species"
# mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "ensembl.org")
# mart <- useDataset("hsapiens_gene_ensembl", mart)
# annotLookup <- getBM(mart = mart,
#                      attributes = c("external_gene_name","ensembl_gene_id", "gene_biotype"),
#                      filters = "ensembl_gene_id", values = cluster5$Gene.stable.ID,
#                      uniqueRows=TRUE)
# 
# names(annotLookup)[2] = "Gene.stable.ID"
# cluster5 = dplyr::left_join(cluster5, annotLookup)#6429

#write.csv(cluster5,file = "/home/panxl/CRC2/cluster5_p>0.05.csv",row.names = F)
H0M_logFC1= read.csv("/home/panxl/CRC2/expression_conservation_and_divergence_p<0.05.csv")
all_hom_cluater=rbind(H0M_logFC1,cluster5)#13660
#write.csv(all_hom_cluater,file = "/home/panxl/CRC2/all_cluater",row.names = F)
########################################8.2.1  CMS分类  ##############################
library(CMScaller)
human_res=read.csv("/home/panxl/CRC2/human439_count.csv",check.names = F,row.names = 1)
human_info=read.csv("/home/panxl/CRC2/filter_human_info.csv",check.names = F,row.names = 1)
human_info_tumor=human_info[human_info$Group=="Tumor",]$sample#335
tumor_res=human_res[,colnames(human_res)%in%human_info_tumor]#335
tumor_cms <- CMScaller(tumor_res, RNAseq=TRUE, doPlot=TRUE, rowNames = "ensg", seed = 1)
human_cmstypes <- data.frame("barcode"=row.names(tumor_cms), "type"=tumor_cms$prediction)
human_cmstypes=na.omit(human_cmstypes)#305
#write.csv(human_cmstypes,file = "/home/panxl/CRC2/human_cmstypes.csv",row.names = F)
unique
cmstypes=read.csv("/home/panxl/CRC2/human_cmstypes.csv")

unique(cmstypes$type)
cmstypes$type <- as.vector(cmstypes$type)
cmstypes$type[is.na(cmstypes$type)] = "None"
cmstypes$type = factor(cmstypes$type)
table(cmstypes$type)
kk=as.data.frame(table(cmstypes$type))
names(kk) = c("type","number")
p1<-ggplot(data=kk, aes(x=type, y=number, fill = type)) +
  geom_bar(stat="identity", width=0.5)+
  labs(x="",y="",title="TCGA Human CMS Type")+
  theme_classic(base_line_size = 1) +
  theme(title = element_text(size = 22),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 20 ,# 修改X轴上字体大小，
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 20,  
                                   vjust = 0.4))

############################# 8.3.1 lim CMS DEG###########################################
library(limma)
human_data <- read.csv("/home/shimw/project/conservation/filter_human_data.csv",row.names = 1, check.names = F)
human_data
human_info <- read.csv("/home/panxl/CRC2/filter_human_info.csv",row.names = 1, check.names = F)
cmstypes=read.csv("/home/panxl/CRC2/human_cmstypes.csv")
human_info[human_info$Group=="Normal",]$sample
cmstypes[cmstypes$type=="CMS1",]$barcode

CMS_DEG=function(x,file){
  CMS_data=human_data[,colnames(human_data)%in%c(cmstypes[cmstypes$type==x,]$barcode,human_info[human_info$Group=="Normal",]$sample)]
  CMS_info=human_info[human_info$sample%in%c(cmstypes[cmstypes$type==x,]$barcode,human_info[human_info$Group=="Normal",]$sample),]
  fit=lmFit(CMS_data,design=model.matrix(~CMS_info$Group))  
  fit=eBayes(fit) 
  options(digits = 4) 
  topTable(fit,coef=2,adjust='BH') 
  deg=topTable(fit,coef=2,adjust='BH',number = Inf)
  write.csv(deg,paste("/home/panxl/CRC2/",file,"_DEG.csv",sep=""),quote = F)
}

CMS_DEG("CMS4", "CMS4")
############################ 8.3.2 lim CMS 同源基因logfc在cluster分布####################
HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]#16550
#mouse
library(tidyverse)
all_mouse_DEG=read.csv("/home/panxl/CRC2/mouse_DEG.csv")
mouse_logFC1 =all_mouse_DEG[all_mouse_DEG$X%in%HOM_one2one$Mouse.gene.stable.ID & all_mouse_DEG$adj.P.Val<=0.05,][,c(1,2)]
a=HOM_one2one[HOM_one2one$Mouse.gene.stable.ID%in%mouse_logFC1$X,]#10980
df=c("CMS1","CMS2","CMS3","CMS4")#purrr::map
result1=lapply(df,function(x){
 y=read.csv(paste("/home/panxl/CRC2/",x,"_DEG.csv",sep=""))
  cmstype_logFC1  =y[y$X%in%HOM_one2one$Gene.stable.ID&y$adj.P.Val<=0.05,]#9667
  b=HOM_one2one[HOM_one2one$Gene.stable.ID%in%cmstype_logFC1$X,]#9667
  #H0M_logFC1_data
  c=unique(rbind(a,b))#
  hh=na.omit(left_join(y,c,c("X"="Gene.stable.ID" )))#32831   12078
  mm=na.omit(left_join(all_mouse_DEG,c,c("X"="Mouse.gene.stable.ID")))#28255 11902
  hom=na.omit(left_join(hh,mm,c("X"="Gene.stable.ID")))#12078  11801
  hom$type=x
  H0M_logFC1=hom[,c(1,2,length(hom),length(hh)-1,length(hh)+2)]
  colnames(H0M_logFC1)=c("Gene.stable.ID","Human_logFC1","type","Mouse.gene.stable.ID","Mouse_logFC1")#11801
  H0M_logFC1$cluster=ifelse(abs(H0M_logFC1$Mouse_logFC1)<=1&abs(H0M_logFC1$Human_logFC1)<=1,"Cluster 5",case_when(
    H0M_logFC1$Mouse_logFC1>=0&H0M_logFC1$Human_logFC1>=0 ~ "Cluster 1",
    H0M_logFC1$Mouse_logFC1>=0&H0M_logFC1$Human_logFC1<=0 ~ "Cluster 2",
    H0M_logFC1$Mouse_logFC1<=0&H0M_logFC1$Human_logFC1<=0 ~ "Cluster 3",
    H0M_logFC1$Mouse_logFC1<=0&H0M_logFC1$Human_logFC1>=0 ~ "Cluster 4"))
  re_df <- data.frame("type"=H0M_logFC1$type,"Gene.stable.ID"=H0M_logFC1$Gene.stable.ID,"Human_logFC1"=H0M_logFC1$Human_logFC1,
                      "Mouse.gene.stable.ID"=H0M_logFC1$Mouse.gene.stable.ID,"Mouse_logFC1"=H0M_logFC1$Mouse_logFC1,
                      "cluster"=H0M_logFC1$cluster)
  return(re_df)
})%>%dplyr::bind_rows()
#write.csv(result1,file = "/home/panxl/CRC2/lim_CMS_expression_conservation_and_divergence_p<0.05.csv",row.names = F)#64912
all_mouse_DEG=read.csv("/home/panxl/CRC2/mouse_DEG.csv")
a=HOM_one2one[HOM_one2one$Mouse.gene.stable.ID%in%all_mouse_DEG$X,]
mouse_logFC1 =all_mouse_DEG[all_mouse_DEG$X%in%HOM_one2one$Mouse.gene.stable.ID & all_mouse_DEG$adj.P.Val<=0.05,][,c(1,2)]#10669
a1=HOM_one2one[HOM_one2one$Mouse.gene.stable.ID%in%mouse_logFC1$X,]
result2=lapply(df,function(x){
  y=read.csv(paste("/home/panxl/CRC2/",x,"_DEG.csv",sep=""))
  cmstype_logFC1  =y[y$X%in%HOM_one2one$Gene.stable.ID&y$adj.P.Val<=0.05,]#9667
  b1=HOM_one2one[HOM_one2one$Gene.stable.ID%in%cmstype_logFC1$X,]#9667
  #H0M_logFC1_data
  c1=unique(rbind(a1,b1))#
  hh1=na.omit(left_join(y,c1,c("X"="Gene.stable.ID" )))#32831   12078
  mm1=na.omit(left_join(all_mouse_DEG,c1,c("X"="Mouse.gene.stable.ID")))#28255 11902
  hom1=na.omit(left_join(hh1,mm1,c("X"="Gene.stable.ID")))#12078  11801
  hom1$type=x
  H0M_logFC1=hom1[,c(1,2,length(hom),length(hh)-1,length(hh)+2)]
  b=HOM_one2one[HOM_one2one$Gene.stable.ID%in%y$X,]#9667
  #H0M_logFC1_data
  c=unique(rbind(a,b))#15497
  hh=na.omit(left_join(y,c,c("X"="Gene.stable.ID" )))#32831   12078
  mm=na.omit(left_join(all_mouse_DEG,c,c("X"="Mouse.gene.stable.ID")))#28255 11902
  hom=na.omit(left_join(hh,mm,c("X"="Gene.stable.ID")))#12078  11801
  hom$type=x
  cluster5=hom[!(hom$X%in%H0M_logFC1$X),c(1,2,length(hom),length(hh)-1,length(hh)+2)]#865
  colnames(cluster5)=c("Gene.stable.ID","Human_logFC1","type","Mouse.gene.stable.ID","Mouse_logFC1")
  re_df <- data.frame("type"=cluster5$type,"Gene.stable.ID"=cluster5$Gene.stable.ID,"Human_logFC1"=cluster5$Human_logFC1,
                      "Mouse.gene.stable.ID"=cluster5$Mouse.gene.stable.ID,"Mouse_logFC1"=cluster5$Mouse_logFC1
  )
  return(re_df)
})%>%dplyr::bind_rows()#2402
result2$cluster="Cluster 5"

#write.csv(result2,file = "/home/panxl/CRC2/lim_CMS_cluster5_p>0.05.csv",row.names = F)#  6169 38158
all_cms_hom_cluater=rbind(result1,result2)#68532
#write.csv(all_cms_hom_cluater,file = "/home/panxl/CRC2/lim_CMS_all_cluater.csv",row.names = F)#72577
###CMS Cluster number 
df=c("CMS1","CMS2","CMS3","CMS4")
cms_DEG=read.csv("/home/panxl/CRC2/lim_CMS_all_cluater.csv")#72577
lim_CMS_cluster_number=lapply(df,function(x){
  ll=cms_DEG[cms_DEG$type==x,]
  oo=as.data.frame(table(ll$cluster))
  names(oo)=c("cluster","number")
  oo$type=x
  return(oo)
})%>%dplyr::bind_rows()
save(lim_CMS_cluster_number,file = "/home/panxl/CRC2/lim_CMS_cluster_number.Rdata")








############################# 8.4.1 dds CMS DEG###########################################

library(DESeq2)
human_count=read.csv("/home/panxl/CRC2/human439_count.csv",check.names = F,row.names = 1)
human_info=read.csv("/home/panxl/CRC2/filter_human_info.csv",check.names = F,row.names = 1)
human_info_tumor=human_info[human_info$Group=="Tumor",]$sample#335
human_info_normal=human_info[human_info$Group=="Normal",]$sample#104
tumor_count=human_count[,colnames(human_count)%in%human_info_tumor]#335
normal_count=human_count[,colnames(human_count)%in%human_info_normal]#104

deseq_func <- function(cmstypes){
  countdata <- cbind(tumor_count[,cmstypes1$barcode[cmstypes1$type==cmstypes]], normal_count)
  coldata = data.frame(row.names = colnames(countdata),
                       "type" = c(rep(cmstypes,length(colnames(countdata)) - length(colnames(normal_count))),
                                  rep("normal",length(colnames(normal_count))))
  )
  dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~DataSet+type)
  dds_obj <- DESeq(dds)
  res <- results(dds_obj, contrast=c("type", cmstypes, "normal"))
  res <- as.data.frame(res)
  res <- res[order(res$padj),]
  res <- tibble::rownames_to_column(res, var = "gene")
  return(tibble::tibble("type"=cmstypes, "res"=list(res)))
}
result <- purrr::map(as.vector(unique(cmstypes$type)), deseq_func)
result2 <- dplyr::bind_rows(result)

for (cmstype in as.vector(unique(cmstypes$type))){
  write_csv(result2$res[result2$type==cmstype][[1]],
            paste("/home/panxl/CRC2/human_", cmstype, ".csv", sep=""))
}
library(org.Hs.eg.db)
cms1=read.csv('/home/panxl/CRC2/human_CMS1.csv')
cms2=read.csv('/home/panxl/CRC2/human_CMS2.csv')
cms3=read.csv('/home/panxl/CRC2/human_CMS3.csv')
cms4=read.csv('/home/panxl/CRC2/human_CMS4.csv')
data1=rbind(cms1,cms2)
data2=rbind(data1,cms3)
data=rbind(data2,cms4)
write.csv(data,file = "/home/panxl/CRC2/human_CMS_DEG.csv",row.names = F)
gene_df = bitr(data$gene,
                         fromType="ENSEMBL",
                         toType="SYMBOL",
                           OrgDb = "org.Hs.eg.db", drop = TRUE)
DATA <- left_join(data,gene_df,c("gene"="ENSEMBL"))#176
DATA =na.omit(DATA )
#write.csv(DATA,file = "/home/panxl/CRC2/human_CMS_SYMBOL_DEG.csv",row.names = F)
human_CMS_DEG=read.csv("/home/panxl/CRC2/human_CMS_DEG.csv" )
#AXIN2 Ensembl:ENSG00000168646
#HNF4A ENSG00000101076





############################ 8.4.2 dds CMS 同源基因logfc在cluster分布####################
HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]#16550
#mouse
library(tidyverse)
all_mouse_DEG=read.csv("/home/panxl/CRC2/mouse_DEG.csv")
mouse_logFC1 =all_mouse_DEG[all_mouse_DEG$X%in%HOM_one2one$Mouse.gene.stable.ID & all_mouse_DEG$adj.P.Val<=0.05,][,c(1,2)]
a=HOM_one2one[HOM_one2one$Mouse.gene.stable.ID%in%mouse_logFC1$X,]#10980
df=c("CMS1","CMS2","CMS3","CMS4")#purrr::map
# cms4=read.csv('/home/panxl/CRC2/human_CMS4.csv')
# cms4$type="CMS4"
# write.csv(cms4,file = "/home/panxl/CRC2/human_CMS4.csv")

result1=lapply(df,function(x){
  x=read.csv(paste("/home/panxl/CRC2/human_",x,".csv",sep=""),row.names=1)
  cmstype_logFC1  =x[x$gene%in%HOM_one2one$Gene.stable.ID&x$padj<=0.05,]#9667
  b=HOM_one2one[HOM_one2one$Gene.stable.ID%in%cmstype_logFC1$gene,]#9667
  #H0M_logFC1_data
  c=unique(rbind(a,b))#
  hh=na.omit(left_join(x,c,c("gene"="Gene.stable.ID" )))#32831   12078
  mm=na.omit(left_join(all_mouse_DEG,c,c("X"="Mouse.gene.stable.ID")))#28255 11902
  hom=na.omit(left_join(hh,mm,c("gene"="Gene.stable.ID")))#12078  11801
  H0M_logFC1=hom[,c(1,3,8,length(hh)-1,length(hh)+2)]
  colnames(H0M_logFC1)=c("Gene.stable.ID","Human_logFC1","type","Mouse.gene.stable.ID","Mouse_logFC1")#11801
  H0M_logFC1$cluster=ifelse(abs(H0M_logFC1$Mouse_logFC1)<=1&abs(H0M_logFC1$Human_logFC1)<=1,"Cluster 5: Not-significant in both species",case_when(
    H0M_logFC1$Mouse_logFC1>=0&H0M_logFC1$Human_logFC1>=0 ~ "Cluster 1: Up-regulated in both species",
    H0M_logFC1$Mouse_logFC1>=0&H0M_logFC1$Human_logFC1<=0 ~ "Cluster 2: Up-regulated in mouse reverse in human",
    H0M_logFC1$Mouse_logFC1<=0&H0M_logFC1$Human_logFC1<=0 ~ "Cluster 3: Down-regulated in both species",
    H0M_logFC1$Mouse_logFC1<=0&H0M_logFC1$Human_logFC1>=0 ~ "Cluster 4: Up-regulated in human reverse in mouse"))
  re_df <- data.frame("type"=H0M_logFC1$type,"Gene.stable.ID"=H0M_logFC1$Gene.stable.ID,"Human_logFC1"=H0M_logFC1$Human_logFC1,
                      "Mouse.gene.stable.ID"=H0M_logFC1$Mouse.gene.stable.ID,"Mouse_logFC1"=H0M_logFC1$Mouse_logFC1,
                      "cluster"=H0M_logFC1$cluster)
  return(re_df)
})%>%dplyr::bind_rows()
#write.csv(result1,file = "/home/panxl/CRC2/CMS_expression_conservation_and_divergence_p<0.05.csv",row.names = F)#64912

all_mouse_DEG=read.csv("/home/panxl/CRC2/mouse_DEG.csv")
a=HOM_one2one[HOM_one2one$Mouse.gene.stable.ID%in%all_mouse_DEG$X,]
mouse_logFC1 =all_mouse_DEG[all_mouse_DEG$X%in%HOM_one2one$Mouse.gene.stable.ID & all_mouse_DEG$adj.P.Val<=0.05,][,c(1,2)]#10669
a1=HOM_one2one[HOM_one2one$Mouse.gene.stable.ID%in%mouse_logFC1$X,]
result2=lapply(df,function(x){
  x=read.csv(paste("/home/panxl/CRC2/human_",x,".csv",sep=""),row.names = 1)
  cmstype_logFC1  =x[x$gene%in%HOM_one2one$Gene.stable.ID&x$padj<=0.05,]#9667
  b1=HOM_one2one[HOM_one2one$Gene.stable.ID%in%cmstype_logFC1$gene,]#9667
  #H0M_logFC1_data
  c1=unique(rbind(a1,b1))#
  hh1=na.omit(left_join(x,c1,c("gene"="Gene.stable.ID" )))#32831   12078
  mm1=na.omit(left_join(all_mouse_DEG,c1,c("X"="Mouse.gene.stable.ID")))#28255 11902
  hom1=na.omit(left_join(hh1,mm1,c("gene"="Gene.stable.ID")))#12078  11801
  H0M_logFC1=hom1[,c(1,3,8,length(hh1)-1,length(hh1)+2)]
  b=HOM_one2one[HOM_one2one$Gene.stable.ID%in%x$gene,]#9667
  #H0M_logFC1_data
  c=unique(rbind(a,b))#15497
  hh=na.omit(left_join(x,c,c("gene"="Gene.stable.ID" )))#32831   12078
  mm=na.omit(left_join(all_mouse_DEG,c,c("X"="Mouse.gene.stable.ID")))#28255 11902
  hom=na.omit(left_join(hh,mm,c("gene"="Gene.stable.ID")))#12078  11801
  cluster5=hom[!(hom$gene%in%H0M_logFC1$gene),c(1,3,8,length(hh)-1,length(hh)+2)]#865
  colnames(cluster5)=c("Gene.stable.ID","Human_logFC1","type","Mouse.gene.stable.ID","Mouse_logFC1")
  re_df <- data.frame("type"=cluster5$type,"Gene.stable.ID"=cluster5$Gene.stable.ID,"Human_logFC1"=cluster5$Human_logFC1,
                      "Mouse.gene.stable.ID"=cluster5$Mouse.gene.stable.ID,"Mouse_logFC1"=cluster5$Mouse_logFC1
  )
  return(re_df)
})%>%dplyr::bind_rows()#5199
result2$cluster="Cluster 5: Not-significant in both species"

#write.csv(result2,file = "/home/panxl/CRC2/CMS_cluster5_p>0.05.csv",row.names = F)#  6169 38158
all_cms_hom_cluater=rbind(result1,result2)#68532
#write.csv(all_cms_hom_cluater,file = "/home/panxl/CRC2/CMS_all_cluater.csv",row.names = F)#72577
###CMS Cluster number 
df=c("CMS1","CMS2","CMS3","CMS4")
cms_DEG=read.csv("/home/panxl/CRC2/CMS_all_cluater.csv")#72577
CMS_cluster_number=lapply(df,function(x){
  ll=cms_DEG[cms_DEG$type==x,]
  oo=as.data.frame(table(ll$cluster))
  names(oo)=c("cluster","number")
  oo$type=x
  return(oo)
})%>%dplyr::bind_rows()
save(CMS_cluster_number,file = "/home/panxl/CRC2/CMS_cluster_number.Rdata")
##################################9.1 human Cluster 散点图################
#画图
library(ggplot2)
library("ggsci")
library(dplyr)
library(ggrepel)
H0M_logFC1= read.csv("/home/panxl/CRC2/all_cluater.csv")
p=ggplot( H0M_logFC1, aes(x = Human_logFC1, y = Mouse_logFC1, color = cluster)) +
  geom_point(size=1) +
  #geom_smooth(method=lm)+
  #scale_color_uchicago() +
  scale_colour_manual(
    labels=c("Cluster 1","Cluster 2",
             "Cluster 3","Cluster 4",
             "Cluster 5"), 
    values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+
  
  labs(x="Human RNA log2 fold change",y="Mouse RNA log2 fold change",title="Correlation of CRC RNA fold change")+
  annotate("text", x = -7.8, y =6, label = "r = 0.34",
           color="#350E20FF",size = 7 )+
  annotate("text", x = -6.8, y =7, label = " p-value <2e-16",
           color="#350E20FF",size = 7 )+
 
  theme_classic() +
  theme(title = element_text(size = 22,
                             color = "black"),
        axis.title = element_text(size=28),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=17), legend.key.size = unit(0.9,"cm"),#图例大小，间隙
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )+ guides(colour = guide_legend(override.aes = list(size=2.5)))#图例点大小
#cor.test(H0M_logFC1$Human_logFC1,H0M_logFC1$Mouse_logFC1)#  0.3382
p
#ggsave(p, file='/home/panxl/CRC2/human_scatter.pdf', width=8.3, height=6.7) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
################################## 9.2 lim CMS  Cluster 散点图################
df=c("CMS1","CMS2","CMS3","CMS4")
cms_DEG=read.csv("/home/panxl/CRC2/lim_CMS_all_cluater.csv")  
data_text<-data.frame(label=c("r=0.349 ","r=0.303","r=0.308","r=0.342"),
                      type=c("CMS1","CMS2","CMS3","CMS4"),
                      x=c(7,7,7,7),
                      y=c(-5,-5,-5,-5))

p2=ggplot( cms_DEG, aes(x = Human_logFC1, y = Mouse_logFC1, color = cluster,group=type)) +
  geom_point(size=0.9)+geom_smooth(method=lm,color="black",size=0.7)+
  scale_colour_manual(
    labels=c("Cluster 1","Cluster 2",
             "Cluster 3","Cluster 4",
             "Cluster 5"), 
    values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+
  
  labs(x="Human RNA log2 fold change",y="Mouse RNA log2 fold change",title="Correlation of CRC RNA fold change")+
  facet_wrap(~type, ncol = 2) +
  geom_text(data=data_text, mapping=aes(x=x,y=y,label=label,colour=NULL),nudge_x=0.1,nudge_y=0.1)+
  theme_bw() +
  theme(title = element_text(size = 22,
                             color = "black"),
        axis.title = element_text(size=10),
        #axis.text = element_text(size=12, color="black"),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.title.x= element_text(size = 25,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.title.y= element_text(size = 25,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        
        legend.title = element_blank(),
        legend.text = element_text(size=15), legend.key.size = unit(0.9,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )+ guides(colour = guide_legend(override.aes = list(size=2.5)))#图例点大小
p2
cor.test(cms_DEG$Human_logFC1[cms_DEG$type=="CMS1"],cms_DEG$Mouse_logFC1[cms_DEG$type=="CMS1"])# 0.349     0.3498          0.3295  0.3525   0.329
 cor.test(cms_DEG$Human_logFC1[cms_DEG$type=="CMS2"],cms_DEG$Mouse_logFC1[cms_DEG$type=="CMS2"])# 0.3029    0.3034          0.2926 0.3066    0.285
 cor.test(cms_DEG$Human_logFC1[cms_DEG$type=="CMS3"],cms_DEG$Mouse_logFC1[cms_DEG$type=="CMS3"])# 0.3081    0.3086             0.2936 0.3093    0.280
 cor.test(cms_DEG$Human_logFC1[cms_DEG$type=="CMS4"],cms_DEG$Mouse_logFC1[cms_DEG$type=="CMS4"])# 0.3422     0.3426          0.3036 0.3561     0.345
#ggsave(p2, file='/home/panxl/CRC2/cms_scatter.pdf',  width=8.3, height=6.7) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径






















######################################### 10.Human  CMS堆积柱状图#####################################################################
#################### 10.1 堆积柱状图   CMS 4
load("/home/panxl/CRC2/lim_CMS_cluster_number.Rdata")
library(ggplot2)
library("ggsci")
library(dplyr)
human_data=lim_CMS_cluster_number
colnames(human_data)=c("Cluster" ,"Number" , "Type" )
ggplot(data=human_data,aes(x=Type,y=Number))+
  geom_col(color="white",aes(fill=Cluster),position="fill")+
  scale_fill_manual(
    labels=c("Cluster 1","Cluster 2",
             "Cluster 3","Cluster 4",
             "Cluster 5"), #图例标签
    values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
  labs(x="",y="",title="Cluster Distribution in CMS")+
  theme_classic()+
  theme(title = element_text(size = 20,
                             color = "black"),
        
        plot.title = element_text(hjust = 0.5),
        
        axis.text.x = element_text(size = 15 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 20,  
                                   color = "black",
                                   vjust = 0.4))#+
theme_void()
#################### 10.2 堆积柱状图  Human  CMS
H0M_logFC=read.csv("/home/panxl/CRC2/all_cluater.csv")
table(H0M_logFC$cluster)
H0M_logFC_df=as.data.frame(table(H0M_logFC$cluster))
names(H0M_logFC_df) = c("Cluster","Number")
Human=H0M_logFC_df
Human$Type="Human"
data=rbind(human_data,Human)
data$Type=factor(data$Type,levels = c("Human","CMS1","CMS2","CMS3","CMS4"))
p=ggplot(data=data,aes(x=Type,y=Number))+
  geom_col(color="white",aes(fill=Cluster),position="fill")+
  scale_fill_manual(
    labels=c("Cluster 1","Cluster 2",
             "Cluster 3","Cluster 4",
             "Cluster 5"), #图例标签
    values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
  labs(x="",y="",title="Cluster Distribution in CMS")+
  theme_classic()+
  theme(title = element_text(size = 22,
                             color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 15 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 45), #角度
        axis.text.y = element_text(size = 20,  
                                   color = "black",
                                   vjust = 0.4))
  
p
#ggsave(p, file='/home/panxl/CRC2/cms_histogram.pdf', width=9, height=8) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
######################### 10.2 分化cluster2,4 柱状图 ##########################################
load("/home/panxl/CRC2/lim_CMS_cluster_number.Rdata")
library(ggplot2)
library("ggsci")
library(dplyr)
human_data=lim_CMS_cluster_number
colnames(human_data)=c("Cluster" ,"Number" , "Type" )
human_data1=human_data[human_data$Cluster%in%c("Cluster 2" ,"Cluster 4"),]
data=data.frame(Type=c("CMS1","CMS2","CMS3","CMS4"),number=c(1660,1996,1793,1488))
ggplot(data=data, aes(x=Type, y=number, fill = Type)) +
  geom_bar(stat="identity", width=0.5)+
  scale_fill_manual(
    labels=c("Cluster 1","Cluster 2",
             "Cluster 3","Cluster 4",
             "Cluster 5"), #图例标签
    values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
  labs(x="",y="",title="")+
  theme_classic(base_line_size = 1) +
  theme(title = element_text(size = 22),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 17 ,# 修改X轴上字体大小，
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 45), #角度
        axis.text.y = element_text(size = 17,  
                                   vjust = 0.4))


######################################### 11.1 all pathway cor(human-mouse logfc)##############################################################
##############################11.1 ovelap pathway gene
# 载入包
library(clusterProfiler)
library(Rgraphviz)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
# # 转换为ENTREZID
# ##human kegg
# human_data= read.csv("/home/panxl/CRC2/human_DEG.csv")#39440
# resdata<-human_data[abs(human_data$logFC)>=1&human_data$adj.P.Val<=0.05,]#12272
# DEG.gene_ENSEMBL = as.character(resdata$X)
# DEG.entrez_id = mapIds(x = org.Hs.eg.db,
#                        keys = DEG.gene_ENSEMBL,
#                        keytype = "ENSEMBL",
#                        column = "ENTREZID")
# DEG.entrez_id = na.omit(DEG.entrez_id)#去除NA 6250
# kegg <- enrichKEGG(gene       = DEG.entrez_id,
#                    keyType       = "kegg",
#                    organism      = "hsa",
#                    pAdjustMethod = "BH",
#                    pvalueCutoff  = 0.05)
# human_kegg_result=kegg@result#331
# #write.csv(human_kegg_result,"/home/panxl/CRC2/TCGA_kegg_result1.csv")
# 
# human_kegg_result1=human_kegg_result[human_kegg_result$pvalue<=0.05,]#92 26
# human_kegg_result_id=human_kegg_result1$ID#92
# human_kegg_result_id=gsub("hsa","",human_kegg_result_id)#331 92
# ##mouse  kegg
# all_mouse_DEG=read.csv("/home/panxl/CRC2/mouse_DEG.csv",row.names = 1)#21235
# res_mouse<-all_mouse_DEG[ abs(all_mouse_DEG$logFC)>=1&all_mouse_DEG$P.Value<=0.05,]#3908
# DEG.gene_ENSEMBL_mouse = as.character(rownames(res_mouse))
# DEG.entrez_id_mouse = mapIds(x = org.Mm.eg.db,
#                              keys = DEG.gene_ENSEMBL_mouse,
#                              keytype = "ENSEMBL",
#                              column = "ENTREZID")
# DEG.entrez_id_mouse = na.omit(DEG.entrez_id_mouse)#去除NA 3381
# kegg_mouse <- enrichKEGG(gene       = DEG.entrez_id_mouse,
#                          keyType       = "kegg",
#                          organism      = "mmu",
#                          pAdjustMethod = "BH",
#                          pvalueCutoff  = 0.05)
# kegg_mouse_result=kegg_mouse@result #318
# kegg_mouse_result1=kegg_mouse_result[kegg_mouse_result$pvalue<=0.05,]#98
# mouse_kegg_id=kegg_mouse_result1$ID
# mouse_kegg_id=gsub("mmu","",mouse_kegg_id)#310 98
# ##overlap pathway
# overlap_pathway=intersect(human_kegg_result_id,mouse_kegg_id)#64
# save(overlap_pathway, file = "/home/panxl/CRC2/overlap_pathway.RData")
# mouse_pathway_id=paste("mmu",overlap_pathway,sep = "")
# mouse_pathway=kegg_mouse_result[kegg_mouse_result$ID%in%mouse_pathway_id,]#63
# # [1] "04080" "04020" "04974" "04514" "04061" "04976" "04640" "04060" "04713" "00982" "04972" "04672"
# # [13] "05032" "04978" "04726" "00983" "04062" "04024" "05033" "04512" "04975" "04151" "04724" "00830"
# # [25] "04970" "05323" "05207" "00053" "00040" "05204" "04911" "00140" "04610" "00592" "04727" "00590"
# # [37] "05144" "05150" "04923" "02010" "04657" "05224" "04014" "04924" "00565" "05146" "04310" "05226"
# # [49] "04261" "01523" "05321" "05412" "00330" "05414" "00591" "04270" "04961" "00980" "05217" "04964"
# # [61] "04750" "04925" "04950" "04022"##################################11.2  pathway cor
# load(file="/home/panxl/CRC2/overlap_pathway.RData")
library(AnnotationDbi)
library(dplyr)
x <- org.Hs.egPATH #数据库
mapped_genes <- mappedkeys(x)#基因对应的通路
xx <- as.list(x[mapped_genes])
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
k <- keys(org.Hs.eg.db, keytype = "PATH")#全部通路
AnnotationDbi::select(org.Hs.eg.db,
                      keys = k,
                      columns = c("ENTREZID","ENSEMBL"),
                      keytype="PATH")
kegg_annotation <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = overlap_pathway,#stringr::str_replace(mouse_pathway$ID,"hsa",""),
                                         columns = c("ENSEMBL"),
                                         keytype="PATH")#3932

human_data= read.csv("/home/panxl/CRC2/human_DEG.csv",row.names = 1)
human_data$Row.names=rownames(human_data)
all_mouse_DEG=read.csv("/home/panxl/CRC2/mouse_DEG.csv",row.names = 1)
all_mouse_DEG$Row.names=rownames(all_mouse_DEG)
HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]

names(kegg_annotation)[2]="Gene.stable.ID"
hom_kegg = dplyr::left_join(kegg_annotation, HOM_one2one)
hom_kegg <- na.omit(hom_kegg)#2456
colnames(all_mouse_DEG)[1] = "Mouse_logFC1"
hom_kegg = dplyr::left_join(hom_kegg, all_mouse_DEG, c("Mouse.gene.stable.ID" = "Row.names"))
colnames(human_data)[3] = "Human_logFC1"
hom_kegg = dplyr::left_join(hom_kegg, human_data, c("Gene.stable.ID" = "Row.names"))
hom_kegg <- dplyr::select(hom_kegg, c("PATH","Gene.stable.ID","Mouse.gene.stable.ID","Mouse_logFC1","Human_logFC1"))
hom_kegg=na.omit(hom_kegg)#2180
save(hom_kegg,file ="/home/panxl/CRC2/hom_kegg.RData")  
cor_kegg <- purrr::map(unique(hom_kegg$PATH),function(x){
  tmp_df <- hom_kegg[hom_kegg$PATH==x,]
  tmp_df <- na.omit(tmp_df)
  a = cor.test(tmp_df$Mouse_logFC1,tmp_df$Human_logFC1)
  #a$p.value
  re_df <- data.frame("kegg_id"=x,"cor"=a$estimate,"pvalue"=a$p.value)
  row.names(re_df)<-x
  return(re_df)
})%>%dplyr::bind_rows()
#write.csv(cor_kegg,"/home/panxl/CRC2/cor_kegg_human_pathway.csv")

# ###05332 foldchange plot 
# pathway05332_data=hom_kegg[hom_kegg$PATH=="04640",]
# ggplot(pathway05332_data, aes( x = Human_logFC1,y = Mouse_logFC1)) +
#   geom_point(size=1) +
#   #geom_abline(slope=1,intercept=0,color="red",size=1)+
#   #scale_color_uchicago() +
#   labs(x="Human gene log2 fold change",y="Mouse gene log2 fold change",title="04640 Pathway Correlation")+
#   #xlab("Mouse gene log2 fold change") +
#   #ylab("Human gene log2 fold change") + 
#   theme_classic() +
#   theme(axis.title = element_text(size=28),
#         #axis.text = element_text(size=24, color="black"),
#         axis.text.x = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
#         axis.text.y = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
#         legend.title = element_blank(),
#         legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
#         axis.ticks.length.x = unit(2,"mm"),
#         axis.ticks.length.y = unit(2,"mm")
#   )#+geom_text(data=ann,aes(label=Gene.name))
# #correlation
# cor(pathway05332_data$Mouse_logFC1, pathway05332_data$Human_logFC1) # 0.7929954
######################################### 11.2  four cluster pathway ##########
## four cluster
df_data = read.csv("/home/panxl/CRC2/all_cluater.csv")
df=c("Cluster 1","Cluster 2",
     "Cluster 3","Cluster 4")
cluster4_pathway.kegg=lapply(df, function(y){
  z=df_data[df_data$cluster ==y,] 
  cluster_df = df_data[df_data$cluster ==y, "Gene.stable.ID"] %>% as.character()
  DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                         keys = cluster_df,
                         keytype = "ENSEMBL",
                         column = "ENTREZID")
  DEG.entrez_id = na.omit(DEG.entrez_id)#去除NA 6574
  kegg <- enrichKEGG(gene       = DEG.entrez_id,
                     keyType       = "kegg",
                     organism      = "hsa",
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05)
  human_kegg_result=kegg@result
  human_kegg_result$type=y
  return( human_kegg_result)  
})%>%dplyr::bind_rows()
save(cluster4_pathway.kegg, file = "/home/panxl/CRC2/cluster4_pathway.kegg.RData")
#write.csv(result, "/home/panxl/CRC_Conservation/RNA_Seq/TCGA_cluster4_pathway.kegg.csv")
#########################12.1四个cluster各取前5的通路条形图########################
load(file = "/home/panxl/CRC2/cluster4_pathway.kegg.RData")
cluster4_pathway.kegg=cluster4_pathway.kegg[cluster4_pathway.kegg$pvalue<=0.05,]
#write.csv(cluster4_pathway.kegg, "/home/panxl/CRC2/cluster4_pathway.kegg.csv")
df=c("Cluster 1","Cluster 2",
     "Cluster 3","Cluster 4")
cluster_pathway=cluster4_pathway.kegg[,c(1,2,6,9,10)]
cluster1_pathway=cluster_pathway[cluster_pathway$type=="Cluster 1"&cluster_pathway$p.adjust<=0.05,][c(1:25),]# 25
cluster2_pathway=cluster_pathway[cluster_pathway$type=="Cluster 2"&cluster_pathway$p.adjust<=0.05,][c(1:25),]#25 
cluster3_pathway=cluster_pathway[cluster_pathway$type=="Cluster 3"&cluster_pathway$p.adjust<=0.05,][c(1:25),]#25
cluster4_pathway=cluster_pathway[cluster_pathway$type=="Cluster 4"&cluster_pathway$p.adjust<=0.05,]#6

#对富集结果按照p.adjust进行从小到大排序，保证最显著的通路在前
cluster1_pathway <-cluster1_pathway [order(cluster1_pathway $p.adjust),]
cluster1_top5 <-cluster1_pathway[1:5,]
cluster2_pathway <-cluster2_pathway [order(cluster2_pathway $p.adjust),]
cluster2_top5 <-cluster2_pathway[1:5,]
cluster3_pathway <-cluster3_pathway [order(cluster3_pathway $p.adjust),]#%>%filter(ONTOLOGY=="BP")
cluster3_top5 <-cluster3_pathway[1:5,]
cluster4_pathway <-cluster4_pathway [order(cluster4_pathway $p.adjust),]#%>%filter(ONTOLOGY=="BP")
cluster4_top5 <-cluster4_pathway[1:5,]
cluster_top20=rbind(cluster1_top5,cluster2_top5,cluster3_top5,cluster4_top5)

cluster_top20 <- data.frame(cluster_top20$Description,cluster_top20$Count,cluster_top20$p.adjust)
colnames(cluster_top20) <- c("Description","count","p.adjust")
cluster_top20$type=c(rep("cluster1",5),rep("cluster2",5),rep("cluster3",5),rep("cluster4",5))
cluster_top20$type= factor(cluster_top20$type, levels = paste("cluster",c(1:4),sep=""))
# cluster_top20=cluster_top20[cluster_top20$type,]
cluster_top20$Description=factor(cluster_top20$Description,levels = unique(c(cluster_top20$Description)))
cluster_top20$Description=factor(cluster_top20$Description,levels =unique( sort(cluster_top20$Description,decreasing = T)))


# cluster_top20%>%mutate(Description = fct_reorder(Description, -p.adjust))


p=ggplot(data=cluster_top20,aes(x=Description,y=-log10(p.adjust),fill=type))+#fill=padj fill颜色填充，使用连续值padj 
  geom_bar(stat="identity") + coord_flip()+#coord_flip()颠倒坐标轴
  labs(x="",y="-log10 (P_value)",title="Cluster Pathway")+
  scale_fill_manual(
    labels=c("Cluster 1","Cluster 2",
             "Cluster 3","Cluster 4",
             "Cluster 5"), #图例标签
    values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
  theme_classic(base_line_size = 0.61) + 
  theme(axis.title.x = element_text(size =29, 
                                    color = "black"),
        #face = "italic"),
        axis.title.y = element_text(size = 23,
                                    color = "black",
                                    #face = "italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        title = element_text(size = 20,
                             color = "black"),
        legend.title = element_blank(),
        legend.text= element_text(size=12), 
        legend.key.size = unit(0.9,"cm"),
        
        axis.text.x = element_text(size = 25 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   #face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 17,  
                                   color = "black",
                                   vjust = 0.4))
p
#ggsave(p, file='/home/panxl/CRC2/cluster_histogram.pdf', width=10, height=6) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
#################12.2 四个cluster各取前5的通路点图################
load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")
cluster_pathway =cluster4_pathway.kegg
colnames(cluster_pathway)
cluster_pathway=cluster_pathway[,c(1,2,6,9,10)]
cluster1_pathway=cluster_pathway[cluster_pathway$type=="Cluster 1"&cluster_pathway$p.adjust<=0.05,][c(1:5),]#26
cluster2_pathway=cluster_pathway[cluster_pathway$type=="Cluster 2"&cluster_pathway$p.adjust<=0.05,][c(1:5),]#[18
cluster3_pathway=cluster_pathway[cluster_pathway$type=="Cluster 3"&cluster_pathway$p.adjust<=0.05,][c(1:5),]#[c(1:25),]#68  25
cluster4_pathway=cluster_pathway[cluster_pathway$type=="Cluster 4"&cluster_pathway$p.adjust<=0.05,][c(1:5),]#6
c12=rbind(cluster1_pathway,cluster2_pathway)  
c123=rbind(c12,cluster3_pathway)  
c1234=unique(rbind(c123,cluster4_pathway)  )
#cluster4_pathway.kegg=c1234[,c(2,3,4,5)]
cluster4_pathway.kegg=cluster_pathway[cluster_pathway$ID%in%c1234$ID,][,c(2,6,9,10)]

colnames(cluster4_pathway.kegg)
colnames(cluster4_pathway.kegg)=c( "paths", "p.adjust", "count" , "cluster"  )
cluster4_pathway.kegg$cluster=gsub("Cluster 1: Up-regulated in both species","Cluster1",cluster4_pathway.kegg$cluster)    
cluster4_pathway.kegg$cluster=gsub("Cluster 2: Up-regulated in mouse reverse in human","Cluster2",cluster4_pathway.kegg$cluster)    
cluster4_pathway.kegg$cluster=gsub("Cluster 3: Down-regulated in both species","Cluster3",cluster4_pathway.kegg$cluster)    
cluster4_pathway.kegg$cluster=gsub("Cluster 4: Up-regulated in human reverse in mouse","Cluster4",cluster4_pathway.kegg$cluster)    
cluster4_pathway.kegg=cluster4_pathway.kegg[cluster4_pathway.kegg$p.adjust<=0.05,]
#cluster4_pathway.kegg$paths=factor(cluster4_pathway.kegg$paths,levels = c(cluster4_pathway.kegg$paths))
#cluster4_pathway.kegg$paths=factor(cluster4_pathway.kegg$paths,levels =unique( sort(c1234$Description,decreasing = T)))
cluster4_pathway.kegg$paths=factor(cluster4_pathway.kegg$paths,levels =c1234$Description)

p=ggplot(cluster4_pathway.kegg,aes(x=cluster,y=paths))+
  theme_classic() +
  geom_point(aes(size=count,color=p.adjust,alpha=0.5))+
  scale_size(range = c(1, 5))+
  
  #scale_color_gradient2(low ="#D3673C",mid="#67B738",high="#8AC1ED")+
  scale_color_gradient(low ="#D3673C",high="#8AC1ED")+
  
  # labs(x="Cluster category",y="KEGG pathways")+
  labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 14, hjust = 1,
                                   face = "plain",colour = "#101010c2"))+
  theme(text = element_text(size=15,colour = "black",face = "plain"),
        legend.position = "right",
        legend.key.height = unit(30, "pt"),
        legend.key.width = unit(20, "pt"),
        axis.line=element_line(size=1),
        axis.text.y = element_text(size=15,face = "plain",colour = "#101010c2"),
        plot.margin = unit(c(1,1,1,1),"cm"))
p
#ggsave(p, file='/home/panxl/CRC2/pathway20.pdf', width=8, height=6.3) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
#################12.2 四个cluster各取5个的通路点图################
df=data.frame(id=c("hsa03030","hsa03430","hsa04110","hsa04310","hsa03440",
     "hsa05330","hsa05340","hsa04672","hsa04060","hsa04062",
     "hsa04020","hsa04080","hsa04024","hsa04713","hsa04725",
     "hsa00591","hsa00564","hsa04929","hsa00561","hsa00051"
))

load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")
cluster_pathway1 =cluster4_pathway.kegg[cluster4_pathway.kegg$pvalue<=0.05,]
cluster_pathway =cluster_pathway1[cluster_pathway1$ID%in%df$id,]
dd=unique(cluster_pathway[,c(1,2)])
ll=left_join(df,dd,c("id"="ID"))
ll$Description
cluster_pathway2=left_join(ll,cluster_pathway,c("id"="ID"))[,c(1,2,7,10,11)]
cluster_pathway2$Description.x=factor(cluster_pathway2$Description.x,levels = unique(c(cluster_pathway2$Description.x)))

ggplot(cluster_pathway2,aes(x=type,y=Description.x))+
  theme_classic() +
  geom_point(aes(size=Count,color=p.adjust,alpha=0.5))+
  scale_size(range = c(1, 5))+
  
  #scale_color_gradient2(low ="#D3673C",mid="#67B738",high="#8AC1ED")+
  scale_color_gradient(low ="#D3673C",high="#8AC1ED")+
  
  # labs(x="Cluster category",y="KEGG pathways")+
  labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 14, hjust = 1,
                                   face = "plain",colour = "#101010c2"))+
  theme(text = element_text(size=15,colour = "black",face = "plain"),
        legend.position = "right",
        legend.key.height = unit(30, "pt"),
        legend.key.width = unit(20, "pt"),
        axis.line=element_line(size=1),
        axis.text.y = element_text(size=15,face = "plain",colour = "#101010c2"),
        plot.margin = unit(c(1,1,1,1),"cm"))

####################################12.3 四个cluster 唯一通路以及包含基因##############
load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")
cluster_pathway =cluster4_pathway.kegg
cluster_pathway=cluster_pathway[,c(1,2,6,10)]
cluster1_pathway=cluster_pathway[cluster_pathway$type=="Cluster 1"&cluster_pathway$p.adjust<=0.05,]#[c(1:25),]#26
cluster2_pathway=cluster_pathway[cluster_pathway$type=="Cluster 2"&cluster_pathway$p.adjust<=0.05,]#[18
cluster3_pathway=cluster_pathway[cluster_pathway$type=="Cluster 3"&cluster_pathway$p.adjust<=0.05,][c(1:25),]#[c(1:25),]#68  25
cluster4_pathway=cluster_pathway[cluster_pathway$type=="Cluster 4"&cluster_pathway$p.adjust<=0.05,]#7   10
intersect(cluster1_pathway$ID,cluster2_pathway$ID)
l12=c("hsa04060","hsa04151", "hsa04610")
cluster1_pathway=cluster1_pathway[!(cluster1_pathway$ID%in%cluster2_pathway$ID),]# 23"hsa04060"cluste2 --1 19
cluster2_pathway=cluster2_pathway[!(cluster2_pathway$ID%in%c("hsa04060","hsa04061")),]#23
intersect(cluster1_pathway$ID,cluster3_pathway$ID)
intersect(cluster1_pathway$ID,cluster4_pathway$ID)
intersect(cluster2_pathway$ID,cluster3_pathway$ID)
l23= "hsa04020"
cluster2_pathway=cluster2_pathway[!(cluster2_pathway$ID%in%c("hsa04020")),]#15
cluster3_pathway=cluster3_pathway[!(cluster3_pathway$ID%in%c("hsa04020")),]#24
intersect(cluster2_pathway$ID,cluster4_pathway$ID)
intersect(cluster3_pathway$ID,cluster4_pathway$ID)
l34=c("hsa04080" ,"hsa04911" ,"hsa04724")
cluster3_pathway=cluster3_pathway[!(cluster3_pathway$ID%in%l34),]#49
cluster4_pathway=cluster4_pathway[!(cluster4_pathway$ID%in%l34),]#5
unique_cluster_pathway1=rbind(cluster1_pathway,cluster2_pathway)#49
unique_cluster_pathway2=rbind(unique_cluster_pathway1,cluster3_pathway)#70
unique_cluster_pathway=rbind(unique_cluster_pathway2,cluster4_pathway)#73
unique_cluster_pathway=na.omit(unique_cluster_pathway)#73
#write.csv(unique_cluster_pathway, "/home/panxl/CRC2/unique_cluster4_pathway.kegg.csv")
########################## unique_cluster_pathway gene
unique_cluster_pathway=read.csv("/home/panxl/CRC2/unique_cluster4_pathway.kegg.csv")
pathway=unique(unique_cluster_pathway$ID)
unique_cluster_pathway=na.omit(unique_cluster_pathway)#73
pathway=unique(unique_cluster_pathway$ID[!(unique_cluster_pathway$ID=="hsa01523")])#hsa01523

org=keggList("organism")
result=lapply(pathway, function(x){
  gs<-keggGet(x)
  gene=unlist(lapply(gs[[1]]$GENE,function(l){strsplit(l,";")}))
  geneList=gene[1:length(gene)%%3==2]
  geneList=data.frame(geneList)
  geneList$type=x
  return(geneList)
})dplyr::bind_rows()#7594
#write.csv(result,file = "/home/panxl/CRC2/kegg_annotation.csv",row.names = F)

############################# 13.1 Tission species HK None 堆积柱状图  ##############################################
##################### 13.1柱状图
df_var=read.csv("/home/panxl/CRC/RNA_seq/var.csv")
df_species=df_var[df_var$state=="High across species",]#744
df_intestine = read.csv("/home/shimw/project/enhancer_map/conserved/colon-specific.csv")#174
df_tissues=df_var[df_var$state=="High across tissues",]#1884
df_none=df_var[df_var$state=="None",]#11862
##Tission species
library(dplyr)
H0M_logFC1=read.csv("/home/panxl/CRC2/all_cluater")#13736
table(H0M_logFC1$cluster)
H0M_logFC1=unique(H0M_logFC1)#13736
intestine_gene=H0M_logFC1[H0M_logFC1$Mouse.gene.stable.ID%in%df_intestine$Gene,]#170
intestine_gene_df=as.data.frame(table(intestine_gene$cluster))
names(intestine_gene_df) = c("cluster","number")
intestine_gene_df$type="Intestine"

#change id 
ensembl_ID_change=read.csv("/home/panxl/ensembl_37gene.txt")#28255
a=ensembl_ID_change[ensembl_ID_change$Gene.stable.ID%in%H0M_logFC1$Gene.stable.ID,]#13668
b=H0M_logFC1[H0M_logFC1$Gene.stable.ID%in%ensembl_ID_change$Gene.stable.ID,]#13668
H0M_logFC1_changeID=left_join(a,b,by="Gene.stable.ID")#13668
df_species = read.csv("/home/shimw/project/enhancer_map/conserved/species-specific.csv")#744
species_gene=H0M_logFC1_changeID[H0M_logFC1_changeID$Gene.name%in%df_species$gene_id,]#615
species_gene_df=as.data.frame(table(species_gene$cluster))
names(species_gene_df) = c("cluster","number")
species_gene_df$type="Species"

###NONE
none_gene=H0M_logFC1_changeID[H0M_logFC1_changeID$Gene.name%in%df_none$gene_id,]#9178
none_gene_df=as.data.frame(table(none_gene$cluster))
names(none_gene_df) = c("cluster","number")
none_gene_df$type="None"


###########tissues-no intestine
library("org.Mm.eg.db")
library(clusterProfiler) 
library(Hmisc)
df_var=read.csv("/home/panxl/CRC/RNA_seq/var.csv")
df_tissues=df_var[df_var$state=="High across tissues",]#1884
y = tolower(df_tissues$gene_id)
h=capitalize(y)#首字母大写
df_tissues$SYMBOL=h
df_intestine = read.csv("/home/shimw/project/enhancer_map/conserved/colon-specific.csv")#174
gene_df_intestine = bitr(df_intestine$Gene, 
                         fromType="ENSEMBL", 
                         toType="SYMBOL", 
                         OrgDb = "org.Mm.eg.db", drop = TRUE)
gene_intestine <- left_join(df_intestine,gene_df_intestine,c("Gene"="ENSEMBL"))#176

df_no_intestine=df_tissues[!(df_tissues$SYMBOL%in%gene_intestine$SYMBOL),]# 1831 53

no_intestine_gene=H0M_logFC1_changeID[H0M_logFC1_changeID$Gene.name%in%df_no_intestine$gene_id,]#1230
no_intestine_gene_df=as.data.frame(table(no_intestine_gene$cluster))
names(no_intestine_gene_df) = c("cluster","number")
no_intestine_gene_df$type="Tissues"


H0M_logFC1=read.csv("/home/panxl/CRC2/all_cluater")
H0M_logFC1_df=as.data.frame(table(H0M_logFC1$cluster))
names(H0M_logFC1_df) = c("cluster","number")
CRC=H0M_logFC1_df
CRC$type="Human"


data_1=rbind(CRC,intestine_gene_df)
data2=rbind(data_1,species_gene_df)
data3=rbind(data2,no_intestine_gene_df)

data4=rbind(data3,none_gene_df)
data4$type=factor(data4$type,levels=c("Human","None","Intestine","Tissues","Species"))
data4$cluster=gsub("Cluster 1: Up-regulated in both species","Cluster 1",data4$cluster)
data4$cluster=gsub("Cluster 2: Up-regulated in mouse reverse in human","Cluster 2",data4$cluster)
data4$cluster=gsub("Cluster 3: Down-regulated in both species","Cluster 3",data4$cluster)

data4$cluster=gsub("Cluster 4: Up-regulated in human reverse in mouse","Cluster 4",data4$cluster)
data4$cluster=gsub("Cluster 5: Not-significant in both species","Cluster 5",data4$cluster)
colnames(data4)=c("Cluster" ,"Number" , "Type")
homekeeping=data.frame(Cluster=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5"),
                       Number=c("194","25","118","39","2963"))
homekeeping$Type="HK"
homekeeping$Number=as.integer(homekeeping$Number)

data5=rbind(data4,homekeeping)
data5$Type=factor(data5$Type,levels=c("Human","None","Intestine","Tissues","Species","HK"))

p=ggplot(data=data5,aes(x=Type,y=Number))+
  geom_col(color="white",aes(fill=Cluster),position="fill")+
  scale_fill_manual(
    labels=c("Cluster 1","Cluster 2",
              "Cluster 3","Cluster 4",
              "Cluster 5"), 
    values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
  labs(x="",y="",title="")+
  theme_classic()+
  theme(title = element_text(size = 22,
                             color = "black"),
        text = element_text(size=20),#标签大小
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 20 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 45), #角度
        axis.text.y = element_text(size = 17,  
                                   color = "black",
                                   vjust = 0.4))#+
p

ggsave(p, file='/home/panxl/CRC2/tission_histogram.pdf', width=9, height=8) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
################################## 13.2 intestine specific gene fisher.test###############################################
###################19.1.mouse intestine gene
df_intestine = read.csv("/home/shimw/project/enhancer_map/conserved/colon-specific.csv")
H0M_logFC1= read.csv("/home/panxl/CRC2/all_cluater")
H0M_logFC1$cluster=gsub("Cluster 1: Up-regulated in both species","Cluster 1",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 2: Up-regulated in mouse reverse in human","Cluster 2",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 3: Down-regulated in both species","Cluster 3",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 4: Up-regulated in human reverse in mouse","Cluster 4",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 5: Not-significant in both species","Cluster 5",H0M_logFC1$cluster)
df_data=H0M_logFC1
colnames(df_data)
cgroups =unique(df_data$cluster)
#cgroups ="Cluster 1"
odds_ratios <- purrr::map(cgroups,function(cgroups){
  genes_fisher_cluster = matrix(c(nrow(df_data[df_data$cluster==cgroups&df_data$Mouse.gene.stable.ID%in%df_intestine$Gene,]),
                                  nrow(df_data[df_data$cluster==cgroups&!(df_data$Mouse.gene.stable.ID%in%df_intestine$Gene),]),
                                  nrow(df_data[df_data$cluster!=cgroups&df_data$Mouse.gene.stable.ID%in%df_intestine$Gene,]),
                                  nrow(df_data[df_data$cluster!=cgroups&!(df_data$Mouse.gene.stable.ID%in%df_intestine$Gene),])),
                                nrow = 2,
                                dimnames =list(c("housekeeping", "non_housekeeping"),
                                               c("cluster1", "Not cluster1")))
  
  
  OR_1 = fisher.test(genes_fisher_cluster, conf.level = 0.95)
  #return(tibble("boxLabels"=stringr::str_replace(cgroups," ",""),"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
  return(tibble("boxLabels"=cgroups,"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
})dplyr::bind_rows()


#"#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"

p1=ggplot(odds_ratios, aes(x = boxOdds, y = boxLabels)) + 
  geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 1, height = 
                   .15, color = c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95")) +
  geom_point(size = 2, color = c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95")) +
  labs(x="Odds ratio",y="Gene's clusters",title="Intestine specific")+
  # xlab("Odds ratio") +
  # ylab("Gene's clusters") + 
  xlim(0,11)+
  scale_color_uchicago() +
  theme(title = element_text(size = 20,
                             color = "black"),
        
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size=24),
        #axis.title.y = element_blank(),
        axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 17,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 17,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )  
p1
#ggsave(p1, file='/home/panxl/CRC2/mouse_intestin.pdf', width=6.5, height=3.6) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径

################ 19.2 mouse specific gene
df_species = read.csv("/home/shimw/project/enhancer_map/conserved/species-specific.csv")#744
ensembl_ID_change=read.csv("/home/panxl/ensembl_37gene.txt")#28255
H0M_logFC2=left_join(H0M_logFC1,ensembl_ID_change,c("Gene.stable.ID"="Gene.stable.ID"))
H0M_logFC3=na.omit(H0M_logFC2)
df_data=H0M_logFC3
colnames(df_species)

cgroups =unique(df_data$cluster)
#cgroups ="Cluster 1"
odds_specific_ratios <- purrr::map(cgroups,function(cgroups){
  genes_fisher_cluster = matrix(c(nrow(df_data[df_data$cluster==cgroups&df_data$Gene.name%in%df_species$gene_id,]),
                                  nrow(df_data[df_data$cluster==cgroups&!(df_data$Gene.name%in%df_species$gene_id),]),
                                  nrow(df_data[df_data$cluster!=cgroups&df_data$Gene.name%in%df_species$gene_id,]),
                                  nrow(df_data[df_data$cluster!=cgroups&!(df_data$Gene.name%in%df_species$gene_id),])),
                                nrow = 2,
                                dimnames =list(c("housekeeping", "non_housekeeping"),
                                               c("cluster1", "Not cluster1")))
  
  
  OR_1 = fisher.test(genes_fisher_cluster, conf.level = 0.95)
  #return(tibble("boxLabels"=stringr::str_replace(cgroups," ",""),"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
  return(tibble("boxLabels"=cgroups,"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
})dplyr::bind_rows()


#"#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"

p2=ggplot(odds_specific_ratios , aes(x = boxOdds, y = boxLabels)) + 
  geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 1, height = 
                   .15, color = c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95")) +
  geom_point(size = 2, color = c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95")) +
  scale_color_uchicago() +
  labs(x="Odds ratio",y="Gene's clusters",title="Species specific")+
  theme(title = element_text(size = 20,
                             color = "black"),
        
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size=24),
        #axis.title.y = element_blank(),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 18,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 18,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )  
p2
#ggsave(p2, file='/home/panxl/CRC2/mouse_specific.pdf', width=6.5, height=3.6) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径


######################################13.3 Tission species None 分布散点图##############################################
library("org.Mm.eg.db")
library(clusterProfiler) 
library(Hmisc)
df_var=read.csv("/home/panxl/CRC/RNA_seq/var.csv")##14490
#df_tissues=df_var[df_var$state=="High across tissues",]#1884
y = tolower(df_var$gene_id)
h=capitalize(y)#首字母大写
df_var$SYMBOL=h
df_intestine = read.csv("/home/shimw/project/enhancer_map/conserved/colon-specific.csv")#174
gene_df_intestine = bitr(df_intestine$Gene, 
                         fromType="ENSEMBL", 
                         toType="SYMBOL", 
                         OrgDb = "org.Mm.eg.db", drop = TRUE)
gene_intestine <- left_join(df_intestine,gene_df_intestine,c("Gene"="ENSEMBL"))#176

intestine=df_var[df_var$SYMBOL%in%gene_intestine$SYMBOL,c(1,2,4,5)]#150
intestine$state="Intestinal specificity"
df_no_intestine=df_var[!(df_var$SYMBOL%in%gene_intestine$SYMBOL),c(1,2,4,5,3)]# 14340
df_data=rbind(df_no_intestine,intestine)#14490
df_data$state <- factor(df_data$state, levels = c("High across tissues","High across species","Intestinal specificity","None"))
library(ggplot2)
library("ggsci")
library(RColorBrewer)
#display.brewer.all()
ggplot( df_data, aes(x = tissues, y = species, color = state)) +
  geom_point(size=2,) +
  scale_color_manual(values = c("#D3A46E","#787B87","#8E6A4B","#99A192"))+
  labs(x="Tissues",y="Species",title="")+
  theme_classic() +
  theme(title = element_text(size = 22,
                             color = "black"),
        axis.title = element_text(size=28),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=24), legend.key.size = unit(1,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )+theme(legend.position = c(0.8,0.9))

####################### 13.3 Tission species logfc分布散点图 ##########################
####intestine
df_intestine = read.csv("/home/shimw/project/enhancer_map/conserved/colon-specific.csv")#174
library(dplyr)
H0M_logFC1=read.csv("/home/panxl/CRC2/all_cluater.csv")#13736
table(H0M_logFC1$cluster)
H0M_logFC1=unique(H0M_logFC1)#13736
intestine_gene=H0M_logFC1[H0M_logFC1$Mouse.gene.stable.ID%in%df_intestine$Gene,]#170
intestine_gene_df=as.data.frame(table(intestine_gene$cluster))
names(intestine_gene_df) = c("cluster","number")
intestine_gene_df$type="Intestine"

p1=ggplot( intestine_gene, aes(x = Human_logFC1, y = Mouse_logFC1, color= cluster)) +
  geom_point(size=1) +
  scale_colour_manual(
    labels=c("Cluster 1","Cluster 2",
             "Cluster 3","Cluster 4",
             "Cluster 5"), 
    values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+
  
  labs(x="Human gene log2 fold change",y="Mouse gene log2 fold change",title="Intestine gene Correlation")+
  annotate("text", x = 1.5, y =-3.5, label = "r = 0.403 p-value = 5e-08",
           color="#350E20FF",size = 6 )+
  theme_classic() +
  theme(title = element_text(size = 22,
                             color = "black"),
        axis.title = element_text(size=20),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 20,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 20,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=15), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )
p1#786x647
cor(intestine_gene$Human_logFC1,intestine_gene$Mouse_logFC1)#0.4033 p-value = 5e-08
cor.test(intestine_gene$Human_logFC1,intestine_gene$Mouse_logFC1)
#ggsave(p1, file='/home/panxl/CRC2/intestine_scatter .pdf', width=9, height=6) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
  
####Species gene Correlation
df_species=df_var[df_var$state=="High across species",]#744
H0M_logFC1=read.csv("/home/panxl/CRC2/all_cluater.csv")#13736
table(H0M_logFC1$cluster)
H0M_logFC1=unique(H0M_logFC1)#13736

#change id 
ensembl_ID_change=read.csv("/home/panxl/ensembl_37gene.txt")#28255
a=ensembl_ID_change[ensembl_ID_change$Gene.stable.ID%in%H0M_logFC1$Gene.stable.ID,]#13668
b=H0M_logFC1[H0M_logFC1$Gene.stable.ID%in%ensembl_ID_change$Gene.stable.ID,]#13668
H0M_logFC1_changeID=left_join(a,b,by="Gene.stable.ID")#13668
df_species = read.csv("/home/shimw/project/enhancer_map/conserved/species-specific.csv")#744
species_gene=H0M_logFC1_changeID[H0M_logFC1_changeID$Gene.name%in%df_species$gene_id,]#615
species_gene_df=as.data.frame(table(species_gene$cluster))
names(species_gene_df) = c("cluster","number")
species_gene_df$type="Species"

p2=ggplot( species_gene, aes(x = Human_logFC1, y = Mouse_logFC1, color= cluster)) +
  geom_point(size=1) +
  scale_colour_manual(
    labels=c("Cluster 1","Cluster 2",
             "Cluster 3","Cluster 4",
             "Cluster 5"), 
    values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+
  
  labs(x="Human gene log2 fold change",y="Mouse gene log2 fold change",title="Species gene Correlation")+
  annotate("text", x = 0.5, y =-3, label = "r = 0.202 p-value = 4e-07",
           color="#350E20FF",size = 7 )+
  theme_classic() +
  theme(title = element_text(size = 22,
                             color = "black"),
        axis.title = element_text(size=20),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 20,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 20,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=15), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )
p2
#ggsave(p2, file='/home/panxl/CRC2/species_scatter .pdf', width=9, height=6) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径

cor(species_gene$Human_logFC1,species_gene$Mouse_logFC1)# 0.2019
cor.test(species_gene$Human_logFC1,species_gene$Mouse_logFC1)
#################################### 14. gseKEGG NES  ###############################################################
################################ 14.1同源基因 human GSEA_KEGG#################################
hom=read.csv("/home/panxl/CRC2/all_cluater")
human=read.csv(paste("/home/panxl/CRC2/human_DEG.csv",sep=""))
human_gene=human[human$X%in%hom$Gene.stable.ID,]
gene_res_ENSEMBL = as.character(human_gene$X)
gene_res_df = bitr(gene_res_ENSEMBL, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_df <- left_join(human_gene,gene_res_df,c("X"="ENSEMBL"))#
gene_df1=na.omit(gene_df)
gene_res0.01= as.vector(gene_df1$ENTREZID)#

geneList<-gene_res0.01
data(geneList,package = "DOSE")
geneList = gene_df$logFC
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList,decreasing = T)
geo_list<- list()
pvalue_Cutoff=0.05
human_gse_kegg<-gseKEGG(geneList= geneList,
                        organism = "hsa",
                        keyType  = 'kegg',
                        nPerm  = 1000,
                        minGSSize = 10,
                        maxGSSize = 500,
                        pvalueCutoff = 1,
                        pAdjustMethod     = "BH") ## 0.05有的没有结果 
#save(human_gse_kegg,file = "/home/panxl/CRC2/human_gse_kegg.Rdata")
################################ 14.2同源基因 CMS GSEA_KEGG#################################
df=c("CMS1","CMS2","CMS3","CMS4")  
hom=read.csv("/home/panxl/CRC2/all_cluater")
resdata=read.csv(paste("/home/panxl/CRC2/CMS4_DEG.csv"))#亚型重复四次
resdata=resdata[resdata$X%in%hom$Gene.stable.ID,]
gene_res_ENSEMBL = as.character(resdata$X)
gene_res_df = bitr(gene_res_ENSEMBL, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_df <- left_join(resdata,gene_res_df,c("X"="ENSEMBL"))#
gene_df1=na.omit(gene_df)
gene_res0.01= as.vector(gene_df1$ENTREZID)#
geneList<-gene_res0.01
data(geneList,package = "DOSE")
geneList = gene_df$logFC
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList,decreasing = T)
geo_list<- list()
pvalue_Cutoff=0.05
CMS4_gse_kegg<-gseKEGG(geneList= geneList,
                       organism = "hsa",
                       keyType  = 'kegg',
                       nPerm  = 1000,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 1,
                       pAdjustMethod     = "BH") ## 0.05有的没有结果 

# save(CMS1_gse_kegg,file = "/home/panxl/CRC2/CMS1_gse_kegg.Rdata")
# save(CMS2_gse_kegg,file = "/home/panxl/CRC2/CMS2_gse_kegg.Rdata")
# save(CMS3_gse_kegg,file = "/home/panxl/CRC2/CMS3_gse_kegg.Rdata")
# save(CMS4_gse_kegg,file = "/home/panxl/CRC2/CMS4_gse_kegg.Rdata")
#############################14.3同源基因 Mouse GSEA_KEGG#############
library(clusterProfiler) 
library(org.Mm.eg.db)
library(patchwork)
hom=read.csv("/home/panxl/CRC2/all_cluater")
mouse=read.csv(paste("/home/panxl/CRC2/mouse_DEG.csv",sep=""))
mouse=mouse[mouse$X%in%hom$Mouse.gene.stable.ID,]
gene_res_ENSEMBL = as.character(mouse$X)
gene_res_df = bitr(gene_res_ENSEMBL, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Mm.eg.db", drop = TRUE)
gene_df <- left_join(mouse,gene_res_df,c("X"="ENSEMBL"))#
gene_df1=na.omit(gene_df)
gene_res0.01= as.vector(gene_df1$ENTREZID)#
geneList<-gene_res0.01
data(geneList,package = "DOSE")
geneList = gene_df$logFC
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList,decreasing = T)
geo_list<- list()
pvalue_Cutoff=0.05
mouse_gse_kegg<-gseKEGG(geneList= geneList,
                        organism = "mmu",
                        keyType  = 'kegg',
                        nPerm  = 1000,
                        minGSSize = 10,
                        maxGSSize = 500,
                        pvalueCutoff = 1,
                        pAdjustMethod     = "BH") ## 0.05有的没有结果

save(mouse_gse_kegg,file = "/home/panxl/CRC2/mouse_gse_kegg.Rdata")
###############################14.4 同源基因 NES ##################################
load("/home/panxl/CRC2/mouse_gse_kegg.Rdata")
load("/home/panxl/CRC2/human_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS1_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS2_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS3_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS4_gse_kegg.Rdata")
####### unique cluster nes 
CMS1_gse=CMS1_gse_kegg@result[,c(2,5)]
colnames(CMS1_gse)=c("Description","CMS1")
CMS2_gse=CMS2_gse_kegg@result[,c(2,5)]
colnames(CMS2_gse)=c("Description","CMS2")
CMS3_gse=CMS3_gse_kegg@result[,c(2,5)]
colnames(CMS3_gse)=c("Description","CMS3")
CMS4_gse=CMS4_gse_kegg@result[,c(2,5)]
colnames(CMS4_gse)=c("Description","CMS4")
mouse_gse=mouse_gse_kegg@result[,c(2,5)]
colnames(mouse_gse)=c("Description","Mouse")
human_gse=human_gse_kegg@result[,c(2,5)]
colnames(human_gse)=c("Description","Human")
a=left_join(CMS1_gse,CMS2_gse,by = "Description")
b=left_join(CMS3_gse,a,by = "Description")
c=left_join(CMS4_gse,b,by = "Description")
d=left_join(mouse_gse,c,by = "Description")
NES=left_join(human_gse,d,by = "Description")
#####four cluster
unique_cluster_pathway=read.csv("/home/panxl/CRC2/unique_cluster4_pathway.kegg.csv")
load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")
cluster_pathway =cluster4_pathway.kegg
cluster_pathway=cluster_pathway[,c(1,2,6,10)]
cluster1_pathway=cluster_pathway[cluster_pathway$type=="Cluster 1"&cluster_pathway$p.adjust<=0.05,]#[c(1:25),]#26
cluster2_pathway=cluster_pathway[cluster_pathway$type=="Cluster 2"&cluster_pathway$p.adjust<=0.05,]#[18
cluster3_pathway=cluster_pathway[cluster_pathway$type=="Cluster 3"&cluster_pathway$p.adjust<=0.05,][c(1:25),]#[c(1:25),]#68  25
cluster4_pathway=cluster_pathway[cluster_pathway$type=="Cluster 4"&cluster_pathway$p.adjust<=0.05,]#7   10
intersect(cluster1_pathway$ID,cluster2_pathway$ID)
l12=c("hsa04060","hsa04151", "hsa04610")
cluster1_pathway=cluster1_pathway[!(cluster1_pathway$ID%in%cluster2_pathway$ID),]# 23"hsa04060"cluste2 --1 19
cluster2_pathway=cluster2_pathway[!(cluster2_pathway$ID%in%c("hsa04060","hsa04061")),]#23
intersect(cluster1_pathway$ID,cluster3_pathway$ID)
intersect(cluster1_pathway$ID,cluster4_pathway$ID)
intersect(cluster2_pathway$ID,cluster3_pathway$ID)
l23= "hsa04020"
cluster2_pathway=cluster2_pathway[!(cluster2_pathway$ID%in%c("hsa04020")),]#15
cluster3_pathway=cluster3_pathway[!(cluster3_pathway$ID%in%c("hsa04020")),]#24
intersect(cluster2_pathway$ID,cluster4_pathway$ID)
intersect(cluster3_pathway$ID,cluster4_pathway$ID)
l34=c("hsa04080" ,"hsa04911" ,"hsa04724")
cluster3_pathway=cluster3_pathway[!(cluster3_pathway$ID%in%l34),]#49
cluster4_pathway=cluster4_pathway[!(cluster4_pathway$ID%in%l34),]#5

NES$cluster=case_when(
  NES$Description%in%cluster1_pathway$Description~"Cluster1",
  NES$Description%in%cluster2_pathway$Description~"Cluster2",
  NES$Description%in%cluster3_pathway$Description~"Cluster3",
  NES$Description%in%cluster4_pathway$Description~"Cluster4",
)
nes=NES[NES$Description%in%unique_cluster_pathway$Description,]
nes=na.omit(nes)
table(nes$cluster)
#write.csv(nes, "/home/panxl/CRC2/pathway_cms_nes.csv")
################################14.5  NES 热图###############################################################  
library(tidyverse)
library(pheatmap)
NES=read.csv("/home/panxl/CRC2/4cms_nes.csv",row.names = 1)%>%arrange(cluster)
colnames(NES)
data=NES[,c(9,8,2,3,4,5)]
rownames(data)=NES$Description
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

ha1 = HeatmapAnnotation(Group=factor(c("Mouse", "Human",  "CMS1" , "CMS2" , "CMS3" , "CMS4" ),
                                     levels = c("Mouse", "Human",  "CMS1" , "CMS2" , "CMS3" , "CMS4" )),
                        col = list(Group = c("Mouse"="#817936","Human"="#78a355",'CMS1' = '#b97573',
                                             'CMS2' = '#686954', 'CMS3' = '#5b5e66', 'CMS4'= '#816147')),show_legend = F)


ha2 = rowAnnotation(cluster=factor(c(rep("Cluster1",25),rep("Cluster2",24),rep("Cluster3",21),rep("Cluster4",3)),
                                   levels = c("Cluster1","Cluster2","Cluster3","Cluster4")),
                    col = list(cluster=c("Cluster1"="#CB9C7A","Cluster2"="#8696a7",
                                         "Cluster3"="#CDB97D","Cluster4"="#7b8b6f")),show_legend = F)

#df_sub1=t(scale(t(df_sub)))"#A59B95"
#library()
Heatmap(data, cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        col=c("#87CEFA", "white", "#CC2121"),
        top_annotation = ha1,
        left_annotation = ha2,
        row_names_gp = gpar(fontsize = 9), 
        show_heatmap_legend = T)
ht_list=Heatmap(data, cluster_rows = F,
                cluster_columns = F,
                show_column_names = T,
                show_row_names = T,
                heatmap_legend_param = list(
                  title = "NES", at = c(-2, 0,4),
                  labels = c("-2","0", "4")
                ),
            
                col=c("#87CEFA", "white", "#CC2121"),
                top_annotation = ha1,
                left_annotation = ha2,
                column_names_gp  =gpar(
                  col = "black",
                  fontsize = 10
                ),
                row_names_gp = gpar(fontsize = 10),
                column_names_rot =45,
                show_heatmap_legend = T)
draw(ht_list, heatmap_legend_side = "left")
#################################### 14.  fgseaRes NES  ###############################################################
################################ 14.1同源基因 human +cms fgseaRes NES #################################
unique_cluster_pathway=read.csv("/home/panxl/CRC2/unique_cluster4_pathway.kegg.csv")
load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")

#####examplePathways
library(clusterProfiler) 
library(org.Hs.eg.db)
library(KEGGREST)
pathway=unique(unique_cluster_pathway$ID)
unique_cluster_pathway=na.omit(unique_cluster_pathway)#65
pathway=unique(unique_cluster_pathway$ID[!(unique_cluster_pathway$ID=="hsa01523")])#hsa01523

org=keggList("organism")
gs<-keggGet("hsa04110")
gene=unlist(lapply(gs[[1]]$GENE,function(l){strsplit(l,";")}))
geneList=gene[1:length(gene)%%3==2]
geneList=data.frame(geneList)
geneList$type=x
return(geneList)
org=keggList("organism")
unique_cluster_pathway_kegg_annotation=lapply(pathway, function(x){
  gs<-keggGet(x)
  gene=unlist(lapply(gs[[1]]$GENE,function(l){strsplit(l,";")}))
  geneList=gene[1:length(gene)%%3==2]
  geneList=data.frame(geneList)
  geneList$type=x
  return(geneList)
})%>%dplyr::bind_rows()#7594
#save(unique_cluster_pathway_kegg_annotation,file = "/home/panxl/CRC2/unique_cluster_pathway_kegg_annotation.RData")
##########################
load("/home/panxl/CRC2/unique_cluster_pathway_kegg_annotation.RData")
kegg.gene_symbol = as.character(unique(unique_cluster_pathway_kegg_annotation$geneList))#2983
kegg.gene_df = bitr(kegg.gene_symbol,
                    fromType="SYMBOL",
                    toType="ENTREZID",
                    OrgDb = "org.Hs.eg.db", drop = TRUE) #2614

kegg.gene_df=na.omit(kegg.gene_df)
gene_df <-left_join(unique_cluster_pathway_kegg_annotation,kegg.gene_df,c("geneList"="SYMBOL"))#
gene_df=na.omit(gene_df)[,c(2,3)]#7229
gene_df$ENTREZID
examplePathways=lapply(unique(gene_df$type), function(x){
  gene_df[gene_df$type==x,]$ENTREZID
})
names(examplePathways)=unique(gene_df$type)#66
#save(examplePathways,file = "/home/panxl/CRC2/human_examplePathways.RData")
#############################
#human=read.csv(paste("/home/panxl/CRC2/human_DEG.csv",sep=""))
#colnames(human)[2]="log2FoldChange"
#colnames(human)[1]="gene"
#write.csv(human,file ="/home/panxl/CRC2/human_Human.csv" )
hom=read.csv("/home/panxl/CRC2/all_cluater")
resdata=read.csv(paste("/home/panxl/CRC2/CMS1_DEG.csv",sep=""))#亚型重复四次
load("/home/panxl/CRC2/human_examplePathways.RData")
#resdata$X resdata$logFC
library(fgsea)
library(ggplot2)
df=c("CMS1","CMS2","CMS3","CMS4","human")

pathway_cms_nes=lapply(df, function(x){
  y=read.csv(paste("/home/panxl/CRC2/",x,"_DEG.csv",sep=""))
  hom_gene=y[y$X%in%hom$Gene.stable.ID,]
  cms_gene_df=as.character(hom_gene$X)
   cms_gene_df=bitr(cms_gene_df, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  gene_df <- left_join(y,cms_gene_df,c("X"="ENSEMBL"))#
  gene_df =na.omit(gene_df)
  gene_res= as.vector(gene_df$ENTREZID)#
  geneList<-gene_res
  data(geneList,package = "DOSE")
  geneList = gene_df$logFC
  names(geneList) = gene_df$ENTREZID
  geneList = sort(geneList,decreasing = T)
  fgseaRes <- fgsea(pathways = examplePathways, 
                    stats    = geneList,
                    minSize  = 15,
                    maxSize  = 500)
  pathway_nes=fgseaRes[,c(1,6)]
  pathway_nes$type=x
  return(pathway_nes)
})%>%dplyr::bind_rows()

load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")
cluster_pathway =cluster4_pathway.kegg
cluster_pathway=cluster_pathway[,c(1,2,6,10)]
cluster1_pathway=cluster_pathway[cluster_pathway$type=="Cluster 1"&cluster_pathway$p.adjust<=0.05,]#[c(1:25),]#26
cluster2_pathway=cluster_pathway[cluster_pathway$type=="Cluster 2"&cluster_pathway$p.adjust<=0.05,]#[18
cluster3_pathway=cluster_pathway[cluster_pathway$type=="Cluster 3"&cluster_pathway$p.adjust<=0.05,][c(1:25),]#[c(1:25),]#68  25
cluster4_pathway=cluster_pathway[cluster_pathway$type=="Cluster 4"&cluster_pathway$p.adjust<=0.05,]#7   10
intersect(cluster1_pathway$ID,cluster2_pathway$ID)
l12=c("hsa04060","hsa04151", "hsa04610")
cluster1_pathway=cluster1_pathway[!(cluster1_pathway$ID%in%cluster2_pathway$ID),]# 23"hsa04060"cluste2 --1 19
cluster2_pathway=cluster2_pathway[!(cluster2_pathway$ID%in%c("hsa04060","hsa04061")),]#23
intersect(cluster1_pathway$ID,cluster3_pathway$ID)
intersect(cluster1_pathway$ID,cluster4_pathway$ID)
intersect(cluster2_pathway$ID,cluster3_pathway$ID)
l23= "hsa04020"
cluster2_pathway=cluster2_pathway[!(cluster2_pathway$ID%in%c("hsa04020")),]#15
cluster3_pathway=cluster3_pathway[!(cluster3_pathway$ID%in%c("hsa04020")),]#24
intersect(cluster2_pathway$ID,cluster4_pathway$ID)
intersect(cluster3_pathway$ID,cluster4_pathway$ID)
l34=c("hsa04080" ,"hsa04911" ,"hsa04724")
cluster3_pathway=cluster3_pathway[!(cluster3_pathway$ID%in%l34),]#49
cluster4_pathway=cluster4_pathway[!(cluster4_pathway$ID%in%l34),]#5


pathway_cms_nes1=pathway_cms_nes
pathway_cms_nes1$cluster=case_when(
  pathway_cms_nes1$pathway%in%cluster1_pathway$ID~"Cluster1",
  pathway_cms_nes1$pathway%in%cluster2_pathway$ID~"Cluster2",
  pathway_cms_nes1$pathway%in%cluster3_pathway$ID~"Cluster3",
  pathway_cms_nes1$pathway%in%cluster4_pathway$ID~"Cluster4",
)#260
pathway_cms_nes1=na.omit(pathway_cms_nes1)
#write.csv(pathway_cms_nes1, "/home/panxl/CRC2/lim_pathway_cms_nes1.csv")

################################ 14.2同源基因 mouse fgseaRes NES #################################
library(clusterProfiler) 
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(KEGGREST)
unique_cluster_pathway=read.csv("/home/panxl/CRC2/unique_cluster4_pathway.kegg.csv")
pathway=gsub("hsa","mmu",unique_cluster_pathway$ID)
pathway=unique(pathway[!(pathway=="mmu01523")])#hsa01523
org=keggList("organism")
unique_mouse_pathway_kegg_annotation=lapply(pathway, function(x){
  gs<-keggGet(x)
  gene=unlist(lapply(gs[[1]]$GENE,function(l){strsplit(l,";")}))
  geneList=gene[1:length(gene)%%3==2]
  geneList=data.frame(geneList)
  geneList$type=x
  return(geneList)
})%>%dplyr::bind_rows()#7812
#save(unique_mouse_pathway_kegg_annotation,file = "/home/panxl/CRC2/unique_mouse_pathway_kegg_annotation.RData")

kegg_annotation=unique_mouse_pathway_kegg_annotation
kegg.gene_symbol = as.character(unique(kegg_annotation$geneList))#3094
kegg.gene_df = bitr(kegg.gene_symbol,
                    fromType="SYMBOL",
                    toType="ENTREZID",
                    OrgDb = "org.Mm.eg.db", drop = TRUE) #2888
kegg.gene_df=na.omit(kegg.gene_df)#2913
gene_df <-left_join(kegg_annotation,kegg.gene_df,c("geneList"="SYMBOL"))#7812
gene_df=na.omit(gene_df)[,c(2,3)]#7812
gene_df$ENTREZID
examplePathways=lapply(unique(gene_df$type), function(x){
  gene_df[gene_df$type==x,]$ENTREZID
})
names(examplePathways)=unique(gene_df$type)#66
#save(examplePathways,file = "/home/panxl/CRC2/mouse_examplePathways.RData")
##########################
load("/home/panxl/CRC2/mouse_examplePathways.RData")
hom=read.csv("/home/panxl/CRC2/all_cluater")

mouse_deg=read.csv("/home/panxl/CRC2/mouse_DEG.csv")
  hom_gene=mouse_deg[mouse_deg$X%in%hom$Mouse.gene.stable.ID,]
  cms_gene_df=as.character(hom_gene$X)
  cms_gene_df=bitr(cms_gene_df, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Mm.eg.db", drop = TRUE)
  gene_df <- left_join(mouse_deg,cms_gene_df,c("X"="ENSEMBL"))#
  gene_df =na.omit(gene_df)
  gene_res= as.vector(gene_df$ENTREZID)#
  geneList<-gene_res
  data(geneList,package = "DOSE")
  geneList = gene_df$logFC
  names(geneList) = gene_df$ENTREZID
  geneList = sort(geneList,decreasing = T)
  fgseaRes <- fgsea(pathways = examplePathways, 
                    stats    = geneList,
                    minSize  = 15,
                    maxSize  = 500)
  pathway_nes=fgseaRes[,c(1,6)]
  load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")
  
      pathway_nes$type="mouse"
      cluster_pathway =cluster4_pathway.kegg
      cluster_pathway=cluster_pathway[,c(1,2,6,10)]
      cluster1_pathway=cluster_pathway[cluster_pathway$type=="Cluster 1"&cluster_pathway$p.adjust<=0.05,]#[c(1:25),]#26
      cluster2_pathway=cluster_pathway[cluster_pathway$type=="Cluster 2"&cluster_pathway$p.adjust<=0.05,]#[18
      cluster3_pathway=cluster_pathway[cluster_pathway$type=="Cluster 3"&cluster_pathway$p.adjust<=0.05,][c(1:25),]#[c(1:25),]#68  25
      cluster4_pathway=cluster_pathway[cluster_pathway$type=="Cluster 4"&cluster_pathway$p.adjust<=0.05,]#7   10
      intersect(cluster1_pathway$ID,cluster2_pathway$ID)
      l12=c("hsa04060","hsa04151", "hsa04610")
      cluster1_pathway=cluster1_pathway[!(cluster1_pathway$ID%in%cluster2_pathway$ID),]# 23"hsa04060"cluste2 --1 19
      cluster2_pathway=cluster2_pathway[!(cluster2_pathway$ID%in%c("hsa04060","hsa04061")),]#23
      intersect(cluster1_pathway$ID,cluster3_pathway$ID)
      intersect(cluster1_pathway$ID,cluster4_pathway$ID)
      intersect(cluster2_pathway$ID,cluster3_pathway$ID)
      l23= "hsa04020"
      cluster2_pathway=cluster2_pathway[!(cluster2_pathway$ID%in%c("hsa04020")),]#15
      cluster3_pathway=cluster3_pathway[!(cluster3_pathway$ID%in%c("hsa04020")),]#24
      intersect(cluster2_pathway$ID,cluster4_pathway$ID)
      intersect(cluster3_pathway$ID,cluster4_pathway$ID)
      l34=c("hsa04080" ,"hsa04911" ,"hsa04724")
      cluster3_pathway=cluster3_pathway[!(cluster3_pathway$ID%in%l34),]#49
      cluster4_pathway=cluster4_pathway[!(cluster4_pathway$ID%in%l34),]#5
      pathway_nes$pathway=gsub( "mmu","hsa",pathway_nes$pathway)
      pathway_nes$cluster=case_when(
        pathway_nes$pathway%in%cluster1_pathway$ID~"Cluster1",
        pathway_nes$pathway%in%cluster2_pathway$ID~"Cluster2",
        pathway_nes$pathway%in%cluster3_pathway$ID~"Cluster3",
        pathway_nes$pathway%in%cluster4_pathway$ID~"Cluster4",
      )#260  
      pathway_nes=na.omit(pathway_nes)    
#save(pathway_nes,file = "/home/panxl/CRC2/mouse_pathway_nes.RData")
##############NES
load("/home/panxl/CRC2/mouse_pathway_nes.RData")
human_nes=read.csv("/home/panxl/CRC2/lim_pathway_cms_nes1.csv")[,-1]
nes=rbind(human_nes,pathway_nes)
table(nes$pathway)
NES=nes[!(nes$pathway%in%c("hsa04514","hsa05150")),]
df=c("CMS1","CMS2","CMS3","CMS4","Human","mouse")

NES_DATA=data.frame(Pathway=NES[NES$type=="CMS1",]$pathway,CMS1=NES[NES$type=="CMS1",]$NES,
                    CMS2=NES[NES$type=="CMS2",]$NES
                    ,CMS3=NES[NES$type=="CMS3",]$NES
                    ,CMS4=NES[NES$type=="CMS4",]$NES
                    ,Human=NES[NES$type=="human",]$NES
                    ,Mouse=NES[NES$type=="mouse",]$NES
                    ,Cluster=NES[NES$type=="CMS1",]$cluster)
load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")
cluster4_pathway.kegg=cluster4_pathway.kegg[,c(1,2)]
NES3=left_join(NES_DATA,cluster4_pathway.kegg,c("Pathway"="ID"))
NES4=na.omit(NES3)
NES5=unique(NES4)



# NES1=lapply(df, function(x){
#   l=NES[NES$type==x,]
#  return(l)
# })%>%dplyr::bind_cols()
# colnames(NES1)
# NES2=NES1[,c(1,2,6,10,14,18,22,4)]
# colnames(NES2)=c("Pathway","CMS1","CMS2","CMS3","CMS4","Human","mouse","Cluster")
# load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")
# cluster4_pathway.kegg=cluster4_pathway.kegg[,c(1,2)]
# NES3=left_join(NES2,cluster4_pathway.kegg,c("Pathway"="ID"))
# NES4=na.omit(NES3)
# NES5=unique(NES4)
#write.csv(NES5,file = "/home/panxl/CRC2/lim_NES.csv")
############################14.3 fgseaRes NES热图#######################      
NES = read.csv("/home/panxl/CRC2/lim_NES.csv",row.names = 1)%>%arrange(Cluster)
rownames(NES)=NES$Description
colnames(NES)
data=NES[,c(2:7)]
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

ha1 = HeatmapAnnotation(Group=factor(c("Mouse", "Human",  "CMS1" , "CMS2" , "CMS3" , "CMS4" ),
                                     levels = c("Mouse", "Human",  "CMS1" , "CMS2" , "CMS3" , "CMS4" )),
                        col = list(Group = c("Mouse"="#817936","Human"="#78a355",'CMS1' = '#b97573',
                                             'CMS2' = '#686954', 'CMS3' = '#5b5e66', 'CMS4'= '#816147')),show_legend = F)

table(NES$Cluster)
ha2 = rowAnnotation(Cluster=factor(c(rep("Cluster1",23),rep("Cluster2",21),rep("Cluster3",21),rep("Cluster4",3)),
                                   levels = c("Cluster1","Cluster2","Cluster3","Cluster4")),
                    #annotation_name_rot = 45,
                    col = list(Cluster=c("Cluster1"="#CB9C7A","Cluster2"="#8696a7",
                                         "Cluster3"="#CDB97D","Cluster4"="#7b8b6f")),show_legend = F)

#df_sub1=t(scale(t(df_sub)))"#A59B95"
#library()
Heatmap(data, cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        col=c("#87CEFA", "white", "#CC2121"),
        top_annotation = ha1,
        left_annotation = ha2,
        row_names_gp = gpar(fontsize = 9), 
        show_heatmap_legend = T)
ht_list=Heatmap(data, cluster_rows = F,
                cluster_columns = F,
                show_column_names = T,
                show_row_names = T,
                heatmap_legend_param = list(
                  title = "NES", at = c(-4, 0,4),
                  labels = c("-4","0", "4")
                ),
                
                col=c("#87CEFA", "white", "#CC2121"),
                top_annotation = ha1,
                left_annotation = ha2,
                column_names_gp  =gpar(
                  col = "black",
                  fontsize = 10
                ),
                row_names_gp = gpar(fontsize = 10),
                #column_names_rot =45,
                show_heatmap_legend = T)
draw(ht_list, heatmap_legend_side = "left")


############################14.4 all gene NES 热图 fgsea ######################################
load("/home/panxl/CRC2/human_examplePathways.RData")
library(fgsea)
library(ggplot2)
df=c("CMS1","CMS2","CMS3","CMS4","human")

pathway_cms_nes=lapply(df, function(x){
  y=read.csv(paste("/home/panxl/CRC2/",x,"_DEG.csv",sep=""))
  hom_gene=y
  cms_gene_df=as.character(hom_gene$X)
  cms_gene_df=bitr(cms_gene_df, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  gene_df <- left_join(y,cms_gene_df,c("X"="ENSEMBL"))#
  gene_df =na.omit(gene_df)
  gene_res= as.vector(gene_df$ENTREZID)#
  geneList<-gene_res
  data(geneList,package = "DOSE")
  geneList = gene_df$logFC
  names(geneList) = gene_df$ENTREZID
  geneList = sort(geneList,decreasing = T)
  fgseaRes <- fgsea(pathways = examplePathways, 
                    stats    = geneList,
                    minSize  = 15,
                    maxSize  = 500)
  pathway_nes=fgseaRes[,c(1,6)]
  pathway_nes$type=x
  return(pathway_nes)
})%>%dplyr::bind_rows()

load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")
cluster_pathway =cluster4_pathway.kegg
cluster_pathway=cluster_pathway[,c(1,2,6,10)]
cluster1_pathway=cluster_pathway[cluster_pathway$type=="Cluster 1"&cluster_pathway$p.adjust<=0.05,]#[c(1:25),]#26
cluster2_pathway=cluster_pathway[cluster_pathway$type=="Cluster 2"&cluster_pathway$p.adjust<=0.05,]#[18
cluster3_pathway=cluster_pathway[cluster_pathway$type=="Cluster 3"&cluster_pathway$p.adjust<=0.05,][c(1:25),]#[c(1:25),]#68  25
cluster4_pathway=cluster_pathway[cluster_pathway$type=="Cluster 4"&cluster_pathway$p.adjust<=0.05,]#7   10
intersect(cluster1_pathway$ID,cluster2_pathway$ID)
l12=c("hsa04060","hsa04151", "hsa04610")
cluster1_pathway=cluster1_pathway[!(cluster1_pathway$ID%in%cluster2_pathway$ID),]# 23"hsa04060"cluste2 --1 19
cluster2_pathway=cluster2_pathway[!(cluster2_pathway$ID%in%c("hsa04060","hsa04061")),]#23
intersect(cluster1_pathway$ID,cluster3_pathway$ID)
intersect(cluster1_pathway$ID,cluster4_pathway$ID)
intersect(cluster2_pathway$ID,cluster3_pathway$ID)
l23= "hsa04020"
cluster2_pathway=cluster2_pathway[!(cluster2_pathway$ID%in%c("hsa04020")),]#15
cluster3_pathway=cluster3_pathway[!(cluster3_pathway$ID%in%c("hsa04020")),]#24
intersect(cluster2_pathway$ID,cluster4_pathway$ID)
intersect(cluster3_pathway$ID,cluster4_pathway$ID)
l34=c("hsa04080" ,"hsa04911" ,"hsa04724")
cluster3_pathway=cluster3_pathway[!(cluster3_pathway$ID%in%l34),]#49
cluster4_pathway=cluster4_pathway[!(cluster4_pathway$ID%in%l34),]#5


pathway_cms_nes1=pathway_cms_nes
pathway_cms_nes1$cluster=case_when(
  pathway_cms_nes1$pathway%in%cluster1_pathway$ID~"Cluster1",
  pathway_cms_nes1$pathway%in%cluster2_pathway$ID~"Cluster2",
  pathway_cms_nes1$pathway%in%cluster3_pathway$ID~"Cluster3",
  pathway_cms_nes1$pathway%in%cluster4_pathway$ID~"Cluster4",
)#260
pathway_cms_nes1=na.omit(pathway_cms_nes1)

load("/home/panxl/CRC2/mouse_pathway_nes.RData")
human_nes=pathway_cms_nes1
nes=rbind(human_nes,pathway_nes)
table(nes$pathway)
NES=human_nes
NES=nes[!(nes$pathway%in%c("hsa04514","hsa00670","hsa05310","")),]
df=c("CMS1","CMS2","CMS3","CMS4","Human")
# which(NES$type=="CMS1")
CMS1=filter(NES,type=="CMS1")[,c(1,2)]
CMS2=filter(NES,type=="CMS2")[,c(1,2)]
CMS3=filter(NES,type=="CMS3")[,c(1,2)]
CMS4=filter(NES,type=="CMS4")[,c(1,2)]
human=filter(NES,type=="human")[,c(1,2,4)]
l=left_join(CMS1,CMS2,by = c("pathway"))
ll=left_join(l,CMS3,by = c("pathway"))
lll=left_join(ll,CMS4,by = c("pathway"))
llll=left_join(lll,human,by = c("pathway"))
colnames(llll)=c("Pathway","CMS1","CMS2","CMS3","CMS4","Human","Cluster")
load("/home/panxl/CRC2/cluster4_pathway.kegg.RData")
cluster4_pathway.kegg=cluster4_pathway.kegg[,c(1,2)]
NES3=left_join(llll,cluster4_pathway.kegg,c("Pathway"="ID"))
NES4=na.omit(NES3)
NES5=unique(NES4)
NES =NES5%>%arrange(Cluster)
rownames(NES)=NES$Description
colnames(NES)
data=NES[,c(2:6)]
rownames(data)=NES$Description
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

ha1 = HeatmapAnnotation(Group=factor(c( "Human",  "CMS1" , "CMS2" , "CMS3" , "CMS4" ),
                                     levels = c( "Human",  "CMS1" , "CMS2" , "CMS3" , "CMS4" )),
                        col = list(Group = c("Human"="#78a355",'CMS1' = '#b97573',
                                             'CMS2' = '#686954', 'CMS3' = '#5b5e66', 'CMS4'= '#816147')),show_legend = F)

table(NES$Cluster)
ha2 = rowAnnotation(Cluster=factor(c(rep("Cluster1",23),rep("Cluster2",22),rep("Cluster3",21),rep("Cluster4",3)),
                                   levels = c("Cluster1","Cluster2","Cluster3","Cluster4")),
                    #annotation_name_rot = 45,
                    col = list(Cluster=c("Cluster1"="#CB9C7A","Cluster2"="#8696a7",
                                         "Cluster3"="#CDB97D","Cluster4"="#7b8b6f")),show_legend = F)

ht_list=Heatmap(data, cluster_rows = F,
                cluster_columns = F,
                show_column_names = T,
                show_row_names = T,
                heatmap_legend_param = list(
                  title = "NES", at = c(-4, 0,4),
                  labels = c("-4","0", "4")
                ),
                
                col=c("#87CEFA", "white", "#CC2121"),
                top_annotation = ha1,
                left_annotation = ha2,
                column_names_gp  =gpar(
                  col = "black",
                  fontsize = 10
                ),
                row_names_gp = gpar(fontsize = 10),
                #column_names_rot =45,
                show_heatmap_legend = T)
draw(ht_list, heatmap_legend_side = "left")







######################################15.CMS correlation  pheatmap##############################################
######################15.1 cor_TCGA_CRC
kegg_annotation=read.csv("/home/panxl/CRC2/kegg_annotation.csv")
kegg.gene_symbol = as.character(unique(kegg_annotation$geneList))#2954
kegg.gene_df = bitr(kegg.gene_symbol,
                    fromType="SYMBOL",
                    toType="ENSEMBL",
                    OrgDb = "org.Hs.eg.db", drop = TRUE) #2888
kegg.gene_df=na.omit(kegg.gene_df)
gene_df <-left_join(kegg_annotation,kegg.gene_df,c("geneList"="SYMBOL"))#4192 71
#write.csv(gene_df,file = "/home/panxl/CRC2/unique_cluster_gene.csv",row.names = F)
gene_df=read.csv("/home/panxl/CRC2/unique_cluster_gene.csv")
table(gene_df$type)
H0M_logFC1=read.csv("/home/panxl/CRC2/all_cluater")
cor_CRC_kegg <- purrr::map(unique(gene_df$type),function(x){
  gene=gene_df$ENSEMBL[gene_df$type==x]
  ll=H0M_logFC1[H0M_logFC1$Gene.stable.ID%in%gene,]
  a = cor.test(ll$Mouse_logFC1,ll$Human_logFC1)
  #a$p.value
  re_df <- data.frame("pathway_id"=x,"TCGA_CRC"=a$estimate,"pvalue"=a$p.value)
  row.names(re_df)<-x
  return(re_df)
})dplyr::bind_rows()
write.csv(cor_CRC_kegg,file = "/home/panxl/CRC2/cor_CRC_pathway_log.csv",row.names = F)
######################15.2 cor_TCGA_CMS
gene_df=read.csv("/home/panxl/CRC2/unique_cluster_gene.csv")
H0M_logFC1=read.csv("/home/panxl/CRC2/CMS_all_cluater.csv")
H0M_logFC=H0M_logFC1[H0M_logFC1$type=="CMS3",]#每个亚型
cor_cms_CRC_kegg <- purrr::map(unique(gene_df$type),function(x){
  gene=gene_df$ENSEMBL[gene_df$type==x]#通路基因
  ll=H0M_logFC[H0M_logFC$Gene.stable.ID%in%gene,]
  a = cor.test(ll$Mouse_logFC1,ll$Human_logFC1)
  re_df <- data.frame("pathway_id"=x,"CMS3"=a$estimate,"pvalue"=a$p.value)
  row.names(re_df)<-x
  return(re_df)
})dplyr::bind_rows()
write.csv(cor_cms_CRC_kegg,file = "/home/panxl/CRC2/cor_CMS3_pathway_log.csv",row.names = F)
######################15.3 cor_TCGA 热图
TCGA_CRC_pathway=read.csv("/home/panxl/CRC2/cor_CRC_pathway_log.csv")[,c(1,2)]
TCGA_CMS3_pathway=read.csv("/home/panxl/CRC2/cor_CMS3_pathway_log.csv")[,c(1,2)]
TCGA_CMS1_pathway=read.csv("/home/panxl/CRC2/cor_CMS1_pathway_log.csv")[,c(1,2)]
TCGA_CMS2_pathway=read.csv("/home/panxl/CRC2/cor_CMS2_pathway_log.csv")[,c(1,2)]
TCGA_CMS4_pathway=read.csv("/home/panxl/CRC2/cor_CMS4_pathway_log.csv")[,c(1,2)]
a=merge(TCGA_CRC_pathway,TCGA_CMS1_pathway)
b=merge(a,TCGA_CMS2_pathway)
c=merge(b,TCGA_CMS3_pathway)
d=merge(c,TCGA_CMS4_pathway)
NES=read.csv("/home/panxl/CRC2/4cms_nes.csv",row.names = 1)
colnames(NES)
data=NES[,c(1,6,7)]
e=merge(d,data)arrange(cluster)#64
data=e[,c(2,3,4,5,6)]
colnames(data)=c("Human","CMS1","CMS2",
                 "CMS3","CMS4")
rownames(data)=e$Description
ha1 = HeatmapAnnotation(Group=factor(c("Human","CMS1","CMS2",
                                       "CMS3","CMS4"),
                                     levels =c("Human","CMS1","CMS2",
                                               "CMS3","CMS4")),
                        col = list(Group = c("Human"="#78a355","CMS1"="#abc88b","CMS2"="#a3cf62",
                                             "CMS3"="#b2d235","CMS4"="#bed742")),show_legend = F)#c("CY"=color_list1[4], "SG"="#65c294","DC"="#e0861a")))#,"CY"= color_list1[1], "SG" = color_list1[2],"DC"= color_list1[6]

table(e$cluster)
ha2 = rowAnnotation(cluster=factor(c(rep("Cluster1",24),rep("Cluster2",23),rep("Cluster3",21),rep("Cluster4",3)),
                                   levels = c("Cluster1","Cluster2","Cluster3","Cluster4")),
                    col = list(cluster=c("Cluster1"="#800000FF","Cluster2"="#767676FF","Cluster3"="#FFA319FF","Cluster4"="#8A9045FF")),show_legend = F)

Heatmap(data, cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        #show_column_names= F,
        # col=c("#87CEFA", "white", "#CC2121"),
        col=c("#4764A6", "white", "#aa2116"),
        top_annotation = ha1,
        left_annotation = ha2,
        row_names_gp = gpar(fontsize = 9), 
        show_heatmap_legend = F)#"navy","white","#ae545e"
################################## 16.cor #####################################################
###########################16.1 crc pathway --gene
library(clusterProfiler) 
library(org.Hs.eg.db)
library(KEGGREST)
load(file = "/home/panxl/CRC2/cluster4_pathway.kegg.RData")#1175
cluster_pathway=cluster4_pathway.kegg[cluster4_pathway.kegg$p.adjust<=0.05,]$ID#125
ll=c("hsa01523","hsa01212")
#ll=c("hsa04110", "hsa03030","hsa03008", "hsa03430","hsa03440", "hsa03460")# no gene pathway
pathway1=cluster_pathway[!(cluster_pathway%in%ll)]
#org=keggList("organism")
result=lapply(pathway1, function(x){
  gs<-keggGet(x)
  gene=unlist(lapply(gs[[1]]$GENE,function(l){strsplit(l,";")}))
  geneList=gene[1:length(gene)%%3==2]
  geneList=data.frame(geneList)
  geneList$type=x
  return(geneList)
})dplyr::bind_rows()#7594
#write.csv(result,file = "/home/panxl/CRC2/crc_cluster_kegg_annotation.csv",row.names = F)
##############################16.2 
kegg_annotation=read.csv( "/home/panxl/CRC2/crc_cluster_kegg_annotation.csv")
kegg.gene_symbol = as.character(unique(kegg_annotation$geneList))#2954
kegg.gene_df = bitr(kegg.gene_symbol,
                    fromType="SYMBOL",
                    toType="ENSEMBL",
                    OrgDb = "org.Hs.eg.db", drop = TRUE) #2888
kegg.gene_df=na.omit(kegg.gene_df)#4632
gene_df <-left_join(kegg_annotation,kegg.gene_df,c("geneList"="SYMBOL"))#17079 123  4192 71
#write.csv(gene_df,file = "/home/panxl/CRC2/crc_4cluster_gene.csv",row.names = F)
gene_df=read.csv("/home/panxl/CRC2/crc_4cluster_gene.csv")
table(gene_df$type)
H0M_logFC1=read.csv("/home/panxl/CRC2/all_cluater")
cor_CRC_kegg <- purrr::map(unique(gene_df$type),function(x){
  gene=gene_df$ENSEMBL[gene_df$type==x]
  ll=H0M_logFC1[H0M_logFC1$Gene.stable.ID%in%gene,]
  a = cor.test(ll$Mouse_logFC1,ll$Human_logFC1)
  #a$p.value
  re_df <- data.frame("pathway_id"=x,"TCGA_CRC"=a$estimate,"pvalue"=a$p.value)
  row.names(re_df)<-x
  return(re_df)
})dplyr::bind_rows()
#write.csv(cor_CRC_kegg,file = "/home/panxl/CRC2/cor_CRC_all_pathway.csv",row.names = F)
#############################16.3
cor_CRC= read.csv("/home/panxl/CRC2/cor_CRC_all_pathway.csv")
load(file = "/home/panxl/CRC2/cluster4_pathway.kegg.RData")#1175
cluster4_pathway.kegg=cluster4_pathway.kegg[cluster4_pathway.kegg$p.adjust<=0.05,]

cluster1_pathway=cor_CRC[cor_CRC$pathway_id%in%cluster4_pathway.kegg[cluster4_pathway.kegg$type=="Cluster 1: Up-regulated in both species",]$ID,]##[c(1:25),]#26
cluster1_pathway$Cluster="Cluster1"
cluster2_pathway=cor_CRC[cor_CRC$pathway_id%in%cluster4_pathway.kegg[cluster4_pathway.kegg$type=="Cluster 2: Up-regulated in mouse reverse in human",]$ID,]##[c(1:25),]#26
cluster2_pathway$Cluster="Cluster2"
cluster3_pathway=cor_CRC[cor_CRC$pathway_id%in%cluster4_pathway.kegg[cluster4_pathway.kegg$type=="Cluster 3: Down-regulated in both species",]$ID,]##[c(1:25),]#26
cluster3_pathway$Cluster="Cluster3"
cluster4_pathway=cor_CRC[cor_CRC$pathway_id%in%cluster4_pathway.kegg[cluster4_pathway.kegg$type=="Cluster 4: Up-regulated in human reverse in mouse",]$ID,]##[c(1:25),]#26
cluster4_pathway$Cluster="Cluster4"
c12=rbind(cluster1_pathway,cluster2_pathway)
c123=rbind(c12,cluster3_pathway)
c1234=rbind(c123,cluster4_pathway)[,c(2,4)]
colnames(c1234)=c("Cor","Cluster")
library(ggsci)
library(ggsignif)

ggplot(c1234,aes(x=Cluster,y=Cor))+
  geom_boxplot(aes(fill=Cluster))#+
  geom_signif(comparisons = list(c("Cluster1","Cluster2"),c("Cluster1","Cluster3")),
              y_position = c(0.6,3),#显著性比较
              map_signif_level = T)
########################17.1 GSEA  PI3K-Akt signaling pathway###########################
load("/home/panxl/CRC2/mouse_gse_kegg.Rdata")
load("/home/panxl/CRC2/human_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS1_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS2_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS3_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS4_gse_kegg.Rdata")


## 画GSEA图，加标签
library(circlize)
library(enrichplot )
library(ComplexHeatmap)
########Human
lgd1 = Legend(
  labels = c(paste("NES:",round(human_gse_kegg@result["hsa04151",]$NES,3)),
             paste("pvalue:",round(human_gse_kegg@result["hsa04151",]$pvalue,3)),
             paste("padj:",round(human_gse_kegg@result["hsa04151",]$p.adjust,3))),
   row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)
 NM''
p11<-gseaplot2(human_gse_kegg, 
              geneSetID = 1, 
              base_size=20,  ## 字体大小
              title = paste("Human",human_gse_kegg@result["hsa04151",]$Description))

grid.newpage()
pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
grid.rect()

p11  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例
########CMS1
## width=1080&height=754   10.8    7.54
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS1_gse_kegg@result["hsa04151",]$NES,3)),
             paste("pvalue:",round(CMS1_gse_kegg@result["hsa04151",]$pvalue,3)),
             paste("padj:",round(CMS1_gse_kegg@result["hsa04151",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p12=gseaplot2(
  CMS1_gse_kegg,
  "hsa04151",
  title = "CMS1 PI3K-Akt signaling pathway",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p12  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例



########CMS2
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS2_gse_kegg@result["hsa04151",]$NES,3)),
             paste("pvalue:",round(CMS2_gse_kegg@result["hsa04151",]$pvalue,3)),
             paste("padj:",round(CMS2_gse_kegg@result["hsa04151",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p13=gseaplot2(
  CMS2_gse_kegg,
  "hsa04151",
  title = "CMS2 PI3K-Akt signaling pathway",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p13  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例

########CMS3
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS3_gse_kegg@result["hsa04151",]$NES,3)),
             paste("pvalue:",round(CMS3_gse_kegg@result["hsa04151",]$pvalue,3)),
             paste("padj:",round(CMS3_gse_kegg@result["hsa04151",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p14=gseaplot2(
  CMS3_gse_kegg,
  "hsa04151",
  title = "CMS3 PI3K-Akt signaling pathway",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p14  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例

########CMS4
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS4_gse_kegg@result["hsa04151",]$NES,3)),
             paste("pvalue:",round(CMS4_gse_kegg@result["hsa04151",]$pvalue,3)),
             paste("padj:",round(CMS4_gse_kegg@result["hsa04151",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p15=gseaplot2(
  CMS4_gse_kegg,
  "hsa04151",
  title = "CMS4 PI3K-Akt signaling pathway",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p15  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例
#####mouse
lgd1 = Legend(
  labels = c(paste("NES:",round(mouse_gse_kegg@result["mmu04151",]$NES,3)),
             paste("pvalue:",round(mouse_gse_kegg@result["mmu04151",]$pvalue,3)),
             paste("padj:",round(mouse_gse_kegg@result["mmu04151",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p16=gseaplot2(
  mouse_gse_kegg,
  "mmu04151",
  title = "Mouse PI3K-Akt signaling pathway",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p16  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例


########################17.2 GSEA  JAK-STAT signal pathway###########################
load("/home/panxl/CRC2/mouse_gse_kegg.Rdata")
load("/home/panxl/CRC2/human_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS1_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS2_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS3_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS4_gse_kegg.Rdata")


## 画GSEA图，加标签
library(circlize)
library(enrichplot )
library(ComplexHeatmap)
########Human
lgd1 = Legend(
  labels = c(paste("NES:",round(human_gse_kegg@result["hsa04630",]$NES,3)),
             paste("pvalue:",round(human_gse_kegg@result["hsa04630",]$pvalue,3)),
             paste("padj:",round(human_gse_kegg@result["hsa04630",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p11<-gseaplot2(human_gse_kegg, 
               geneSetID = 1, 
               base_size=20,  ## 字体大小
               title = paste("Human",human_gse_kegg@result["hsa04630",]$Description))

grid.newpage()
pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
grid.rect()

p11  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例
########CMS1
## width=1080&height=754   10.8    7.54
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS1_gse_kegg@result["hsa04630",]$NES,3)),
             paste("pvalue:",round(CMS1_gse_kegg@result["hsa04630",]$pvalue,3)),
             paste("padj:",round(CMS1_gse_kegg@result["hsa04630",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p12=gseaplot2(
  CMS1_gse_kegg,
  "hsa04630",
  title = "CMS1 JAK-STAT signal pathway",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p12  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例



########CMS2
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS2_gse_kegg@result["hsa04630",]$NES,3)),
             paste("pvalue:",round(CMS2_gse_kegg@result["hsa04630",]$pvalue,3)),
             paste("padj:",round(CMS2_gse_kegg@result["hsa04630",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p13=gseaplot2(
  CMS2_gse_kegg,
  "hsa04630",
  title = "CMS2 JAK-STAT signal pathway",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p13  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例

########CMS3
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS3_gse_kegg@result["hsa04630",]$NES,3)),
             paste("pvalue:",round(CMS3_gse_kegg@result["hsa04630",]$pvalue,3)),
             paste("padj:",round(CMS3_gse_kegg@result["hsa04630",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p14=gseaplot2(
  CMS3_gse_kegg,
  "hsa04630",
  title = "CMS3 JAK-STAT signal pathway",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p14  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例

########CMS4
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS4_gse_kegg@result["hsa04630",]$NES,3)),
             paste("pvalue:",round(CMS4_gse_kegg@result["hsa04630",]$pvalue,3)),
             paste("padj:",round(CMS4_gse_kegg@result["hsa04630",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p15=gseaplot2(
  CMS4_gse_kegg,
  "hsa04630",
  title = "CMS4 JAK-STAT signal pathway",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p15  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例
#####mouse
lgd1 = Legend(
  labels = c(paste("NES:",round(mouse_gse_kegg@result["mmu04630",]$NES,3)),
             paste("pvalue:",round(mouse_gse_kegg@result["mmu04630",]$pvalue,3)),
             paste("padj:",round(mouse_gse_kegg@result["mmu04630",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p16=gseaplot2(
  mouse_gse_kegg,
  "mmu04630",
  title = "Mouse JAK-STAT signal pathway",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p16  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例
########################17.3 GSEA  Th17 cell differential###########################
load("/home/panxl/CRC2/mouse_gse_kegg.Rdata")
load("/home/panxl/CRC2/human_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS1_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS2_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS3_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS4_gse_kegg.Rdata")


## 画GSEA图，加标签
library(circlize)
library(enrichplot )
library(ComplexHeatmap)
########Human
lgd1 = Legend(
  labels = c(paste("NES:",round(human_gse_kegg@result["hsa04659",]$NES,3)),
             paste("pvalue:",round(human_gse_kegg@result["hsa04659",]$pvalue,3)),
             paste("padj:",round(human_gse_kegg@result["hsa04659",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p11<-gseaplot2(human_gse_kegg, 
               geneSetID = 1, 
               base_size=20,  ## 字体大小
               title = paste("Human",human_gse_kegg@result["hsa04659",]$Description))

grid.newpage()
pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
grid.rect()

p11  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例
########CMS1
## width=1080&height=754   10.8    7.54
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS1_gse_kegg@result["hsa04659",]$NES,3)),
             paste("pvalue:",round(CMS1_gse_kegg@result["hsa04659",]$pvalue,3)),
             paste("padj:",round(CMS1_gse_kegg@result["hsa04659",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p12=gseaplot2(
  CMS1_gse_kegg,
  "hsa04659",
  title = "CMS1 Th17 cell differential",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p12  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例



########CMS2
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS2_gse_kegg@result["hsa04659",]$NES,3)),
             paste("pvalue:",round(CMS2_gse_kegg@result["hsa04659",]$pvalue,3)),
             paste("padj:",round(CMS2_gse_kegg@result["hsa04659",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p13=gseaplot2(
  CMS2_gse_kegg,
  "hsa04659",
  title = "CMS2 Th17 cell differential",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p13  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例

########CMS3
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS3_gse_kegg@result["hsa04659",]$NES,3)),
             paste("pvalue:",round(CMS3_gse_kegg@result["hsa04659",]$pvalue,3)),
             paste("padj:",round(CMS3_gse_kegg@result["hsa04659",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p14=gseaplot2(
  CMS3_gse_kegg,
  "hsa04659",
  title = "CMS3 Th17 cell differential",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p14  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例

########CMS4
lgd1 = Legend(
  labels = c(paste("NES:",round(CMS4_gse_kegg@result["hsa04659",]$NES,3)),
             paste("pvalue:",round(CMS4_gse_kegg@result["hsa04659",]$pvalue,3)),
             paste("padj:",round(CMS4_gse_kegg@result["hsa04659",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p15=gseaplot2(
  CMS4_gse_kegg,
  "hsa04659",
  title = "CMS4 Th17 cell differential",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p15  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例
#####mouse
lgd1 = Legend(
  labels = c(paste("NES:",round(mouse_gse_kegg@result["mmu04659",]$NES,3)),
             paste("pvalue:",round(mouse_gse_kegg@result["mmu04659",]$pvalue,3)),
             paste("padj:",round(mouse_gse_kegg@result["mmu04659",]$p.adjust,3))),
  row_gap = unit(1.5, "mm"),  ## 控制文字行与行之间的间隔
  ## 控制文字字体大小
  labels_gp = gpar(fontsize = 20)  )

pd = packLegend(lgd1)

p16=gseaplot2(
  mouse_gse_kegg,
  "mmu04659",
  title = "Mouse Th17 cell differential",
  color = "green",
  base_size = 19,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)
p16  ## 画热图
draw(pd, x = unit(0.27, "npc"), y = unit(0.65, "npc"))   # 画图例

#########################18.cluster3(5)-intestine(species) specific gene pathway######
############18.1 cluster3-intestine specific gene pathway#######################
df_intestine = read.csv("/home/shimw/project/enhancer_map/conserved/colon-specific.csv")
H0M_logFC1= read.csv("/home/panxl/CRC2/all_cluater")
H0M_logFC1$cluster=gsub("Cluster 1: Up-regulated in both species","Cluster 1",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 2: Up-regulated in mouse reverse in human","Cluster 2",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 3: Down-regulated in both species","Cluster 3",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 4: Up-regulated in human reverse in mouse","Cluster 4",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 5: Not-significant in both species","Cluster 5",H0M_logFC1$cluster)
df_data=H0M_logFC1[H0M_logFC1$cluster=="Cluster 3",]
##mouse
intestine_gene=df_intestine[df_intestine$Gene%in%df_data$Mouse.gene.stable.ID,]$Gene
library(clusterProfiler) 
library(org.Mm.eg.db)
library(patchwork)
gene_res_ENSEMBL = intestine_gene
gene_res_df = bitr(gene_res_ENSEMBL, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Mm.eg.db", drop = TRUE)
gene_df <- left_join(df_data,gene_res_df,c("Mouse.gene.stable.ID"="ENSEMBL"))#
gene_df1=na.omit(gene_df)
mouse_intestine_kegg <- enrichKEGG(gene= gene_df1$ENTREZID,
                                   keyType       = "kegg",
                                   organism      = "mmu",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05)

mm=mouse_intestine_kegg@result
#save(mouse_intestine_kegg,file = "/home/panxl/CRC2/mouse_intestine_kegg.Rdata")

# gene_res0.01= as.vector(gene_df1$ENTREZID)#
# geneList<-gene_res0.01
# data(geneList,package = "DOSE")
# geneList = gene_df$Mouse_logFC1
# names(geneList) = gene_df$ENTREZID
# geneList = sort(geneList,decreasing = T)
# geo_list<- list()
# pvalue_Cutoff=0.05
# mouse_intestine_gse_kegg<-gseKEGG(geneList= geneList,
#                         organism = "mmu",
#                         keyType  = 'kegg',
#                         nPerm  = 1000,
#                         minGSSize = 10,
#                         maxGSSize = 500,
#                         pvalueCutoff = 0.05,
#                         pAdjustMethod     = "BH") ## 0.05有的没有结果 
# #save(mouse_intestine_gse_kegg,file = "/home/panxl/CRC2/mouse_intestine_gse_kegg.Rdata")

####human
intestine_gene=df_data[df_data$Mouse.gene.stable.ID%in%df_intestine$Gene,]$Gene.stable.ID
library(clusterProfiler) 
library(org.Hs.eg.db)
library(patchwork)
gene_res_ENSEMBL = intestine_gene
gene_res_df = bitr(gene_res_ENSEMBL, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_df <- left_join(df_data,gene_res_df,c("Gene.stable.ID"="ENSEMBL"))#
gene_df1=na.omit(gene_df)
human_intestine_kegg <- enrichKEGG(gene= gene_df1$ENTREZID,
                      keyType       = "kegg",
                      organism      = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)

hh=human_intestine_kegg@result
#save(human_intestine_kegg,file = "/home/panxl/CRC2/human_intestine_kegg.Rdata")
# gene_res0.01= as.vector(gene_df1$ENTREZID)#
# geneList<-gene_res0.01
# data(geneList,package = "DOSE")
# geneList = gene_df$Human_logFC1
# names(geneList) = gene_df$ENTREZID
# geneList = sort(geneList,decreasing = T)
# geo_list<- list()
# pvalue_Cutoff=0.05
# human_intestine_gse_kegg<-gseKEGG(geneList= geneList,
#                                   organism = "hsa",
#                                   keyType  = 'kegg',
#                                   nPerm  = 1000,
#                                   minGSSize = 10,
#                                   maxGSSize = 500,
#                                   pvalueCutoff = 0.05,
#                                   pAdjustMethod     = "BH") ## 0.05有的没有结果 
# #save(human_intestine_gse_kegg,file = "/home/panxl/CRC2/human_intestine_gse_kegg.Rdata")
################### intestine pathway
load("/home/panxl/CRC2/mouse_intestine_kegg.Rdata")
load("/home/panxl/CRC2/human_intestine_kegg.Rdata")
mouse_intestine=mouse_intestine_kegg@result#125
mouse_intestine1=mouse_intestine[mouse_intestine$pvalue<=0.05,]#12 2
mouse_intestine1$Description
human_intestine=human_intestine_kegg@result#126
human_intestine1=human_intestine[human_intestine$pvalue<=0.05,]#12 2
human_intestine1$Description
mouse_intestine1=mouse_intestine1mutate(Description = fct_reorder(Description, -pvalue))
ggplot(data=mouse_intestine1,aes(x=Description,y=-log10(pvalue),fill=p.adjust))+#fill=padj fill颜色填充，使用连续值padj 
  geom_bar(stat="identity") + coord_flip()+#coord_flip()颠倒坐标轴
  labs(x="",y="-log10 (P_value)",title="KEGG Intestine Pathway")+
  scale_fill_gradient(low="#DC143C",high="#8696a7")+
  theme_classic()+
  theme(axis.title.x = element_text(size =22, 
                                    color = "black"),
        #face = "italic"),
        axis.title.y = element_text(size = 25,
                                    color = "black",
                                    #face = "italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        title = element_text(size = 22,
                             color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 17 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   #face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 19,  
                                   color = "black",
                                   vjust = 0.4))

human_intestine1=human_intestine1mutate(Description = fct_reorder(Description, -pvalue))

ggplot(data=human_intestine1,aes(x=Description,y=-log10(pvalue),fill=p.adjust))+#fill=padj fill颜色填充，使用连续值padj 
  geom_bar(stat="identity") + coord_flip()+#coord_flip()颠倒坐标轴
  labs(x="",y="-log10 (P_value)",title="KEGG Intestine Pathway")+
  scale_fill_gradient(low="#DC143C",high="#8696a7")+
  theme_classic()+
  theme(axis.title.x = element_text(size =22, 
                                    color = "black"),
        #face = "italic"),
        axis.title.y = element_text(size = 25,
                                    color = "black",
                                    #face = "italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        title = element_text(size = 22,
                             color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 17 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   #face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 19,  
                                   color = "black",
                                   vjust = 0.4))
############18.2 cluster5-species specific gene pathway###################
df_species = read.csv("/home/shimw/project/enhancer_map/conserved/species-specific.csv")#744
ensembl_ID_change=read.csv("/home/panxl/ensembl_37gene.txt")#28255
H0M_logFC1= read.csv("/home/panxl/CRC2/all_cluater")
H0M_logFC1$cluster=gsub("Cluster 1: Up-regulated in both species","Cluster 1",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 2: Up-regulated in mouse reverse in human","Cluster 2",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 3: Down-regulated in both species","Cluster 3",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 4: Up-regulated in human reverse in mouse","Cluster 4",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 5: Not-significant in both species","Cluster 5",H0M_logFC1$cluster)
H0M_logFC2=left_join(H0M_logFC1,ensembl_ID_change,c("Gene.stable.ID"="Gene.stable.ID"))
H0M_logFC3=na.omit(H0M_logFC2)
df_data=H0M_logFC3[H0M_logFC3$cluster=="Cluster 5",]
species_gene=df_data[df_data$Gene.name%in%df_species$gene_id,]
##mouse
library(clusterProfiler) 
library(org.Mm.eg.db)
library(patchwork)
gene_res_ENSEMBL = species_gene$Mouse.gene.stable.ID
gene_res_df = bitr(gene_res_ENSEMBL, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Mm.eg.db", drop = TRUE)
gene_df <- left_join(species_gene,gene_res_df,c("Mouse.gene.stable.ID"="ENSEMBL"))#
gene_df1=na.omit(gene_df)
mouse_species_kegg <- enrichKEGG(gene= gene_df1$ENTREZID,
                                   keyType       = "kegg",
                                   organism      = "mmu",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05)

mm=mouse_species_kegg@result
#save(mouse_species_kegg,file = "/home/panxl/CRC2/mouse_species_kegg.Rdata")

# gene_res0.01= as.vector(gene_df1$ENTREZID)#
# geneList<-gene_res0.01
# data(geneList,package = "DOSE")
# geneList = gene_df$Mouse_logFC1
# names(geneList) = gene_df$ENTREZID
# geneList = sort(geneList,decreasing = T)
# geo_list<- list()
# pvalue_Cutoff=0.05
# mouse_species_gse_kegg<-gseKEGG(geneList= geneList,
#                         organism = "mmu",
#                         keyType  = 'kegg',
#                         nPerm  = 1000,
#                         minGSSize = 10,
#                         maxGSSize = 500,
#                         pvalueCutoff = 0.05,
#                         pAdjustMethod     = "BH") ## 0.05有的没有结果
# #save(mouse_species_gse_kegg,file = "/home/panxl/CRC2/mouse_species_gse_kegg.Rdata")

####human
library(clusterProfiler) 
library(org.Hs.eg.db)
library(patchwork)
gene_res_ENSEMBL = species_gene$Gene.stable.ID
gene_res_df = bitr(gene_res_ENSEMBL, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_df <- left_join(df_data,gene_res_df,c("Gene.stable.ID"="ENSEMBL"))#
gene_df1=na.omit(gene_df)
human_species_kegg <- enrichKEGG(gene= gene_df1$ENTREZID,
                                   keyType       = "kegg",
                                   organism      = "hsa",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05)

hh=human_species_kegg@result
#save(human_species_kegg,file = "/home/panxl/CRC2/human_species_kegg.Rdata")
gene_res0.01= as.vector(gene_df1$ENTREZID)#
geneList<-gene_res0.01
data(geneList,package = "DOSE")
geneList = gene_df$Human_logFC1
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList,decreasing = T)
geo_list<- list()
pvalue_Cutoff=0.05
human_species_gse_kegg<-gseKEGG(geneList= geneList,
                                  organism = "hsa",
                                  keyType  = 'kegg',
                                  nPerm  = 1000,
                                  minGSSize = 10,
                                  maxGSSize = 500,
                                  pvalueCutoff = 0.05,
                                  pAdjustMethod     = "BH") ## 0.05有的没有结果
# #save(human_species_gse_kegg,file = "/home/panxl/CRC2/human_species_gse_kegg.Rdata")
################### intestine pathway
load("/home/panxl/CRC2/mouse_species_kegg.Rdata")
load("/home/panxl/CRC2/human_species_kegg.Rdata")
mouse_species=mouse_species_kegg@result
mouse_species1=mouse_species[mouse_species$p.adjust<=0.05,]#15
human_species=human_species_kegg@result
human_species1=human_species[human_species$p.adjust<=0.05,]#14
mouse_species1$Description
human_species1$Description
mouse_species1=mouse_species1%>%mutate(Description = fct_reorder(Description, -p.adjust))
ggplot(data=mouse_species1,aes(x=Description,y=-log10(p.adjust),fill=p.adjust))+#fill=padj fill颜色填充，使用连续值padj 
  geom_bar(stat="identity") + coord_flip()+#coord_flip()颠倒坐标轴
  labs(x="",y="-log10 (P_adjust)",title="KEGG Species Pathway")+
  scale_fill_gradient(low="#CB9C7A",high="#8696a7")+
  theme_classic()+
  theme(axis.title.x = element_text(size =22, 
                                    color = "black"),
        #face = "italic"),
        axis.title.y = element_text(size = 25,
                                    color = "black",
                                    #face = "italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        title = element_text(size = 22,
                             color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 17 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   #face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 19,  
                                   color = "black",
                                   vjust = 0.4))

human_species1=human_species1%>%mutate(Description = fct_reorder(Description, -p.adjust))

ggplot(data=human_species1,aes(x=Description,y=-log10(p.adjust),fill=p.adjust))+#fill=padj fill颜色填充，使用连续值padj 
  geom_bar(stat="identity") + coord_flip()+#coord_flip()颠倒坐标轴
  labs(x="",y="-log10 (P_adjust)",title="KEGG Species Pathway")+
  scale_fill_gradient(low="#CB9C7A",high="#8696a7")+
  theme_classic()+
  theme(axis.title.x = element_text(size =22, 
                                    color = "black"),
        #face = "italic"),
        axis.title.y = element_text(size = 25,
                                    color = "black",
                                    #face = "italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        title = element_text(size = 22,
                             color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 17 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   #face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 19,  
                                   color = "black",
                                   vjust = 0.4))

############################## 19 反卷积细胞成分#######################################################
# install.packages("remotes")
# library(remotes)
# local({r <- getOption("repos")  
# r["CRAN"] <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
# options(repos=r)}) 
# remotes::install_github("icbi-lab/")
library(ggplot2)
library(immunedeconv)
library(tidyverse)
###############################19.1 human and mouse TPM#######################################
#########################mouse tpm
# 读入表达数据
library(GenomicFeatures)
## Calculate gene length for 
txdb <- makeTxDbFromGFF("//home/panxl/reference/mouse_reference/gencode.vM20.annotation.gtf",format="gtf")
## 计算基因长度
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
names(exons_gene_lens) = gsub("\\..*","",names(exons_gene_lens))
len_data<-as.data.frame(exons_gene_lens)
exons_gene_lens=unlist(exons_gene_lens)
len_data<-as.data.frame(t(len_data))
len_data$gene_id<-rownames(len_data)
colnames(len_data)[1]<-"length"
#write.csv(len_data,"/home/panxl/CRC2/mouse_gene_length.csv", quote = F)
len_data <- read.csv("/home/panxl/CRC2/mouse_gene_length.csv",header = T,row.names = 1, check.names = F)
exons_gene_lens <- len_data$length
names(exons_gene_lens) <- len_data$gene_id
tpm_calc <- function(countdf, filename){
  kb <- exons_gene_lens / 1000
  #table(rownames(countdf)%in%names(kb))
  #rownames(countdf)[rownames(countdf)%in%names(kb)]
  countdf = countdf[rownames(countdf)%in%names(kb),]
  #kb[match(rownames(countdf),names(kb))]
  rpk <- countdf / kb[match(rownames(countdf),names(kb))]
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  write.csv(tpm, paste("/home/panxl/CRC2/","tpm-",filename,".csv",sep = ""), quote = F)
}
mouse_readcounts=read.csv("/home/panxl/CRC2/mouse78_count.csv")
rownames(mouse_readcounts)=mouse_readcounts$X
mouse_info <- read.csv("/home/panxl/CRC2/filter_mouse_info.csv",row.names = 1, check.names = F)
mouse_colon_tumor=mouse_readcounts[,colnames(mouse_readcounts)%in%mouse_info[mouse_info$Group=="Tumor",]$SRR_id]#44
mouse_colon_normal=mouse_readcounts[,colnames(mouse_readcounts)%in%mouse_info[mouse_info$Group=="Normal",]$SRR_id]#34
tpm_calc(mouse_colon_tumor, "mouse_colon_tumor")
tpm_calc(mouse_colon_normal, "mouse_colon_normal")
########################human tpm
####human
library(GenomicFeatures)
## 载入参考文件
txdb <- makeTxDbFromGFF("/home/shimw/project/cell_line/data/Homo_sapiens.GRCh38.104.gtf",format="gtf")
## 计算基因长度
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
len_data<-as.data.frame(exons_gene_lens)
exons_gene_lens=unlist(exons_gene_lens)
len_data<-as.data.frame(t(len_data))
len_data$gene_id<-rownames(len_data)
colnames(len_data)[1]<-"length"
#write.csv(len_data,"/home/panxl/CRC2/human_gene_length.csv", quote = F)
# genes <- stringr::str_split_fixed(len_data$gene_id, "\\.", 2)[,1]
# len_data$gene_id <- genes
# len_data <- len_data[!duplicated(genes), ]
len_data <- read.csv("/home/panxl/CRC2/human_gene_length.csv",header = T,row.names = 1, check.names = F)
exons_gene_lens <- len_data$length
names(exons_gene_lens) <- len_data$gene_id
tpm_calc <- function(countdf, filename){
  kb <- exons_gene_lens / 1000
  #table(rownames(countdf)%in%names(kb))
  #rownames(countdf)[rownames(countdf)%in%names(kb)]
  countdf = countdf[rownames(countdf)%in%names(kb),]
  #kb[match(rownames(countdf),names(kb))]
  rpk <- countdf / kb[match(rownames(countdf),names(kb))]
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  write.csv(tpm, paste("/home/panxl/CRC2/","tpm-",filename,".csv",sep = "", row.names=F), quote = F)
}

human_readcounts=read.csv("/home/panxl/CRC2/human439_count.csv")
colnames(human_readcounts)=gsub("[.]","-",colnames(human_readcounts))
human_info <- read.csv("/home/panxl/CRC2/filter_human_info.csv",row.names = 1, check.names = F)
human_colon_tumor=human_readcounts[,colnames(human_readcounts)%in%human_info[human_info$Group=="Tumor",]$sample]#335
human_colon_normal=human_readcounts[,colnames(human_readcounts)%in%human_info[human_info$Group=="Normal",]$sample]#104
tpm_calc(human_colon_tumor, "human_colon_tumor")
tpm_calc(human_colon_normal, "human_colon_normal")
#####################19.2 反卷积解析细胞成分#################################
############同源基因
hom_gene=read.csv("/home/panxl/CRC2/all_cluater")
hom_df = bitr (hom_gene$Gene.stable.ID, 
                       fromType="ENSEMBL", 
                       toType="SYMBOL", 
                       OrgDb = "org.Hs.eg.db", drop = TRUE)
hom_RNASeq <- left_join(hom_gene,hom_df ,c("Gene.stable.ID"="ENSEMBL"))#176
hom_RNASeq1 =na.omit(hom_RNASeq)
hom_RNASeq2=hom_RNASeq1[!duplicated(hom_RNASeq1$Gene.stable.ID),]
hom_RNASeq3=hom_RNASeq2[!duplicated(hom_RNASeq2$SYMBOL),]
rownames(hom_RNASeq3)=hom_RNASeq3$SYMBOL
colnames(hom_RNASeq3)
#write.csv(hom_RNASeq3,file = "/home/panxl/CRC2/all_cluater_symbol.csv")

########## mouse 表达矩阵
hom_symbol = read.csv("/home/panxl/CRC2/all_cluater_symbol.csv")
mouse_tpm=read.csv("~/CRC2/tpm-mouse_colon_tumor.csv",row.names = 1)
mouse_tpm$ENSEMBL=rownames(mouse_tpm)
mouse_RNASeq <- left_join(mouse_tpm,hom_symbol ,c("ENSEMBL"="Mouse.gene.stable.ID"))#176
mouse_RNASeq1 =na.omit(mouse_RNASeq)
rownames(mouse_RNASeq1)=mouse_RNASeq1$SYMBOL
colnames(mouse_RNASeq2)
mouse_RNASeq2=mouse_RNASeq1[,-c(45:51)]#tumor 45,46 normal 35,36
#write.table(mouse_RNASeq2, file = "/home/panxl/CRC2/deconvolute/Mouse.txt",sep="\t")
#write.csv(mouse_RNASeq2, file = "/home/panxl/CRC2/mouse_deconvolute.csv")
colnames(mouse_RNASeq2)
mouse_RNASeq2[1:5,]
###### mouse deconvolute
# library(immunedeconv)
# #lim=read.table("/home/panxl/CRC2/deconvolute/LM22.txt",sep = "\t") 
# #ll=intersect(rownames(mouse_RNASeq2),lim$V1)#547  358
# mouse_deconvolute=read.csv("/home/panxl/CRC2/mouse_deconvolute.csv")
# set_cibersort_binary("/home/panxl/CRC2/deconvolute/CIBERSORT.R")
# set_cibersort_mat("/home/panxl/CRC2/deconvolute/LM22.txt")
# res <- deconvolute(mouse_deconvolute, method="cibersort")
setwd("/home/panxl/CRC2/deconvolute")
source('CIBERSORT.R')
#Mouse.txt=read.table("/home/panxl/CRC2/deconvolute/Mouse.txt")
Mouse_result <- CIBERSORT('LM22.txt','Mouse.txt', perm = 1000, QN = T)  #perm置换次数=1000，QN分位数归一化=TRUE
#############human 表达矩阵
hom_symbol = read.csv("/home/panxl/CRC2/all_cluater_symbol.csv")
human_tpm=read.csv("~/CRC2/tpm-human_colon_tumor.csv",row.names = 1)
human_tpm$ENSEMBL=rownames(human_tpm)
human_RNASeq <- left_join(human_tpm,hom_symbol ,c("ENSEMBL"="Gene.stable.ID"))#176
human_RNASeq1 =na.omit(human_RNASeq)
rownames(human_RNASeq1)=human_RNASeq1$SYMBOL
colnames(human_RNASeq2)
human_RNASeq2=human_RNASeq1[,-c(336:342)]#tumor 45,46 normal 35,36
#write.csv(human_RNASeq2, file = "/home/panxl/CRC2/human_deconvolute.csv")
#write.table(human_RNASeq2, file = "/home/panxl/CRC2/deconvolute/Human.txt",sep = "\t")
##########human  deconvolute
# library(immunedeconv)
# human=read.csv("/home/panxl/CRC2/human_deconvolute.csv")
# lim=read.table("/home/panxl/CRC2/deconvolute/LM22.txt",sep = "\t") 
# ll=intersect(human$X,lim$V1)#547  358
# human_deconvolute=read.csv("/home/panxl/CRC2/human_deconvolute.csv")
# set_cibersort_binary("/home/panxl/CRC2/deconvolute/CIBERSORT.R")
# set_cibersort_mat("/home/panxl/CRC2/deconvolute/LM22.txt")
# res <- deconvolute(human_deconvolute, method="cibersort")
setwd("/home/panxl/CRC2/deconvolute")
source('CIBERSORT.R')
#human.txt=read.table("/home/panxl/CRC2/deconvolute/Human.txt")
Human_result <- CIBERSORT('LM22.txt','Human.txt', perm = 1000, QN = T)  #perm置换次数=1000，QN分位数归一化=TRUE
#在同一文件夹下可以得到运算结果（"CIBERSORT-Results.txt"）
#注意Cibersort结果的默认文件名为CIBERSORT-Results.txt，在同一文件夹下进行第二次运算会覆盖第一次得到的文件，建议在每一次运算之后对文件重命名。

# ##########human
# deconvolute_human <- function(tpm, filename){
#   human_tpm=read.csv(paste("/home/panxl/CRC2/",tpm,".csv",sep=""),row.names = 1)
#   #human_tpm <- read.csv("~/CRC2/tpm-human_colon_tumor.csv",row.names = 1)
#   human_tpm$ENSEMBL=rownames(human_tpm)
#   human_RNASeq <- left_join(human_tpm,hom_RNASeq3 ,c("ENSEMBL"="Gene.stable.ID"))#176
#   human_RNASeq1 =na.omit(human_RNASeq)
#   rownames(human_RNASeq1)=human_RNASeq1$SYMBOL
#   colnames(human_RNASeq1)
#   human_RNASeq2=human_RNASeq1[,-c(336:341)]#tumor 45,46 normal 35,36
#   human_res <- deconvolute(human_RNASeq2, method="quantiseq")
#   write.csv(human_RNASeq2, paste("/home/panxl/CRC2/","deconvolute-",filename,".csv",sep = ""), quote = F)
#   
# }
# deconvolute_human("tpm-human_colon_tumor","human-tumor")
#############################19.3.反卷积画图 #######################################################
deconvolute_mouse=read.table("/home/panxl/CRC2/deconvolute/mouse_CIBERSORT-Results.txt",sep="\t")
colnames(deconvolute_mouse)=deconvolute_mouse[1,]
deconvolute_mouse=deconvolute_mouse[-1,-c(24,25,26)]
rownames(deconvolute_mouse)=deconvolute_mouse[,1]
#deconvolute_mouse=deconvolute_mouse[-11,]
list=c(2:23)
deconvolute_mouse1=lapply(list, function(x){
  ll=deconvolute_mouse[,c(1,x)]
  ll$type=colnames(ll)[2]
  ff=data.frame(Cell_type=ll[,3],Number=ll[,2],Sample=ll[,1])
  return(ff)
})%>%dplyr::bind_rows()
deconvolute_mouse1$Number=as.numeric(deconvolute_mouse1$Number)
#save(deconvolute_mouse1,file = "/home/panxl/CRC2/deconvolute_Mouse.RData")
ggplot(data=deconvolute_mouse1,aes(x=Sample,y=Number))+
  #geom_bar(position = pos)+
  geom_col(color="white",aes(fill=Cell_type),position="fill")+
  scale_fill_igv()+
    # labs(x="",y="",title="Cluster Distribution in CMS")+
  theme_classic()+
  theme(title = element_text(size = 20,
                             color = "black"),
        
        plot.title = element_text(hjust = 0.5),
        
        axis.text.x = element_text(size = 15 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 90), #角度
        axis.text.y = element_text(size = 20,  
                                   color = "black",
                                   vjust = 0.4))#+
theme_void()

###human
deconvolute_human=read.table("/home/panxl/CRC2/deconvolute/human_CIBERSORT-Results.txt",sep="\t")
colnames(deconvolute_human)=deconvolute_human[1,]
deconvolute_human=deconvolute_human[-1,-c(24,25,26)]
rownames(deconvolute_human)=deconvolute_human[,1]

list=c(2:23)
deconvolute_human1=lapply(list, function(x){
  ll=deconvolute_human[,c(1,x)]
  ll$type=colnames(ll)[2]
  ff=data.frame(Cell_type=ll[,3],Number=ll[,2],Sample=ll[,1])
  return(ff)
})%>%dplyr::bind_rows()
deconvolute_human1$Number=as.numeric(deconvolute_human1$Number)
#save(deconvolute_human1,file = "/home/panxl/CRC2/deconvolute_Human.RData")
ggplot(data=deconvolute_human1,aes(x=Sample,y=Number))+
  geom_col(color="white",aes(fill=Cell_type),position="fill")+
  scale_fill_igv()+
   # labs(x="",y="",title="Cluster Distribution in CMS")+
  theme_classic()+
  theme(title = element_text(size = 20,
                             color = "black"),
        
        plot.title = element_text(hjust = 0.5),
        
        axis.text.x = element_text(size = 10 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 90), #角度
        axis.text.y = element_text(size = 20,  
                                   color = "black",
                                   vjust = 0.4))#+
theme_void()
cmstypes1=read.csv("/home/panxl/CRC2/human_cmstypes.csv")
cmstypes1$barcode=gsub("-",".",cmstypes1$barcode)
########human CMS1
deconvolute_CMS1=deconvolute_human1[deconvolute_human1$Sample%in%cmstypes1[cmstypes1$type=="CMS1",]$barcode,]
ggplot(data=deconvolute_CMS1,aes(x=Sample,y=Number))+
  geom_col(color="white",aes(fill=Cell_type),position="fill")+
  scale_fill_igv()+
  # scale_fill_manual(
  #   labels=c("Cluster 1","Cluster 2",
  #            "Cluster 3","Cluster 4",
  #            "Cluster 5"), #图例标签
  #   values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
   labs(x="",y="",title="CMS1")+
  theme_classic()+
  theme(title = element_text(size = 20,
                             color = "black"),
        
        plot.title = element_text(hjust = 0.5),
        
        axis.text.x = element_text(size = 15 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 90), #角度
        axis.text.y = element_text(size = 20,  
                                   color = "black",
                                   vjust = 0.4))#+
theme_void()

########human CMS2
deconvolute_CMS2=deconvolute_human1[deconvolute_human1$Sample%in%cmstypes1[cmstypes1$type=="CMS2",]$barcode,]
ggplot(data=deconvolute_CMS2,aes(x=Sample,y=Number))+
  geom_col(color="white",aes(fill=Cell_type),position="fill")+
  scale_fill_igv()+
  # scale_fill_manual(
  #   labels=c("Cluster 1","Cluster 2",
  #            "Cluster 3","Cluster 4",
  #            "Cluster 5"), #图例标签
  #   values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
  labs(x="",y="",title="CMS2")+
  theme_classic()+
  theme(title = element_text(size = 20,
                             color = "black"),
        
        plot.title = element_text(hjust = 0.5),
        
        axis.text.x = element_text(size = 8 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 90), #角度
        axis.text.y = element_text(size = 20,  
                                   color = "black",
                                   vjust = 0.4))#+
theme_void()

########human CMS3
deconvolute_CMS3=deconvolute_human1[deconvolute_human1$Sample%in%cmstypes1[cmstypes1$type=="CMS3",]$barcode,]
ggplot(data=deconvolute_CMS3,aes(x=Sample,y=Number))+
  geom_col(color="white",aes(fill=Cell_type),position="fill")+
  scale_fill_igv()+
  # scale_fill_manual(
  #   labels=c("Cluster 1","Cluster 2",
  #            "Cluster 3","Cluster 4",
  #            "Cluster 5"), #图例标签
  #   values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
 labs(x="",y="",title="CMS3")+
  theme_classic()+
  theme(title = element_text(size = 20,
                             color = "black"),
        
        plot.title = element_text(hjust = 0.5),
        
        axis.text.x = element_text(size = 15 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 90), #角度
        axis.text.y = element_text(size = 20,  
                                   color = "black",
                                   vjust = 0.4))#+
theme_void()
########human CMS4
deconvolute_CMS4=deconvolute_human1[deconvolute_human1$Sample%in%cmstypes1[cmstypes1$type=="CMS4",]$barcode,]
ggplot(data=deconvolute_CMS4,aes(x=Sample,y=Number))+
  geom_col(color="white",aes(fill=Cell_type),position="fill")+
  scale_fill_igv()+
  # scale_fill_manual(
  #   labels=c("Cluster 1","Cluster 2",
  #            "Cluster 3","Cluster 4",
  #            "Cluster 5"), #图例标签
  #   values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
 labs(x="",y="",title="CMS4")+
  theme_classic()+
  theme(title = element_text(size = 20,
                             color = "black"),
        
        plot.title = element_text(hjust = 0.5),
        
        axis.text.x = element_text(size = 10 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 90), #角度
        axis.text.y = element_text(size = 20,  
                                   color = "black",
                                   vjust = 0.4))#+
theme_void()
########################### 19.4 反卷积 mean ###############
###mouse
deconvolute_mouse=read.csv("/home/panxl/CRC2/deconvolute/mouse_CIBERSORT-Results.txt",sep="\t")
rownames(deconvolute_mouse)=deconvolute_mouse[,1]
deconvolute_mouse=deconvolute_mouse[,-c(1,24,25,26)]
deconvolute_mouse1=as.data.frame(t(deconvolute_mouse))
deconvolute_mouse1$Mouse=apply(deconvolute_mouse1, 1, function(x){mean(as.numeric(x))})
###human
deconvolute_human=read.csv("/home/panxl/CRC2/deconvolute/human_CIBERSORT-Results.txt",sep="\t")
rownames(deconvolute_human)=deconvolute_human[,1]
deconvolute_human=deconvolute_human[,-c(1,24,25,26)]
deconvolute_human1=as.data.frame(t(deconvolute_human))
deconvolute_human1$Human=apply(deconvolute_human1,1,function(x){mean(as.numeric(x))})
#######CMS
cmstypes1=read.csv("/home/panxl/CRC2/human_cmstypes.csv")
cmstypes1$barcode=gsub("-",".",cmstypes1$barcode)
deconvolute_CMS1=deconvolute_human1[,colnames(deconvolute_human1)%in%cmstypes1[cmstypes1$type=="CMS1",]$barcode]
deconvolute_CMS1$CMS1=rowMeans(deconvolute_CMS1)
deconvolute_CMS2=deconvolute_human1[,colnames(deconvolute_human1)%in%cmstypes1[cmstypes1$type=="CMS2",]$barcode]
deconvolute_CMS2$CMS2=rowMeans(deconvolute_CMS2)
deconvolute_CMS3=deconvolute_human1[,colnames(deconvolute_human1)%in%cmstypes1[cmstypes1$type=="CMS3",]$barcode]
deconvolute_CMS3$CMS3=rowMeans(deconvolute_CMS3)
deconvolute_CMS4=deconvolute_human1[,colnames(deconvolute_human1)%in%cmstypes1[cmstypes1$type=="CMS4",]$barcode]
deconvolute_CMS4$CMS4=rowMeans(deconvolute_CMS4)
data=data.frame(cell_type=rownames(deconvolute_mouse1),Mouse=deconvolute_mouse1$Mouse,Human=deconvolute_human1$Human,
                CMS1=deconvolute_CMS1$CMS1,CMS2=deconvolute_CMS2$CMS2,CMS3=deconvolute_CMS3$CMS3,CMS4=deconvolute_CMS4$CMS4)
#data=data[-11,]
list=c(2:7)
DATA=lapply(list, function(x){
  ll=data[,c(1,x)]
  ll$type=colnames(ll)[2]
  ff=data.frame(Cell_type=ll[,1],Number=ll[,2],Sample=ll[,3])
  return(ff)
})%>%dplyr::bind_rows()
ggplot(data=DATA,aes(x=Sample,y=Number))+
  geom_col(color="white",aes(fill=Cell_type),position="fill")+
  scale_fill_igv()+
  # scale_fill_manual(
  #   labels=c("Cluster 1","Cluster 2",
  #            "Cluster 3","Cluster 4",
  #            "Cluster 5"), #图例标签
  #   values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
  # labs(x="",y="",title="Cluster Distribution in CMS")+
  theme_classic()+
  theme(title = element_text(size = 20,
                             color = "black"),
        
        plot.title = element_text(hjust = 0.5),
        
        axis.text.x = element_text(size = 15 ,# 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 45), #角度
        axis.text.y = element_text(size = 20,  
                                   color = "black",
                                   vjust = 0.4))#+
theme_void()
#############################20.NACC pathway #####################################################
###########################20.1NACC--Logfc 相关性散点图###############################
library(readxl)
hom=read_excel("~/CRC2/zz.xlsx")
hom1=hom[!(hom$NACC=='NA'),]
list=hom1$`Human ensembl id`
results=lapply(list, function(x){
  hom2=hom1[hom1$`Human ensembl id`==x,]
  hom2$'Human fold change'=as.numeric(hom2$'Human fold change')
  hom2$'Mouse fold change'=as.numeric(hom2$'Mouse fold change')
  hom2$hom=abs(hom2$'Human fold change'-hom2$'Mouse fold change')
  return(hom2)
})%>%dplyr::bind_rows()
results$NACC=as.numeric(results$NACC)
#save(results,file = "/home/panxl/CRC2/NACC_HOM.RData")
cor(results$NACC, results$hom)#0.3529
cor.test(results$NACC, results$hom)#0.3529 p-value <2e-16
ggplot(results, aes(x=NACC,y=log10(hom))) +
  geom_point(size=1) +
  # scale_color_uchicago() +
  annotate("text", x = 6, y = -3.5, label = "r=0.3529, p-value <2e-16",
           color="#350E20FF",size = 7 )+# x = -3, y = -0.3
 # geom_smooth(method=lm,color="red")+
 labs(x="NACC",y="absolute deviation of fold change",title="")+
  theme_classic(base_line_size = 1) +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24, color="black"),
        title = element_text(size = 30,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"))
################################20.2 NACC  Kegg Pathway############################# 
load("/home/panxl/CRC2/NACC_HOM.RData")
NACC_HOM=results
Tp=mean(NACC_HOM$NACC)#2.242
sdp=sd(NACC_HOM$NACC)#0.7933
##########.1kegg pathway
hom_gene=read.csv("/home/panxl/CRC2/all_cluater")
library(org.Hs.eg.db)
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                          keys = hom_gene$Gene.stable.ID,
                          keytype = "ENSEMBL",
                          column = "ENTREZID")#
DEG.entrez_id = na.omit(DEG.entrez_id)#去除NA 13646
hom_kegg <- enrichKEGG(gene       = DEG.entrez_id,
                      keyType       = "kegg",
                      organism      = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)
kegg_result=hom_kegg@result
#save(hom_kegg,file = "/home/panxl/CRC2/hom_kegg.RData")
load("/home/panxl/CRC2/hom_kegg.RData")
kegg_result=hom_kegg@result
kegg_pathway=kegg_result$ID
z_score_kegg=lapply(kegg_pathway, function(x){
  l=kegg_result[kegg_result$ID=="hsa04060",]$geneID
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=kegg_result[kegg_result$ID=="hsa04060",]$Description
  nacc=NACC_HOM[NACC_HOM$`Human ensembl id`%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(NACC_HOM$NACC))/(sd(NACC_HOM$NACC))
  kk=data.frame(pathway_id=x,z_score=z,gene_number=num,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()  
write.csv(z_score_kegg,file = "/home/panxl/CRC2/z_score_kegg.csv")
#save(z_score_kegg,file = "/home/panxl/CRC2/z_score_kegg.RData")
colnames(z_score_kegg)
##########.2 Kegg NACC Pathway 画图
load("/home/panxl/CRC2/z_score_kegg.RData")
Data_pathway=z_score_kegg[z_score_kegg$pathway_id%in%c("hsa00591",
                                                       "hsa00564",
                                                       "hsa00511",
                                                        "hsa00561",
                                                        "hsa00051",
                                                       "hsa03030",
                                                       "hsa03430",
                                                       "hsa03010",
                                                       "hsa05330",
                                                       "hsa05332",
                                                       "hsa04060",
                                                        "hsa04062","hsa05320"

                                                       
),]

library(ggrepel)  #标签用
ggplot(z_score_kegg, aes(x=gene_number,y=z_score,color=z_score)) +
  geom_point(size=2.5) +
  scale_color_gradient2(low ="#8AC1ED",mid="#67B738",high="#D3673C" )+
  ggrepel::geom_label_repel(
    aes(label =Data_pathway$Description),
    data =Data_pathway,
    color="black",
    size = 4,
    box.padding=unit(0.5, "lines"), point.padding=unit(0.5, "lines"), 
    segment.color = "#A9A9A9", segment.size = 1,
     #arrow = arrow(length=unit(0.01, "npc")),force = 1, max.iter = 3e3, 
  )+
 labs(x="Number of genes in KEGG category",y="Z-score",title="KEGG pathway")+
  theme_classic(base_line_size = 1) +
  theme(#axis.title = element_text(size=20),
        #axis.text = element_text(size=14, color="black"),
        title = element_text(size =28,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 28,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 28,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"))#172 x 144
################################20.3 NACC  GO pathway######################################
load("/home/panxl/CRC2/NACC_HOM.RData")
NACC_HOM=results
Tp=mean(NACC_HOM$NACC)#2.242
sdp=sd(NACC_HOM$NACC)#0.7933
###############GO pathway
hom_gene=read.csv("/home/panxl/CRC2/all_cluater")
library(org.Hs.eg.db)
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = hom_gene$Gene.stable.ID,
                       keytype = "ENSEMBL",
                       column = "ENTREZID")#
DEG.entrez_id = na.omit(DEG.entrez_id)#去除NA 13646
hom_GO<- enrichGO(gene = DEG.entrez_id,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.01)
hom_GO=hom_kegg
GO_result=hom_GO@result
#save(hom_GO,file = "/home/panxl/CRC2/hom_GO.RData")
load("/home/panxl/CRC2/hom_GO.RData")
GO_result=hom_GO@result#2731
GO_pathway=GO_result$ID
z_score_GO=lapply(GO_pathway, function(x){
  l=GO_result[GO_result$ID==x,]$geneID
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=GO_result[GO_result$ID==x,]$Description
  nacc=NACC_HOM[NACC_HOM$`Human ensembl id`%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(NACC_HOM$NACC))/(sd(NACC_HOM$NACC))
  kk=data.frame(pathway_id=x,z_score=z,gene_number=num,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()  
#write.csv(z_score_GO,file = "/home/panxl/CRC2/z_score_GO.csv")
#save(z_score_GO,file = "/home/panxl/CRC2/z_score_GO.RData")
colnames(z_score_GO)
############GO画图
load("/home/panxl/CRC2/z_score_GO.RData")
Data_gopathway=z_score_GO[abs(z_score_GO$z_score)>1,]
ggplot(z_score_GO, aes(x=gene_number,y=z_score,color=z_score)) +
  geom_point(size=1) +
  scale_color_gradient2(low ="#8AC1ED",mid="#67B738",high="#D3673C" )+
  ggrepel::geom_label_repel(
    aes(label =Data_gopathway$Description),
    data =Data_gopathway,
    color="black",
    size = 3,
    box.padding=unit(0.5, "lines"), point.padding=unit(0.5, "lines"), 
    segment.color = "#A9A9A9", segment.size = 1,
    #arrow = arrow(length=unit(0.01, "npc")),force = 1, max.iter = 3e3, 
  )+
  # annotate("text", x = -0.2, y = -3.5, label = "r=0.3529, p-value <2e-16",
  #          color="#350E20FF",size = 7 )+# x = -3, y = -0.3
  #geom_smooth(method=lm,color="red")+
  labs(x="Number of genes in GO category",y="Z-score",title="GO pathway")+
  theme_classic(base_line_size = 1) +
  theme(#axis.title = element_text(size=20),
    #axis.text = element_text(size=14, color="black"),
    title = element_text(size =28,
                         color = "black",
                         hjust = 0.5),
    axis.text.x = element_text(size = 28,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
    axis.text.y = element_text(size = 28,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
    legend.title = element_blank(),
    legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
    axis.ticks.length.x = unit(2,"mm"))
   
############################# 21 NACC CMS 热图######################################################
############################# 21.1 计算CMS gene NACC ####################################
##生科院服务器
# cmstypes1=read.csv("/home/panxl/CRC2/human_cmstypes.csv")
# human_data <- read.csv("/home/shimw/project/conservation/filter_human_data.csv",row.names = 1, check.names = F)
# human_ortho=read.csv("/home/panxl/CRC2/all_cluater")
# mouse_data <- read.csv("/home/shimw/project/conservation/filter_mouse_data.csv",row.names = 1, check.names = F)
# mouse_ortho = mouse_data[rownames(mouse_data) %in% human_ortho$Mouse.gene.stable.ID,]
# mouse_pairs_cor=cor(t(mouse_ortho))
# save(mouse_pairs_cor,file = "/home/panxl/CRC2/mouse_pairs_cor.RData")
# load("/home/panxl/CRC2/mouse_pairs_cor.RData")
# load("/home/xlpan/CRC/NACC_CMS/mouse_pairs_cor.RData")
# cmstypes1=read.csv("/home/xlpan/CRC/NACC_CMS/human_cmstypes.csv")
# human_data <- read.csv("/home/xlpan/CRC/NACC_CMS/filter_human_data.csv",row.names = 1, check.names = F)
# human_ortho=read.csv("/home/xlpan/CRC/NACC_CMS/all_cluater")
# #####CMS1
# CMS1_data=human_data[,colnames(human_data)%in%cmstypes1[cmstypes1$type=="CMS1",]$barcode]
# CMS1_data=CMS1_data[rownames(CMS1_data) %in% human_ortho$Gene.stable.ID,]
# CMS1_data_t=t(CMS1_data)
# CMS1_pairs_cor = cor(CMS1_data_t)
# library(dplyr)
# CMS1_pairs_top20_cor <- purrr::map(colnames(CMS1_pairs_cor),function(m){
#   x = CMS1_pairs_cor[,m]
#   x= sort(x, decreasing = TRUE)
#   names(x[2:21])
#   tmp_orth <- human_ortho[match(c(m,names(x[2:21])),human_ortho$Gene.stable.ID),]$Mouse.gene.stable.ID
#   kk = dist(rbind(x[2:21],mouse_pairs_cor[tmp_orth[2:21],tmp_orth[1]]))
#   tibble("gene"=m,
#          "d"=kk[1])
# })%>%bind_rows()
# CMS1_pairs_top20_cor = CMS1_pairs_top20_cor[match(human_ortho$Gene.stable.ID,CMS1_pairs_top20_cor$gene),]
# readr::write_csv(CMS1_pairs_top20_cor,"/home/xlpan/CRC/NACC_CMS/CMS1_pairs_distance.csv")
# CMS1_mouse_pairs_top20_cor <- purrr::map(colnames(mouse_pairs_cor),function(m){
#   x = mouse_pairs_cor[,m]
#   x= sort(x, decreasing = TRUE)
#   names(x[2:21])
#   tmp_orth <- human_ortho[match(c(m,names(x[2:21])),human_ortho$Mouse.gene.stable.ID),]$Gene.stable.ID
#   kk = dist(rbind(x[2:21],CMS1_pairs_cor[tmp_orth[2:21],tmp_orth[1]]))
#   tibble("gene"=m,
#          "d"=kk[1])
# })%>%bind_rows()
# CMS1_mouse_pairs_top20_cor = CMS1_mouse_pairs_top20_cor[match(human_ortho$Mouse.gene.stable.ID,CMS1_mouse_pairs_top20_cor$gene),]
# readr::write_csv(CMS1_mouse_pairs_top20_cor,"/home/xlpan/CRC/NACC_CMS/CMS1_mouse_pairs_distance.csv")
# CMS1_pairs_top20_cor = read.csv("/home/xlpan/CRC/NACC_CMS/CMS1_pairs_distance.csv")
# CMS1_mouse_pairs_top20_cor = read.csv("/home/xlpan/CRC/NACC_CMS/CMS1_mouse_pairs_distance.csv")
# CMS1_pairs_top20_cor = CMS1_pairs_top20_cor[match(human_ortho$Gene.stable.ID,CMS1_pairs_top20_cor$gene),]
# CMS1_mouse_pairs_top20_cor = CMS1_mouse_pairs_top20_cor[match(human_ortho$Mouse.gene.stable.ID,CMS1_mouse_pairs_top20_cor$gene),]
# CMS1_pairs_top20_cor$dd = (CMS1_pairs_top20_cor$d + CMS1_mouse_pairs_top20_cor$d)/2
# CMS1_NACC_score = CMS1_pairs_top20_cor[, c("gene", "dd")]
# colnames(CMS1_NACC_score) = c("gene", "NACC")
# readr::write_csv(CMS1_NACC_score,"/home/xlpan/CRC/NACC_CMS/CMS1_mouse_crc_nacc.csv")
###################################Supplementary table S5.CMS KEGG Z-score#########################################
#######加载数据
CMS1=read.csv("/home/panxl/CRC2/NACC/CMS1_mouse_crc_nacc.csv")
CMS2=read.csv("/home/panxl/CRC2/NACC/CMS2_mouse_crc_nacc.csv")
CMS3=read.csv("/home/panxl/CRC2/NACC/CMS3_mouse_crc_nacc.csv")
CMS4=read.csv("/home/panxl/CRC2/NACC/CMS4_mouse_crc_nacc.csv")
load("/home/panxl/CRC2/NACC_HOM.RData")
load("/home/panxl/CRC2/mouse_gse_kegg.Rdata")
load("/home/panxl/CRC2/human_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS1_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS2_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS3_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS4_gse_kegg.Rdata")
CMS1_gse=CMS1_gse_kegg@result
df=CMS1_gse$ID
CMS1_gse=CMS1_gse_kegg@result
CMS1_z_score=lapply(df, function(x){
  CMS1_gse=CMS1_gse_kegg@result
  l=CMS1_gse[CMS1_gse$ID==x,]$core_enrichment
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=CMS1_gse[CMS1_gse$ID==x,]$Description
  nacc=CMS1[CMS1$gene%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(CMS1$NACC))/(sd(CMS1$NACC))
  kk=data.frame(CMS1=z,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()
CMS1=left_join(CMS1_gse,CMS1_z_score,by="Description")
colnames(CMS1)[12]="NACC"
CMS1$Cluster="CMS1"
CMS2_z_score=lapply(df, function(x){
  CMS2_gse=CMS2_gse_kegg@result
  l=CMS2_gse[CMS2_gse$ID==x,]$core_enrichment
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=CMS2_gse[CMS2_gse$ID==x,]$Description
  nacc=CMS2[CMS2$gene%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(CMS2$NACC))/(sd(CMS2$NACC))
  kk=data.frame(CMS2=z,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()

CMS2=left_join(CMS2_gse,CMS2_z_score,by="Description")
colnames(CMS2)[12]="NACC"
CMS2$Cluster="CMS2"


CMS3_z_score=lapply(df, function(x){
  CMS3_gse=CMS3_gse_kegg@result
  l=CMS3_gse[CMS3_gse$ID==x,]$core_enrichment
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=CMS3_gse[CMS3_gse$ID==x,]$Description
  nacc=CMS3[CMS3$gene%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(CMS3$NACC))/(sd(CMS3$NACC))
  kk=data.frame(CMS3=z,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()

CMS3=left_join(CMS3_gse,CMS3_z_score,by="Description")
colnames(CMS3)[12]="NACC"
CMS3$Cluster="CMS3"


CMS4_z_score=lapply(df, function(x){
  CMS4_gse=CMS4_gse_kegg@result
  l=CMS4_gse[CMS4_gse$ID==x,]$core_enrichment
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=CMS4_gse[CMS4_gse$ID==x,]$Description
  nacc=CMS4[CMS4$gene%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(CMS1$NACC))/(sd(CMS1$NACC))
  kk=data.frame(CMS4=z,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()


CMS4=left_join(CMS4_gse,CMS4_z_score,by="Description")
colnames(CMS4)[12]="NACC"
CMS4$Cluster="CMS4"



human_z_score=lapply(df, function(x){
  human_gse=human_gse_kegg@result
  l=human_gse[human_gse$ID==x,]$core_enrichment
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=human_gse[human_gse$ID==x,]$Description
  nacc=results[results$`Human ensembl id`%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(CMS1$NACC))/(sd(CMS1$NACC))
  kk=data.frame(Human=z,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()
human=left_join(human_gse,human_z_score,by="Description")
colnames(human)[12]="NACC"
human$Cluster="Human"
data1=rbind(CMS1,CMS2)
data2=rbind(data1,CMS3)
data3=rbind(data2,CMS4)
data=rbind(data3,human)
#write.csv(data,file = "/home/panxl/CRC2/CMS_KEGG_Z-score.csv")
############################# 21.2  CMS kegg pathway z_score ####################################
#######加载数据
CMS1=read.csv("/home/panxl/CRC2/NACC/CMS1_mouse_crc_nacc.csv")
CMS2=read.csv("/home/panxl/CRC2/NACC/CMS2_mouse_crc_nacc.csv")
CMS3=read.csv("/home/panxl/CRC2/NACC/CMS3_mouse_crc_nacc.csv")
CMS4=read.csv("/home/panxl/CRC2/NACC/CMS4_mouse_crc_nacc.csv")
load("/home/panxl/CRC2/NACC_HOM.RData")
unique_cluster_pathway=read.csv("/home/panxl/CRC2/unique_cluster4_pathway.kegg.csv")
load("/home/panxl/CRC2/mouse_gse_kegg.Rdata")
load("/home/panxl/CRC2/human_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS1_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS2_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS3_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS4_gse_kegg.Rdata")
df=unique_cluster_pathway$ID

CMS1_z_score=lapply(df, function(x){
  CMS1_gse=CMS1_gse_kegg@result
  l=CMS1_gse[CMS1_gse$ID==x,]$core_enrichment
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=CMS1_gse[CMS1_gse$ID==x,]$Description
  nacc=CMS1[CMS1$gene%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(CMS1$NACC))/(sd(CMS1$NACC))
  kk=data.frame(CMS1=z,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()

CMS2_z_score=lapply(df, function(x){
  CMS2_gse=CMS2_gse_kegg@result
  l=CMS2_gse[CMS2_gse$ID==x,]$core_enrichment
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=CMS2_gse[CMS2_gse$ID==x,]$Description
  nacc=CMS2[CMS2$gene%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(CMS2$NACC))/(sd(CMS2$NACC))
  kk=data.frame(CMS2=z,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()



CMS3_z_score=lapply(df, function(x){
  CMS3_gse=CMS3_gse_kegg@result
  l=CMS3_gse[CMS3_gse$ID==x,]$core_enrichment
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=CMS3_gse[CMS3_gse$ID==x,]$Description
  nacc=CMS3[CMS3$gene%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(CMS3$NACC))/(sd(CMS3$NACC))
  kk=data.frame(CMS3=z,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()



CMS4_z_score=lapply(df, function(x){
  CMS4_gse=CMS4_gse_kegg@result
  l=CMS4_gse[CMS4_gse$ID==x,]$core_enrichment
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=CMS4_gse[CMS4_gse$ID==x,]$Description
  nacc=CMS4[CMS4$gene%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(CMS1$NACC))/(sd(CMS1$NACC))
  kk=data.frame(CMS4=z,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()



human_z_score=lapply(df, function(x){
  human_gse=human_gse_kegg@result
  l=human_gse[human_gse$ID==x,]$core_enrichment
  ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
  ENSEMBL_id =bitr(ll$ENTREZID,
                   fromType="ENTREZID", 
                   toType= "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
  Description=human_gse[human_gse$ID==x,]$Description
  nacc=results[results$`Human ensembl id`%in%ENSEMBL_id$ENSEMBL,]$NACC
  num=length(nacc)
  mean_nacc=mean(nacc)
  z=(mean_nacc-mean(CMS1$NACC))/(sd(CMS1$NACC))
  kk=data.frame(Human=z,Description=Description)
  return(kk)
})%>%dplyr::bind_rows()

save(CMS1_z_score,file = "/home/panxl/CRC2/CMS1_z_score.RData")
save(CMS2_z_score,file = "/home/panxl/CRC2/CMS2_z_score.RData")
save(CMS3_z_score,file = "/home/panxl/CRC2/CMS3_z_score.RData")
save(CMS4_z_score,file = "/home/panxl/CRC2/CMS4_z_score.RData")
save(human_z_score,file = "/home/panxl/CRC2/human_z_score.RData")
####
load("/home/panxl/CRC2/CMS1_z_score.RData")
load("/home/panxl/CRC2/CMS2_z_score.RData")
load("/home/panxl/CRC2/CMS3_z_score.RData")
load("/home/panxl/CRC2/CMS4_z_score.RData")
load("/home/panxl/CRC2/human_z_score.RData")

a=left_join(CMS1_z_score,CMS2_z_score,by = "Description")
b=left_join(CMS3_z_score,a,by = "Description")
c=left_join(CMS4_z_score,b,by = "Description")
d=left_join(human_z_score,c,by = "Description")
unique_cluster_pathway=read.csv("/home/panxl/CRC2/unique_cluster4_pathway.kegg.csv")
NACC=left_join(d,unique_cluster_pathway,by = "Description")
NACC=NACC[,c(1,5,6,4,3,2,10)]
colnames(NACC)
#write.csv(NACC, "/home/panxl/CRC2/NACC_cms.csv")
###################21.3  CMS NACC z_score 热图###################################
library(tidyverse)
library(pheatmap)
NACC=read.csv("/home/panxl/CRC2/NACC_cms.csv",row.names = 1)
table(NACC$type)
rownames(NACC)=NACC$Description
data=NACC[,c(1,2,3,4,5)]
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

ha1 = HeatmapAnnotation(Group=factor(c( "Human",  "CMS1" , "CMS2" , "CMS3" , "CMS4" ),
                                     levels = c( "Human",  "CMS1" , "CMS2" , "CMS3" , "CMS4" )),
                        col = list(Group = c("Human"="#78a355",'CMS1' = '#b97573',
                                             'CMS2' = '#686954', 'CMS3' = '#5b5e66', 'CMS4'= '#816147')),show_legend = F)


ha2 = rowAnnotation(Cluster=factor(c(rep("Cluster1",25),rep("Cluster2",24),rep("Cluster3",21),rep("Cluster4",3)),
                                   levels = c("Cluster1","Cluster2","Cluster3","Cluster4")),
                      annotation_name_rot = 45,
                      col = list(Cluster=c("Cluster1"="#CB9C7A","Cluster2"="#8696a7",
                                         "Cluster3"="#CDB97D","Cluster4"="#7b8b6f")),show_legend = F)

#df_sub1=t(scale(t(df_sub)))"#A59B95"
#library()
  ht_list=Heatmap(data, cluster_rows = F,
        cluster_columns = F,
        column_names_rot = 45,
        show_column_names = T,
        show_row_names = T,
        heatmap_legend_param = list(
          title = "NACC", at = c(-3, 0,1),
          labels = c("-3","", "1")
        ),
    
      col=c("#87CEFA", "white", "#CC2121"),
        top_annotation = ha1,
        left_annotation = ha2,
        row_names_gp = gpar(fontsize = 9), 
        show_heatmap_legend = T)

draw(ht_list, heatmap_legend_side = "left")
############################### 21.4 验证NACC热图---Wnt pathway #############################
HOM=read.csv("/home/panxl/CRC2/all_cluater.csv")
load("/home/panxl/CRC2/human_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS1_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS2_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS3_gse_kegg.Rdata")
load( "/home/panxl/CRC2/CMS4_gse_kegg.Rdata")
CMS1_gse=CMS1_gse_kegg@result
l=CMS1_gse[CMS1_gse$ID=="hsa04310",]$core_enrichment
ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))

Wnt_gene_ENSEMBL= mapIds(x = org.Hs.eg.db,
                         keys = ll$ENTREZID,
                         keytype = "ENTREZID",
                         column ="ENSEMBL")#27
Wnt_gene_ENSEMBL = na.omit(Wnt_gene_ENSEMBL)#去除NA 27


cms1=HOM[HOM$Gene.stable.ID%in%Wnt_gene_ENSEMBL,]
cms_1=CMS1[CMS1$gene%in%Wnt_gene_ENSEMBL,]
colnames(cms_1)=c("gene","CMS1")
cms_2=CMS2[CMS2$gene%in%Wnt_gene_ENSEMBL,]
colnames(cms_2)=c("gene","CMS2")
cms_3=CMS3[CMS3$gene%in%Wnt_gene_ENSEMBL,]
colnames(cms_3)=c("gene","CMS3")
cms_4=CMS4[CMS4$gene%in%Wnt_gene_ENSEMBL,]
colnames(cms_4)=c("gene","CMS4")

load("/home/panxl/CRC2/NACC_HOM.RData")
hh=results[results$`Human ensembl id`%in%Wnt_gene_ENSEMBL,][,c(3,4,9)]
table(cms1)






cmstypes1=read.csv("/home/panxl/CRC2/human_cmstypes.csv")
human_data <- read.csv("/home/panxl/CRC2/filter_human_data.csv",row.names = 1, check.names = F)
human_ortho=read.csv("/home/panxl/CRC2/all_cluater")
CMS2_data=human_data[,colnames(human_data)%in%cmstypes1[cmstypes1$type=="CMS2",]$barcode]
CMS2_data=CMS2_data[rownames(CMS2_data) %in% human_ortho$Gene.stable.ID,]
############################### 21.5 CMS NACC 密度图 #############################

CMS1=read.csv("/home/panxl/CRC2/NACC/CMS1_mouse_crc_nacc.csv")
CMS2=read.csv("/home/panxl/CRC2/NACC/CMS2_mouse_crc_nacc.csv")
CMS3=read.csv("/home/panxl/CRC2/NACC/CMS3_mouse_crc_nacc.csv")
CMS4=read.csv("/home/panxl/CRC2/NACC/CMS4_mouse_crc_nacc.csv")
load("/home/panxl/CRC2/NACC_HOM.RData")
CMS1$Type="CMS1"
CMS2$Type="CMS2"
CMS3$Type="CMS3"
CMS4$Type="CMS4"
data1=rbind(CMS1,CMS2)
data2=rbind(data1,CMS3)
data=rbind(data2,CMS4)
summary(data[data$Type=="CMS1",])
summary(data[data$Type=="CMS2",])
summary(data[data$Type=="CMS3",])
summary(data[data$Type=="CMS4",])
t.test(data[data$Type=="CMS1",]$NACC,data[data$Type=="CMS4",]$NACC)
##全部基因
ggplot(data=data,aes(x=NACC))+
  geom_density(alpha=0.3,#透明度
               aes(fill=Type))+
  scale_fill_manual(
    labels=c("CMS1","CMS2","CMS3","CMS4"), #图例标签
    values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f"))+#颇色
  
  theme_classic()+ 
  theme(title = element_text(size = 22,
                             color = "black"),
        axis.title = element_text(size=28),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )
  theme(legend.position='none')

############################ 22. Chemokine signaling pathway in mouse ############################
load("/home/panxl/CRC2/hom_kegg.RData")
hom=read.csv("/home/panxl/CRC2/all_cluater")
hom_kegg=hom_kegg@result
l=hom_kegg[hom_kegg$ID=="hsa04062",]$geneID
ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
ENSEMBL_id =bitr(ll$ENTREZID,
                 fromType="ENTREZID", 
                 toType= "ENSEMBL", 
                 OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_id=hom[hom$Gene.stable.ID%in%ENSEMBL_id$ENSEMBL,]
data=data.frame(table(gene_id$cluster)) 
colnames(data)=c("Cluster",'Number')
data$Cluster=gsub("Cluster 1: Up-regulated in both species","Cluster 1",data$Cluster)
data$Cluster=gsub("Cluster 2: Up-regulated in mouse reverse in human","Cluster 2",data$Cluster)
data$Cluster=gsub("Cluster 3: Down-regulated in both species","Cluster 3",data$Cluster)
data$Cluster=gsub("Cluster 4: Up-regulated in human reverse in mouse","Cluster 4",data$Cluster)
data$Cluster=gsub("Cluster 5: Not-significant in both species","Cluster 5",data$Cluster)
data$Cluster=as.factor(data$Cluster)
ggplot(data=data,aes(x=Cluster,y=Number, fill=Cluster))+
  geom_bar(stat="identity")+
         # geom_col(color="white",aes(fill=Number),position="fill")+
         scale_fill_manual(
           labels=c("Cluster 1","Cluster 2",
                    "Cluster 3","Cluster 4",
                    "Cluster 5"), #图例标签
           values=c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"))+#颇色
         labs(x="",y="",title="Chemokine signaling pathway in mouse")+
  theme_classic(base_line_size = 1) +
  
  theme(#axis.title = element_text(size=20),
    #axis.text = element_text(size=14, color="black"),
    title = element_text(size = 20,
                         color = "black",
                         hjust = 0.5),
    axis.text.x = element_text(size = 25,colour = "black",  vjust = 0.5, # 位置
                               hjust = 0.5, 
                               angle = 45), #角度)),
    axis.text.y = element_text(size = 25,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
    legend.title = element_blank(),
    legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
    axis.ticks.length.x = unit(2,"mm"))##967 954

#############################23.1 Cluster2 CMS fpkm box 图#################
hom=read.csv("/home/panxl/CRC2/all_cluater.csv")
Cluster2_gene=hom[hom$cluster=="Cluster 2",]$Gene.stable.ID
df=c("CMS1","CMS2","CMS3","CMS4")
box_Cluster2_result=lapply(df, function(x){
  CMS_FPKM=read.csv(paste("/home/panxl/CRC2/fpkm_",x,".csv",sep = ""),row.names = 1)
  CMS_FPKM1=CMS_FPKM[rownames(CMS_FPKM)%in%Cluster2_gene,]
  CMS_FPKM1$mean=rowMeans(CMS_FPKM1)
  box_data=data.frame(gene_ID=rownames(CMS_FPKM1),mean=CMS_FPKM1$mean,Type=x)
  return(box_data)
})%>%dplyr::bind_rows()
#save(box_Cluster2_result,file = "/home/panxl/CRC2/box_Cluster2_result.RData")
############画图
load("/home/panxl/CRC2/box_Cluster2_result.RData")
library(ggsci)
library(ggsignif)
library(ggpubr)
boxplot(log2(mean+1) ~ Type,data=box_Cluster2_result,
        col= c('#b97573','#686954',"#CDB97D",'#5b5e66','#816147'),
        # las = 2,
       ylab="log2(FPKM)",xlab = NULL,outline = F,xaxt = "n",cex.lab=1.3)

ggplot(box_Cluster2_result,aes(x=Type,y=log2(mean+1)))+
geom_boxplot(aes(fill=Type))+
  labs(y="log2(FPKM)",x="",title="Cluster2")+
geom_signif(comparisons = list(c("CMS4","CMS1")),
            #y_position = c(0.6,3),#显著性比较
            map_signif_level = T)
#('CMS1' = '#b97573','CMS2' = '#686954', 'CMS3' = '#5b5e66', 'CMS4'= '#816147')
#############################23.2 Cluster3 CMS fpkm box 图##############
hom=read.csv("/home/panxl/CRC2/all_cluater.csv")
Cluster3_gene=hom[hom$cluster=="Cluster 3",]$Gene.stable.ID
df=c("CMS1","CMS2","CMS3","CMS4")
box_Cluster3_result=lapply(df, function(x){
  CMS_FPKM=read.csv(paste("/home/panxl/CRC2/fpkm_",x,".csv",sep = ""),row.names = 1)
  CMS_FPKM1=CMS_FPKM[rownames(CMS_FPKM)%in%Cluster3_gene,]
  CMS_FPKM1$mean=rowMeans(CMS_FPKM1)
  box_data=data.frame(gene_ID=rownames(CMS_FPKM1),mean=CMS_FPKM1$mean,Type=x)
  return(box_data)
})%>%dplyr::bind_rows()
#save(box_Cluster3_result,file = "/home/panxl/CRC2/box_Cluster3_result.RData")
###############画图
load("/home/panxl/CRC2/box_Cluster3_result.RData")
library(ggsci)
library(ggsignif)

boxplot(log2(mean+1) ~ Type,data=box_Cluster3_result,
        col= c('#b97573','#686954',"#CDB97D",'#5b5e66','#816147'),
        # las = 2,
        ylab="log2(FPKM)",xlab = NULL,outline = F,xaxt = "n",cex.lab=1.3)



ggplot(box_Cluster3_result,aes(x=Type,y=log2(mean)))+
  labs(y="log2(FPKM)",x="",title="Cluster3")+
  geom_boxplot(aes(fill=Type))+
  geom_signif(comparisons = list(c("CMS4","CMS1"),c("CMS4","CMS3")),
              #y_position = c(0.6,3),#显著性比较
              map_signif_level = T)




#########################24.wnt_pathway_nes##################################
#################24.1gseKEGG wnt_pathway_nes#################
df=c("CMS1","CMS2","CMS3","CMS4") 
wnt_pathway_nes=lapply(df, function(x){
  resdata=read.csv(paste("/home/panxl/CRC2/human_",x,".csv",sep=""))#亚型重复四次
  gene_res_ENSEMBL = as.character(resdata$gene)
  gene_res_df = bitr(gene_res_ENSEMBL, 
                     fromType="ENSEMBL", 
                     toType="ENTREZID", 
                     OrgDb = "org.Hs.eg.db", drop = TRUE)
  gene_df <- left_join(resdata,gene_res_df,c("gene"="ENSEMBL"))#
  gene_df1=na.omit(gene_df)
  gene_res0.01= as.vector(gene_df1$ENTREZID)#
  geneList<-gene_res0.01
  data(geneList,package = "DOSE")
  geneList = gene_df$log2FoldChange
  names(geneList) = gene_df$ENTREZID
  geneList = sort(geneList,decreasing = T)
  geo_list<- list()
  pvalue_Cutoff=0.05
  CMS3_gse_kegg<-gseKEGG(geneList= geneList,
                         organism = "hsa",
                         keyType  = 'kegg',
                         nPerm  = 1000,
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff = 1,
                         pAdjustMethod     = "BH") ## 0.05有的没有结果 
  CMS3_gse_kegg=CMS3_gse_kegg@result
  ll=CMS3_gse_kegg[CMS3_gse_kegg$ID=="hsa04310",]$NES
  data=data.frame(Type=x,NES=ll)
  return(data)
  
})%>%dplyr::bind_rows()
#save(wnt_pathway_nes,file = "/home/panxl/CRC2/wnt_pathway_nes.RData")
load("/home/panxl/CRC2/wnt_pathway_nes.RData")
#################24.2fgseaRes NES wnt_pathway_nes#################
load("/home/panxl/CRC2/human_examplePathways.RData")

load("/home/panxl/CRC2/unique_cluster_pathway_kegg_annotation.RData")
gene=unique_cluster_pathway_kegg_annotation[unique_cluster_pathway_kegg_annotation$type=="hsa04310",]
gene$geneList
gene_df=bitr(gene$geneList, 
                 fromType="SYMBOL", 
                 toType="ENSEMBL", 
                 OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_df <- left_join(gene,gene_df,c("geneList"="SYMBOL"))#
wnt_geneList =na.omit(gene_df)
wnt_geneList1 =wnt_geneList [!(duplicated(wnt_geneList$geneList)),]


human_DEG=read.csv("/home/panxl/CRC2/human_DEG.csv")
human=human_DEG[human_DEG$logFC>1,]
human_gene=human[human$X%in%wnt_geneList1$ENSEMBL,]
human_gene_df=bitr(human_gene$X, 
                 fromType="ENSEMBL", 
                 toType="ENTREZID", 
                 OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_df <- left_join(human_gene,human_gene_df,c("X"="ENSEMBL"))#
gene_df =na.omit(gene_df)
gene_res= as.vector(gene_df$ENTREZID)#
geneList<-gene_res
data(geneList,package = "DOSE")
geneList = gene_df$logFC
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList,decreasing = T)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = geneList,
                  minSize  = 15,
                  maxSize  = 500)
human_nes=fgseaRes$NES #0.773 




CMS1_DEG=read.csv("/home/panxl/CRC2/human_CMS1.csv")
CMS1=CMS1_DEG[CMS1_DEG$log2FoldChange>1,]
CMS1_gene=CMS1[CMS1$gene%in%wnt_geneList1$ENSEMBL,]
CMS1_gene_df=bitr(CMS1_gene$gene, 
                   fromType="ENSEMBL", 
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_df <- left_join(CMS1_gene,CMS1_gene_df,c("gene"="ENSEMBL"))#
gene_df =na.omit(gene_df)
gene_res= as.vector(gene_df$ENTREZID)#
geneList<-gene_res
data(geneList,package = "DOSE")
geneList = gene_df$log2FoldChange
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList,decreasing = T)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = geneList,
                  minSize  = 15,
                  maxSize  = 500)
CMS1_nes #0.903   hsa04310




CMS2_DEG=read.csv("/home/panxl/CRC2/human_CMS2.csv")
CMS2=CMS2_DEG[CMS2_DEG$log2FoldChange>1,]
CMS2_gene=CMS2[CMS2$gene%in%wnt_geneList1$ENSEMBL,]
CMS2_gene_df=bitr(CMS2_gene$gene, 
                  fromType="ENSEMBL", 
                  toType="ENTREZID", 
                  OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_df <- left_join(CMS2_gene,CMS2_gene_df,c("gene"="ENSEMBL"))#
gene_df =na.omit(gene_df)
gene_res= as.vector(gene_df$ENTREZID)#
geneList<-gene_res
data(geneList,package = "DOSE")
geneList = gene_df$log2FoldChange
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList,decreasing = T)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = geneList,
                  minSize  = 15,
                  maxSize  = 500)
CMS2_nes #0.815   hsa04310


CMS3_DEG=read.csv("/home/panxl/CRC2/human_CMS3.csv")
CMS3=CMS3_DEG[CMS3_DEG$log2FoldChange>1,]
CMS3_gene=CMS3[CMS3$gene%in%wnt_geneList1$ENSEMBL,]
CMS3_gene_df=bitr(CMS3_gene$gene, 
                  fromType="ENSEMBL", 
                  toType="ENTREZID", 
                  OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_df <- left_join(CMS3_gene,CMS3_gene_df,c("gene"="ENSEMBL"))#
gene_df =na.omit(gene_df)
gene_res= as.vector(gene_df$ENTREZID)#
geneList<-gene_res
data(geneList,package = "DOSE")
geneList = gene_df$log2FoldChange
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList,decreasing = T)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = geneList,
                  minSize  = 15,
                  maxSize  = 500)
CMS3_nes #1.144   hsa04310


CMS4_DEG=read.csv("/home/panxl/CRC2/human_CMS4.csv")
CMS4=CMS4_DEG[CMS4_DEG$log2FoldChange>1,]
CMS4_gene=CMS4[CMS4$gene%in%wnt_geneList1$ENSEMBL,]
CMS4_gene_df=bitr(CMS4_gene$gene, 
                  fromType="ENSEMBL", 
                  toType="ENTREZID", 
                  OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_df <- left_join(CMS4_gene,CMS4_gene_df,c("gene"="ENSEMBL"))#
gene_df =na.omit(gene_df)
gene_res= as.vector(gene_df$ENTREZID)#
geneList<-gene_res
data(geneList,package = "DOSE")
geneList = gene_df$log2FoldChange
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList,decreasing = T)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = geneList,
                  minSize  = 15,
                  maxSize  = 500)
CMS4_nes #0.459
############################## 25.1小鼠亚型##############################################
library(CMScaller)
mouse_res=read.csv("/home/panxl/CRC2/mouse78_count.csv",check.names = F,row.names = 1)
mouse_info=read.csv("/home/panxl/CRC2/filter_mouse_info.csv",check.names = F,row.names = 1)
mouse_info_tumor=mouse_info[mouse_info$Group=="Tumor",]$SRR_id#335
tumor_res=mouse_res[,colnames(mouse_res)%in%mouse_info_tumor]#335
tumor_res$Mouse.gene.stable.ID=rownames(tumor_res)
HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]#16550
tumor_res11=left_join(tumor_res,HOM_one2one,by="Mouse.gene.stable.ID")
tumor_res12=na.omit(tumor_res11)
rownames(tumor_res12)=tumor_res12$Gene.stable.ID
colnames( tumor_res12)
tumor_res2=tumor_res12[,-c(45,46,47)]
tumor_cms <- CMScaller(tumor_res2, RNAseq=TRUE, doPlot=TRUE, rowNames = "ensg", seed = 1)
mouse_cmstypes <- data.frame("barcode"=row.names(tumor_cms), "type"=tumor_cms$prediction)
mouse_cmstypes=na.omit(mouse_cmstypes)#305
#write.csv(mouse_cmstypes,file = "/home/panxl/CRC2/mouse_cmstypes.csv",row.names = F)
unique
cmstypes=read.csv("/home/panxl/CRC2/mouse_cmstypes.csv")

unique(cmstypes$type)
cmstypes$type <- as.vector(cmstypes$type)
cmstypes$type[is.na(cmstypes$type)] = "None"
cmstypes$type = factor(cmstypes$type)
table(cmstypes$type)
kk=as.data.frame(table(cmstypes$type))
names(kk) = c("type","number")
p1<-ggplot(data=kk, aes(x=type, y=number, fill = type)) +
  geom_bar(stat="identity", width=0.5)+
  labs(x="",y="",title="TCGA Human CMS Type")+
  theme_classic(base_line_size = 1) +
  theme(title = element_text(size = 22),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 20 ,# 修改X轴上字体大小，
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 20,  
                                   vjust = 0.4))
p1


