rm(list=ls())

human_data <- read.csv("/home/shimw/project/conservation/filter_human_data.csv",row.names = 1, check.names = F)
human_info <- read.csv("/home/panxl/CRC2/filter_human_info.csv",row.names = 1, check.names = F)
mouse_data <- read.csv("/home/shimw/project/conservation/filter_mouse_data.csv",row.names = 1, check.names = F)
mouse_info <- read.csv("/home/panxl/CRC2/filter_mouse_info.csv",row.names = 1, check.names = F)
#cluster_type = read.csv("/home/panxl/CRC2/all_cluater",check.names = F)
#cluster_type = cluster_type[cluster_type$cluster!="Cluster 5: Not-significant in both species",]
#human_data <- human_data[cluster_type$Gene.stable.ID,]%>%t()%>%
#  scale()%>%t()%>%as.data.frame()
#mouse_data <- mouse_data[cluster_type$Mouse.gene.stable.ID,]%>%t()%>%
#  scale()%>%t()%>%as.data.frame()
#row.names(mouse_data) <- cluster_type$Gene.stable.ID

HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]

HOM_sub = HOM_one2one[HOM_one2one$Gene.stable.ID %in% rownames(human_data),]
HOM_sub = HOM_sub[HOM_sub$Mouse.gene.stable.ID %in% rownames(mouse_data),]
readr::write_rds(HOM_sub,"/home/shimw/project/conservation/HOM_sub.rds")


human_ortho = human_data[rownames(human_data) %in% HOM_sub$Gene.stable.ID,]
nrow(human_ortho)
mouse_ortho = mouse_data[rownames(mouse_data) %in% HOM_sub$Mouse.gene.stable.ID,]
nrow(mouse_ortho)

human_data_t = t(human_ortho)
mouse_data_t = t(mouse_ortho)
rm("human_ortho")
rm("mouse_ortho")
rm(list = c("human_ortho","mouse_ortho","human_data","mouse_data","human_info",
     "mouse_info","HOM_one2one"))

human_pairs_cor = cor(human_data_t)
mouse_pairs_cor = cor(mouse_data_t)
readr::write_rds(human_pairs_cor,"/home/shimw/project/conservation/NACC/human_gene_pairs_cor.rds")
readr::write_rds(mouse_pairs_cor,"/home/shimw/project/conservation/NACC/mouse_gene_pairs_cor.rds")



library(dplyr)
human_pairs_top20_cor <- purrr::map(colnames(human_pairs_cor),function(m){
  x = human_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  names(x[2:21])
  tmp_orth <- HOM_sub[match(c(m,names(x[2:21])),HOM_sub$Gene.stable.ID),]$Mouse.gene.stable.ID
  kk = dist(rbind(x[2:21],mouse_pairs_cor[tmp_orth[2:21],tmp_orth[1]]))
  tibble("gene"=m,
         "d"=kk[1])
})%>%bind_rows()
human_pairs_top20_cor = human_pairs_top20_cor[match(HOM_sub$Gene.stable.ID,human_pairs_top20_cor$gene),]
mouse_pairs_top20_cor <- purrr::map(colnames(mouse_pairs_cor),function(m){
  x = mouse_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  names(x[2:21])
  tmp_orth <- HOM_sub[match(c(m,names(x[2:21])),HOM_sub$Mouse.gene.stable.ID),]$Gene.stable.ID
  kk = dist(rbind(x[2:21],human_pairs_cor[tmp_orth[2:21],tmp_orth[1]]))
  tibble("gene"=m,
         "d"=kk[1])
})%>%bind_rows()

human_pairs_top20_cor = human_pairs_top20_cor[match(HOM_sub$Gene.stable.ID,human_pairs_top20_cor$gene),]
mouse_pairs_top20_cor = mouse_pairs_top20_cor[match(HOM_sub$Mouse.gene.stable.ID,mouse_pairs_top20_cor$gene),]
readr::write_csv(human_pairs_top20_cor,"/home/shimw/project/conservation/NACC/human_pairs_distance.csv")
readr::write_csv(mouse_pairs_top20_cor,"/home/shimw/project/conservation/NACC/mouse_pairs_distance.csv")

human_pairs_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/human_pairs_distance.csv")
mouse_pairs_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/mouse_pairs_distance.csv")


human_pairs_top20_cor$dd = (human_pairs_top20_cor$d + mouse_pairs_top20_cor$d)/2

hist(human_pairs_top20_cor$dd)




###########################
###随机
HOM_sub = readr::read_rds("/home/shimw/project/conservation/HOM_sub.rds")
human_gene_pairs_cor <- readr::read_rds("/home/shimw/project/conservation/NACC/human_gene_pairs_cor.rds")
mouse_gene_pairs_cor <-readr::read_rds("/home/shimw/project/conservation/NACC/mouse_gene_pairs_cor.rds")


set.seed(12345)
human_gene_random_cor  <- purrr::map(colnames(human_gene_pairs_cor),function(m){
  x = human_gene_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  # names(x[2:21])
  rand_num = sample(1:nrow(mouse_gene_pairs_cor),20,replace = T)
  tmp_orth <- HOM_sub[HOM_sub$Gene.stable.ID==m,]$Mouse.gene.stable.ID
  kk = dist(rbind(x[2:21],mouse_gene_pairs_cor[rand_num,tmp_orth[1]]))
  # kk = dist(rbind(x[sample(1:nrow(human_gene_pairs_cor),20,replace = T)],mouse_gene_pairs_cor[rand_num,tmp_orth[1]]))
  tibble("gene"=m,
         "d"=kk[1])
})%>%bind_rows()
# readr::write_csv(human_gene_random_cor,"/home/shimw/project/conservation/NACC/human_gene_random_distance.csv")
human_gene_random_cor = human_gene_random_cor[match(HOM_sub$Gene.stable.ID,human_gene_random_cor$gene),]
set.seed(54321)
####出现了个别情况peak在所有样本中相同，cor为NA值，这里给0，也不影响后续计算，因为我们是选择的前20的peak
mouse_gene_random_cor  <- purrr::map(colnames(mouse_gene_pairs_cor),function(m){
  x = mouse_gene_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  # names(x[2:21])
  rand_num = sample(1:nrow(human_gene_pairs_cor),20,replace = T)
  tmp_orth <- HOM_sub[HOM_sub$Mouse.gene.stable.ID==m,]$Gene.stable.ID
  kk = dist(rbind(x[2:21],human_gene_pairs_cor[rand_num,tmp_orth[1]]))
  # kk = dist(rbind(x[sample(1:nrow(mouse_gene_pairs_cor),20,replace = T)],human_gene_pairs_cor[rand_num,tmp_orth[1]]))
  tibble("gene"=m,
         "d"=kk[1])
})%>%bind_rows()
# readr::write_csv(mouse_gene_random_cor,"/home/shimw/project/conservation/NACC/mouse_gene_random_distance.csv")
mouse_gene_random_cor = mouse_gene_random_cor[match(HOM_sub$Mouse.gene.stable.ID,mouse_gene_random_cor$gene),]

readr::write_csv(human_gene_random_cor,"/home/shimw/project/conservation/NACC/human_gene_random_distance.csv")
readr::write_csv(mouse_gene_random_cor,"/home/shimw/project/conservation/NACC/mouse_gene_random_distance.csv")

human_gene_random_cor = readr::read_csv("/home/shimw/project/conservation/NACC/human_gene_random_distance.csv")
mouse_gene_random_cor = readr::read_csv("/home/shimw/project/conservation/NACC/mouse_gene_random_distance.csv")

human_gene_random_cor$score = (human_gene_random_cor$d + mouse_gene_random_cor$d)/2
hist(human_gene_random_cor$score)



############################################
###直方图
human_pairs_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/human_pairs_distance.csv")
mouse_pairs_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/mouse_pairs_distance.csv")
human_gene_random_cor = readr::read_csv("/home/shimw/project/conservation/NACC/human_gene_random_distance.csv")
mouse_gene_random_cor = readr::read_csv("/home/shimw/project/conservation/NACC/mouse_gene_random_distance.csv")

human_pairs_top20_cor$score = (human_pairs_top20_cor$d + mouse_pairs_top20_cor$d)/2
human_pairs_top20_cor$type= "true"
human_gene_random_cor$score = (human_gene_random_cor$d + mouse_gene_random_cor$d)/2
human_gene_random_cor$type= "random"
human_gene_compare <- rbind(human_pairs_top20_cor,human_gene_random_cor)
library(ggthemes)


stringr::str_wrap("Orthologous neighbourhood genes",width = 10)

# ggplot(human_gene_compare,aes(x = score,fill = type, color = type,)) + geom_density(alpha = 0.3,adjust = 1/2)
p <- ggplot(human_gene_compare,aes(x = score,fill = type)) + 
  geom_histogram(aes(y=..density..),colour="black",position = "identity",alpha = 0.6,bins = 50)+
  annotate("text", x= c(1.2,5.2) , y= c(0.7,0.6) ,
           label= c(stringr::str_wrap("Orthologous neighbourhood genes",width = 10),
                    stringr::str_wrap("Random neighbourhood genes",width = 10)),
           size = 5, lineheight=0.85)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 0.5)+
  scale_fill_manual(values=c("true" = "#d1706d", "random" = "#8ea9cf"))+
  labs(y="Density", x="Distance", fill=NULL)+
  theme(legend.position = "none",
        axis.title=element_text(size=18),
        axis.text=element_text(size=15),
        plot.margin = margin(0.6, 0.1, 0.1, 0.1, "cm"))+
  coord_cartesian(clip = "off")
  
ggsave("/home/shimw/project/conservation/NACC/gene_random_compare.pdf",p,width = 4.5,height = 3.5)


########################################################################
#####四个cluster的NACC的分布图
human_pairs_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/human_pairs_distance.csv")
mouse_pairs_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/mouse_pairs_distance.csv")
human_pairs_top20_cor$score = (human_pairs_top20_cor$d + mouse_pairs_top20_cor$d)/2
df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
df = dplyr::left_join(human_pairs_top20_cor, df_data, by = c("gene"="Gene.stable.ID"))
my_comparisons <- list(c("Cluster 1", "Cluster 5"), c("Cluster 2", "Cluster 5"), c("Cluster 3", "Cluster 5"), c("Cluster 4", "Cluster 5"))

p1 <- ggplot(data = na.omit(df),aes(x = cluster,y = score,fill=cluster)) +
  geom_boxplot(width = 0.6)+
  stat_compare_means(comparisons = my_comparisons,label = "p.format", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 0.6)+
  scale_fill_manual(values=c("Cluster 1" = "#CB9C7A", "Cluster 2" = "#8696a7", "Cluster 3" = "#CDB97D",
                             "Cluster 4" = "#7b8b6f", "Cluster 5" = "#A59B95"))+
  labs( y="Phastcons score", x=NULL, fill=NULL) +
  scale_x_discrete(labels = c("Cluster 1" = "Cluster1", "Cluster 2" = "Cluster2", "Cluster 3" = "Cluster3",
                              "Cluster 4" = "Cluster4", "Cluster 5" = "Cluster5"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.9))







aa =boxplot(score ~ cluster,data=na.omit(df), 
            col= c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"), 
            # las = 2,
            ylab="NACC score",xlab = NULL,outline = F,xaxt = "n",cex.lab=1.3,ylim = c(0,6.5))

Cairo::CairoPDF("/home/shimw/project/conservation/NACC/nacc_cluster_boxplot.pdf",width = 3.6, height=4.4)
par(mar=c(4, 4.5, 1, 1))
boxplot(score ~ cluster,data=na.omit(df), 
        col= c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"), 
        # las = 2,
        ylab="NACC score",xlab = NULL,outline = F,xaxt = "n",cex.lab=1.3,ylim = c(0,6.2))
tick <- seq_along(df$cluster)
axis(1, at = 1:5, labels = F)
text( 1:5 - 0.25, par("usr")[3] - 0.7, c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5"), srt = 45, xpd = T)
text( 1:4, aa$stats[5,1:4] + 0.3, c("N.S.","p<2.2e-16","N.S.","p<2.2e-16"))
dev.off()



#######################################################
####换个计算方法，使用转录本的相似性计算聚类
HOM_sub = readr::read_rds("/home/shimw/project/conservation/HOM_sub.rds")
human_data_t = human_data[rownames(human_data) %in% HOM_sub$Gene.stable.ID,]%>%t()
mouse_data_t = mouse_data[rownames(mouse_data) %in% HOM_sub$Mouse.gene.stable.ID,]%>%t()
human_gene_pairs_cor <- readr::read_rds("/home/shimw/project/conservation/NACC/human_gene_pairs_cor.rds")
mouse_gene_pairs_cor <-readr::read_rds("/home/shimw/project/conservation/NACC/mouse_gene_pairs_cor.rds")
rm("human_data","human_info","mouse_data","mouse_info","HOM_one2one")


human_pairs_exp_distance <- purrr::map(colnames(human_gene_pairs_cor),function(m){
  x = human_gene_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  tmp_orth <- HOM_sub[match(c(m,names(x[2:21])),HOM_sub$Gene.stable.ID),]$Mouse.gene.stable.ID
  mean_d <- purrr::map(tmp_orth[-1],function(k){
    dd = dist(rbind(mouse_data_t[,tmp_orth[1]],mouse_data_t[,k]))
    return(dd)
  })%>%unlist()%>%mean()
  tibble("gene"=m,
         "d"=mean_d)
})%>%bind_rows()
human_pairs_exp_distance = human_pairs_exp_distance[match(HOM_sub$Gene.stable.ID,human_pairs_exp_distance$gene),]
readr::write_csv(human_pairs_exp_distance,"/home/shimw/project/conservation/NACC/human_pairs_exp_distance.csv")


mouse_pairs_exp_distance <- purrr::map(colnames(mouse_gene_pairs_cor),function(m){
  x = mouse_gene_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  tmp_orth <- HOM_sub[match(c(m,names(x[2:21])),HOM_sub$Mouse.gene.stable.ID),]$Gene.stable.ID
  mean_d <- purrr::map(tmp_orth[-1],function(k){
    dd = dist(rbind(human_data_t[,tmp_orth[1]],human_data_t[,k]))
    return(dd)
  })%>%unlist()%>%mean()
  tibble("gene"=m,
         "d"=mean_d)
})%>%bind_rows()
mouse_pairs_exp_distance = mouse_pairs_exp_distance[match(HOM_sub$Mouse.gene.stable.ID,mouse_pairs_exp_distance$gene),]
readr::write_csv(mouse_pairs_exp_distance,"/home/shimw/project/conservation/NACC/mouse_pairs_exp_distance.csv")

human_pairs_exp_distance$score = (human_pairs_exp_distance$d + mouse_pairs_exp_distance$d)/2
hist(human_pairs_exp_distance$score)


##########################################################################
###随机
human_random_exp_distance <- purrr::map(colnames(human_gene_pairs_cor),function(m){
  x = human_gene_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  rand_num = sample(1:nrow(mouse_gene_pairs_cor),20,replace = T)
  tmp_orth <- HOM_sub[HOM_sub$Gene.stable.ID==m,]$Mouse.gene.stable.ID
  mean_d <- purrr::map(rand_num,function(k){
    dd = dist(rbind(mouse_data_t[,tmp_orth[1]],mouse_data_t[,k]))
    return(dd)
  })%>%unlist()%>%mean()
  tibble("gene"=m,
         "d"=mean_d)
})%>%bind_rows()
human_random_exp_distance = human_random_exp_distance[match(HOM_sub$Gene.stable.ID,human_random_exp_distance$gene),]
readr::write_csv(human_random_exp_distance,"/home/shimw/project/conservation/NACC/human_random_exp_distance.csv")


mouse_random_exp_distance <- purrr::map(colnames(mouse_gene_pairs_cor),function(m){
  x = mouse_gene_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  rand_num = sample(1:nrow(human_gene_pairs_cor),20,replace = T)
  tmp_orth <- HOM_sub[HOM_sub$Mouse.gene.stable.ID==m,]$Gene.stable.ID
  mean_d <- purrr::map(tmp_orth[-1],function(k){
    dd = dist(rbind(human_data_t[,tmp_orth[1]],human_data_t[,k]))
    return(dd)
  })%>%unlist()%>%mean()
  tibble("gene"=m,
         "d"=mean_d)
})%>%bind_rows()
mouse_random_exp_distance = mouse_random_exp_distance[match(HOM_sub$Mouse.gene.stable.ID,mouse_random_exp_distance$gene),]
readr::write_csv(mouse_random_exp_distance,"/home/shimw/project/conservation/NACC/mouse_random_exp_distance.csv")

human_random_exp_distance$score = (human_random_exp_distance$d + mouse_random_exp_distance$d)/2
hist(human_pairs_exp_distance$score)

df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
df = left_join(df_data,human_pairs_top20_cor, by = c("Gene.stable.ID" = "gene"))
summary(df[df$cluster=="Cluster 1",]$dd)
summary(df[df$cluster=="Cluster 2",]$dd)
summary(df[df$cluster=="Cluster 3",]$dd)
summary(df[df$cluster=="Cluster 4",]$dd)
summary(df[df$cluster=="Cluster 5",]$dd)

df_data <- readr::read_csv("/home/shimw/project/conservation/phastcons/gene_to_tx.csv")
exon_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/exon_scores.csv")
exon_scores = left_join(exon_scores,df_data)

exon_scores[match(human_pairs_top20_cor$gene,exon_scores$Gene.stable.ID),]$score
cor.test(exon_scores[match(human_pairs_top20_cor$gene,exon_scores$Gene.stable.ID),]$score,
         human_pairs_top20_cor$dd)


tss_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/tss_scores.csv")
tss_scores = left_join(tss_scores,df_data)
cor.test(tss_scores[match(human_pairs_top20_cor$gene,tss_scores$Gene.stable.ID),]$score,
         human_pairs_top20_cor$dd)




library(phastCons100way.UCSC.hg19)
phast <- phastCons100way.UCSC.hg19

calc_score <- function(ob_grange){
  kk = c()
  for (x in names(ob_grange)) {
    print(x)
    y=ob_grange[[x]]
    if (length(y)==0) {
      kk = c(kk,NA)
      next
    }
    tmp_y = gscores(phast,y,summaryFun = sum)
    tmp_y = tmp_y[!is.na(tmp_y$default),]
    tmp_score <- sum(tmp_y$default)/sum(width(tmp_y))
    kk = c(kk,tmp_score)
  }
  kk = unlist(kk)
  return(tibble::tibble("tx_id"=names(ob_grange),"score"=kk))
}

txdb = makeTxDbFromGFF("/NAS/luozh/CRC_conservation/reference/gencode.v32lift37.annotation.gtf", format="gtf")
all_transcript = transcriptsBy(txdb,by = "gene")
genetotx = purrr::map(names(all_transcript),function(x){
  y = all_transcript[[x]]
  tx_id = y[order(width(y),decreasing = TRUE)[1],]$tx_name
  return(tibble::tibble("gene"=x,"tx_id" = tx_id))
})%>%dplyr::bind_rows()
genetotx$gene = gsub("\\..*", "", genetotx$gene)
readr::write_csv(genetotx,"/home/shimw/project/conservation/phastcons/human_genetotx.csv")
# genetotx <- readr::read_csv("/home/shimw/project/conservation/phastcons/genetotx.csv")
genetotx <- readr::read_csv("/home/shimw/project/conservation/phastcons/human_genetotx.csv")
df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data = df_data[,c("Gene.stable.ID","Mouse.gene.stable.ID","cluster")]
names(df_data)[1] = "gene"
df_data = dplyr::left_join(df_data,genetotx)
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
all_exon = exonsBy(txdb, by="tx", use.names = T)
all_exon = all_exon[df_data$tx_id]
exon_scores <- calc_score(all_exon)
exon_scores = left_join(exon_scores,df_data)


## foldchange deviation
df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
df_data$deviation = abs(df_data$Human_logFC1 - df_data$Mouse_logFC1)

aa = boxplot(deviation ~ cluster,data=na.omit(df_data), 
             col= c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"), 
             # las = 2,
             ylab="Absolute deviation in log2FC",xlab = NULL,outline = F,xaxt = "n",cex.lab=1.3,ylim = c(0,6.2))


Cairo::CairoPDF("/home/shimw/project/conservation/NACC/foldchange_deviation_boxplot.pdf",width = 3.6, height=4.4)
par(mar=c(4, 4.5, 1, 1))
boxplot(deviation ~ cluster,data=na.omit(df_data), 
        col= c("#CB9C7A","#8696a7","#CDB97D","#7b8b6f","#A59B95"), 
        # las = 2,
        ylab="Absolute deviation in log2FC",xlab = NULL,outline = F,xaxt = "n",cex.lab=1.3,ylim = c(0,6.2))
tick <- seq_along(df_data$cluster)
axis(1, at = 1:5, labels = F)
text( 1:5 - 0.25, par("usr")[3] - 0.7, c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5"), srt = 45, xpd = T)
text(c(2,4), aa$stats[5,c(2,4)] + 0.3, c("p<2.2e-16","p<2.2e-16"))
dev.off()


df_data[df_data$cluster=="Cluster 2",]$deviation
df_data[df_data$cluster%in%c("Cluster 1", "Cluster 3","Cluster 5"),]$deviation
t.test(df_data[df_data$cluster=="Cluster 4",]$deviation, df_data[df_data$cluster%in%c("Cluster 1", "Cluster 3","Cluster 5"),]$deviation)






