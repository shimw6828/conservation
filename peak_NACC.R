####在生科院服务器上运行
library(dplyr)
human_peak_id <- readr::read_tsv("/home/mwshi/project/conservation/NACC/human_peak_id.bed",
                                 col_names = c("h_chr","h_start","h_end","id"))
human_peak_id <- dplyr::mutate(human_peak_id,"h_peak_id" = stringr::str_c(h_chr,h_start,h_end,sep = "_"))

mouse_peak_id <- readr::read_tsv("/home/mwshi/project/conservation/NACC/mouse_peak_id.bed",
                                 col_names = c("m_chr","m_start","m_end","m_peak_id","middle","gene_id","id"))

# intersect(unique(human_peak_id$id),unique(mouse_peak_id$id))
####id确实是唯一的，算过unique之后，数量和之前没有区别
instersect_id = intersect(human_peak_id$id,mouse_peak_id$id)

human_peak_id = human_peak_id[human_peak_id$id%in%instersect_id,]
mouse_peak_id = mouse_peak_id[mouse_peak_id$id%in%instersect_id,]

orth_peak <- dplyr::left_join(human_peak_id,mouse_peak_id)
kk = group_by(orth_peak,h_peak_id, m_peak_id)
nn = summarise(kk,a=n())
orth_peak <- nn %>%top_n(1, a)%>%group_by(m_peak_id)%>%top_n(1, a)
human_peak_id <- human_peak_id%>% dplyr::mutate(h_length = h_end - h_start)%>%
  dplyr::select(h_peak_id,h_length)%>%dplyr::distinct()
mouse_peak_id<-mouse_peak_id%>% dplyr::mutate(m_length = m_end - m_start)%>%
  dplyr::select(m_peak_id,m_length)%>%dplyr::distinct()
orth_peak
orth_peak<- left_join(orth_peak,human_peak_id)%>%left_join(mouse_peak_id)
orth_peak<- orth_peak%>%group_by(m_peak_id) %>%top_n(1, h_length)%>%group_by(h_peak_id)%>%top_n(1, m_length)
orth_peak = orth_peak[!orth_peak$m_peak_id%in%names(table(orth_peak$m_peak_id)[table(orth_peak$m_peak_id)>1]),]
orth_peak = orth_peak[!orth_peak$h_peak_id%in%names(table(orth_peak$h_peak_id)[table(orth_peak$h_peak_id)>1]),]




readr::write_csv(orth_peak,"/home/mwshi/project/conservation/NACC/peak_orth_map.csv")

# table(orth_peak$h_peak_id)[table(orth_peak$h_peak_id)>1]
human_peak_id <- human_peak_id%>%select(-id)%>%dplyr::distinct()
mouse_peak_id <- mouse_peak_id%>%select(-id,-middle,-gene_id)%>%dplyr::distinct()
orth_peak <- readr::read_csv("/home/mwshi/project/conservation/NACC/peak_orth_map.csv")
human_peak_id[match(orth_peak$h_peak_id,human_peak_id$h_peak_id),]
mouse_peak_id[match(orth_peak$m_peak_id,mouse_peak_id$m_peak_id),]

readr::write_tsv(human_peak_id[match(orth_peak$h_peak_id,human_peak_id$h_peak_id),-4],
                 "/home/mwshi/project/conservation/human_orthpeak.bed",col_names = F)

readr::write_tsv(mouse_peak_id[match(orth_peak$m_peak_id,mouse_peak_id$m_peak_id),-4],
                 "/home/mwshi/project/conservation/mouse_orthpeak.bed",col_names = F)

all = cbind(human_peak_id[match(orth_peak$h_peak_id,human_peak_id$h_peak_id),],mouse_peak_id[match(orth_peak$m_peak_id,mouse_peak_id$m_peak_id),])


#########################################################################
orth_peak <- readr::read_csv("/home/shimw/project/conservation/NACC/peak_orth_map.csv")
###有一个值是chr11_9.6e+07_96000845，手动改了一下
orth_peak[!orth_peak$h_peak_id%in%row.names(human_peak_count),]$h_peak_id = "chr11_96000000_96000845"
human_peak_count <- read.csv("/home/shimw/project/conservation/enhancer_count/human_enhancer_all_counts.csv",row.names = 1)
mouse_peak_count <- read.csv("/home/shimw/project/conservation/enhancer_count/mouse_enhancer_all_counts.csv",row.names = 1)

human_peak_count = human_peak_count[orth_peak$h_peak_id,]
mouse_peak_count = mouse_peak_count[orth_peak$m_peak_id,]

human_peak_count_t = t(human_peak_count)
mouse_peak_count_t = t(mouse_peak_count)

human_peak_pairs_cor = cor(human_peak_count_t)
mouse_peak_pairs_cor = cor(mouse_peak_count_t)

readr::write_rds(human_peak_pairs_cor,"/home/shimw/project/conservation/NACC/human_peak_pairs_cor.rds")
readr::write_rds(mouse_peak_pairs_cor,"/home/shimw/project/conservation/NACC/mouse_peak_pairs_cor.rds")



human_peak_top20_cor  <- purrr::map(colnames(human_peak_pairs_cor),function(m){
  print(m)
  x = human_peak_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  # names(x[2:21])
  tmp_orth <- orth_peak[match(c(m,names(x[2:21])),orth_peak$h_peak_id),]$m_peak_id
  kk = dist(rbind(x[2:21],mouse_peak_pairs_cor[tmp_orth[2:21],tmp_orth[1]]))
  tibble("gene"=m,
         "d"=kk[1])
})%>%bind_rows()
human_peak_top20_cor = human_peak_top20_cor[match(orth_peak$h_peak_id,human_peak_top20_cor$gene),]
readr::write_csv(human_peak_top20_cor,"/home/shimw/project/conservation/NACC/human_peak_distance.csv")

####出现了个别情况peak在所有样本中相同，cor为NA值，这里给0，也不影响后续计算，因为我们是选择的前20的peak
mouse_peak_pairs_cor[is.na(mouse_peak_pairs_cor)] = 0
mouse_peak_top20_cor  <- purrr::map(colnames(mouse_peak_pairs_cor),function(m){
  print(m)
  x = mouse_peak_pairs_cor[,m]
  
  x= sort(x, decreasing = TRUE)
  names(x[2:21])
  tmp_orth <- orth_peak[match(c(m,names(x[2:21])),orth_peak$m_peak_id),]$h_peak_id
  kk = dist(rbind(x[2:21],human_peak_pairs_cor[tmp_orth[2:21],tmp_orth[1]]))
  tibble("gene"=m,
         "d"=kk[1])
})%>%bind_rows()
mouse_peak_top20_cor = mouse_peak_top20_cor[match(orth_peak$m_peak_id,mouse_peak_top20_cor$gene),]
readr::write_csv(mouse_peak_top20_cor,"/home/shimw/project/conservation/NACC/mouse_peak_distance.csv")


human_peak_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/human_peak_distance.csv")
mouse_peak_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/mouse_peak_distance.csv")
human_peak_top20_cor$score = (human_peak_top20_cor$d + mouse_peak_top20_cor$d)/2
hist(human_peak_top20_cor$score)





###########################
###随机
human_peak_pairs_cor <- readr::read_rds("/home/shimw/project/conservation/NACC/human_peak_pairs_cor.rds")
mouse_peak_pairs_cor <-readr::read_rds("/home/shimw/project/conservation/NACC/mouse_peak_pairs_cor.rds")
orth_peak <- readr::read_csv("/home/shimw/project/conservation/NACC/peak_orth_map.csv")
###有一个值是chr11_9.6e+07_96000845，手动改了一下
orth_peak[orth_peak$h_peak_id=="chr11_9.6e+07_96000845",]$h_peak_id = "chr11_96000000_96000845"


set.seed(12345)
human_peak_random_cor  <- purrr::map(colnames(human_peak_pairs_cor),function(m){
  x = human_peak_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  # names(x[2:21])
  rand_num = sample(1:nrow(mouse_peak_pairs_cor),20,replace = T)
  tmp_orth <- orth_peak[orth_peak$h_peak_id==m,]$m_peak_id
  kk = dist(rbind(x[2:21],mouse_peak_pairs_cor[rand_num,tmp_orth[1]]))
  tibble("gene"=m,
         "d"=kk[1])
})%>%bind_rows()
human_peak_random_cor = human_peak_random_cor[match(orth_peak$h_peak_id,human_peak_random_cor$gene),]
readr::write_csv(human_peak_random_cor,"/home/shimw/project/conservation/NACC/human_peak_random_distance.csv")
####出现了个别情况peak在所有样本中相同，cor为NA值，这里给0，也不影响后续计算，因为我们是选择的前20的peak
mouse_peak_pairs_cor[is.na(mouse_peak_pairs_cor)] = 0
mouse_peak_random_cor  <- purrr::map(colnames(mouse_peak_pairs_cor),function(m){
  x = mouse_peak_pairs_cor[,m]
  x= sort(x, decreasing = TRUE)
  # names(x[2:21])
  rand_num = sample(1:nrow(human_peak_pairs_cor),20,replace = T)
  tmp_orth <- orth_peak[orth_peak$m_peak_id==m,]$h_peak_id
  kk = dist(rbind(x[2:21],human_peak_pairs_cor[rand_num,tmp_orth[1]]))
  tibble("gene"=m,
         "d"=kk[1])
})%>%bind_rows()
mouse_peak_random_cor = mouse_peak_random_cor[match(orth_peak$m_peak_id,mouse_peak_random_cor$gene),]
readr::write_csv(mouse_peak_random_cor,"/home/shimw/project/conservation/NACC/mouse_peak_random_distance.csv")




############################################
###画分布图
human_peak_random_cor = readr::read_csv("/home/shimw/project/conservation/NACC/human_peak_random_distance.csv")
mouse_peak_random_cor = readr::read_csv("/home/shimw/project/conservation/NACC/mouse_peak_random_distance.csv")
human_peak_random_cor$score = (human_peak_random_cor$d + mouse_peak_random_cor$d)/2
human_peak_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/human_peak_distance.csv")
mouse_peak_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/mouse_peak_distance.csv")
human_peak_top20_cor$score = (human_peak_top20_cor$d + mouse_peak_top20_cor$d)/2

human_peak_top20_cor$type = "true"
human_peak_random_cor$type = "random"

human_peak_compare = rbind(human_peak_top20_cor,human_peak_random_cor)
library(ggplot2)
p <- ggplot(human_peak_compare,aes(x = score,fill = type)) + 
  geom_histogram(aes(y=..density..),colour="black",position = "identity",alpha = 0.6,bins = 50)+
  annotate("text", x= c(1.2,5.2) , y= c(0.7,0.6) ,
           label= c(stringr::str_wrap("Orthologous neighbourhood enhancers",width = 10),
                    stringr::str_wrap("Random neighbourhood enhancers",width = 10)),
           size = 5, lineheight=0.85)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 0.5)+
  scale_fill_manual(values=c("true" = "#d1706d", "random" = "#8ea9cf"))+
  labs(y="Density", x="Distance", fill=NULL)+
  theme(legend.position = "none",
        axis.title=element_text(size=18),
        axis.text=element_text(size=15),
        plot.margin = margin(0.6, 0.1, 0.1, 0.1, "cm"))+
  coord_cartesian(clip = "off")














