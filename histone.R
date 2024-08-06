#####密度图
#################################################################
##数据位置在/home/shimw/project/conservation/histone
##首先做测试
read_density_data <- function(table_path){
  input_species = strsplit(table_path,"/")[[1]][7]
  input_cluster = strsplit(strsplit(table_path,"/")[[1]][8],"_")[[1]][1]
  df = read.delim(table_path, sep="\t", skip=2, header=FALSE, check.names = FALSE,stringsAsFactors =FALSE)
  df = df[, -c(1,2)]
  rownames(df) = c("Normal", "Tumor")
  df <- mutate_all(df, function(x) as.numeric(x))
  colnames(df) = 1:length(colnames(df))
  df = t(df)%>%as.data.frame()
  df$bin = as.numeric(row.names(df))
  long_table <- tidyr::gather(df, 'Group', 'score', -bin)
  long_table$species = input_species
  long_table$cluster = input_cluster
  return(long_table)
}


human_tables=list.files("/home/shimw/project/conservation/histone/human",full.names = T,pattern ="H3K27ac")
human_table_result = purrr::map(human_tables,read_density_data)%>%dplyr::bind_rows()
mouse_tables=list.files("/home/shimw/project/conservation/histone/mouse",full.names = T,pattern ="H3K27ac")
mouse_table_result = purrr::map(mouse_tables,read_density_data)%>%dplyr::bind_rows()
all_table_result=dplyr::bind_rows(human_table_result,mouse_table_result)

library("ggplot2")
library("ggsci")
library(lemon)
p1 = ggplot(all_table_result,aes(x=bin,y=score, group=Group, color=Group)) +
  #geom_line() +
  geom_smooth(stat = "smooth",method = "loess",formula = y ~ x,se = FALSE,span=0.01,size=0.6)  +
  theme_classic(base_family = "sans",base_size = 14) +
  scale_colour_manual(values = c("#259aff","#ff6e3a"))+
  scale_x_continuous(breaks = seq(0,2000,by=1000), labels=c("-10k", "TSS", "10k"))+
  ylab("H3K27ac intensity")+xlab("Gene promoter region") +
  # facet_grid(species ~ cluster,scales = "free")+
  facet_rep_grid(species~ cluster, scales='free')+
  theme(strip.background = element_blank(),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3),
        strip.text.y = element_blank())
  ##上面是human，下面是mouse

ggsave("/home/shimw/project/conservation/H3K27ac_histone.pdf",p1,width =8 ,height = 4)

################################################################################
####更换不同的marker画图
markers = c("H3K27ac","H3K27me3","H3K4me1","H3K4me3","H3K9me3")
marker = "H3K9me3"
for (marker in markers) {
  tmp_human_tables=list.files("/home/shimw/project/conservation/histone/human",full.names = T,pattern =marker)
  tmp_human_table_result = purrr::map(tmp_human_tables,read_density_data)%>%dplyr::bind_rows()
  tmp_mouse_tables=list.files("/home/shimw/project/conservation/histone/mouse",full.names = T,pattern =marker)
  tmp_mouse_table_result = purrr::map(tmp_mouse_tables,read_density_data)%>%dplyr::bind_rows()
  tmp_all_table_result=dplyr::bind_rows(tmp_human_table_result,tmp_mouse_table_result)
  p1 = ggplot(tmp_all_table_result,aes(x=bin,y=score, group=Group, color=Group)) +
    #geom_line() +
    geom_smooth(stat = "smooth",method = "loess",formula = y ~ x,se = FALSE,span=0.01,size=0.6)  +
    theme_classic(base_family = "sans",base_size = 14) +
    scale_colour_manual(values = c("#259aff","#ff6e3a"))+
    scale_x_continuous(breaks = seq(0,2000,by=1000), labels=c("-10", "TSS", "10"))+
    ylab(paste(marker,"intensity"))+xlab("Distance from promoter(kb)") +
    # facet_grid(species ~ cluster,scales = "free")+
    facet_rep_grid(species~ cluster, scales='free')+
    theme(strip.background = element_blank(),
          axis.line=element_line(size=0.3),
          axis.ticks=element_line(size=0.3),
          strip.text.y = element_blank())
  
  ##上面是human，下面是mouse
  ggsave(paste("/home/shimw/project/conservation/",marker,"_histone.pdf",sep = ""),p1,width =8 ,height = 4)
  
}







