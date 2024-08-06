####调整bed文件的基因顺序，方便计算computeMatrix

####human
human_data <- read.csv("/home/shimw/project/conservation/filter_human_data.csv",row.names = 1, check.names = F)
human_info <- read.csv("/home/panxl/CRC2/filter_human_info.csv",row.names = 1, check.names = F)

human_info[human_info$Group=="Tumor",]$sample
human_expression_level <- apply(human_data[,human_info[human_info$Group=="Tumor",]$sample],1,mean)%>%
  sort(decreasing = TRUE)
bed_file = readr::read_tsv("/home/shimw/hg19.bed",col_names = c("chr","start","end","gene"))
bed_file$gene = stringr::str_split_fixed(bed_file$gene,"\\.",2)[,1]
# bed_file$gene%in%names(human_expression_level)
####只需要这部分的基因
human_expression_level = human_expression_level[names(human_expression_level)[names(human_expression_level)%in%bed_file$gene]]
bed_file <- bed_file[match(names(human_expression_level),bed_file$gene),]
out_chr <- c("GL000241.1", "GL000192.1", "GL000230.1", "GL000220.1", "GL000237.1", "GL000199.1", "GL000193.1",
  "GL000212.1", "GL000202.1", "GL000228.1", "GL000195.1", "GL000205.1", "GL000204.1","chrM")
bed_file <- bed_file[!bed_file$chr%in%out_chr,]
readr::write_tsv(bed_file,"/home/shimw/project/conservation/hg19_rank.bed",col_names =F)








mouse_data <- read.csv("/home/shimw/project/conservation/filter_mouse_data.csv",row.names = 1, check.names = F)
mouse_info <- read.csv("/home/panxl/CRC2/filter_mouse_info.csv",row.names = 1, check.names = F)
