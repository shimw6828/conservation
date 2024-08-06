library('GenomicRanges') 
library('rtracklayer') 
library('dplyr') 
library(GenomicFeatures)
library(ChIPseeker)
bedfile = "/home/zhluo/Project/CRC_conservation/step13_chip_seq_DEG/all_narrow_peak_human/human_H3K27ac_master.sort.merge.bed"
data<-read.table(bedfile,header=F)
bigwig_dir <- "/home/zhluo/Project/CRC_conservation/step13_chip_seq_DEG/all_narrow_peak_human/bigwig_enhancer"
list.files(bigwig_dir,full.names = F)
pr_files = list.files(bigwig_dir,full.names = T)
file_name = list.files(bigwig_dir,full.names = F)%>%stringr::str_remove(".bigwig")
names(pr_files)=file_name
outputfname="/home/mwshi/project/conservation/enhancer_count/human_enhancer_all_counts.csv"
colnames(data)<-c("chr","start","end")
bed<-with(data,GRanges(chr,IRanges(start,end)))

hg19 <-  makeTxDbFromGFF("/home/mwshi/project/conservation/gencode.v32lift37.annotation.gtf", format = "gtf")
peakAnno <- annotatePeak(bed, tssRegion=c(-3000, 3000), TxDb=hg19)
peakAnno <- as.GRanges(peakAnno)
peakAnno$geneId = stringr::str_split_fixed(peakAnno$geneId,"\\.",2)[,1]
bw=BigWigFileList(pr_files) 
counts <- matrix(NA, nrow = length(bed), ncol = length(bw))
colnames(counts) <- names(bw)
chromnames=levels(seqnames(bed))

for(i in seq_len(length(bw))) 
{
  last=1 	 
  coverage <- import(bw[[i]], as = 'RleList')
  for(j in chromnames)
  {
    range_vals=ranges(bed[seqnames(bed)==j])
    cur_coverage=coverage[[j]] 
    if(is.null(cur_coverage))
    {
      counts[last:(last+length(range_vals)-1),i]=matrix(0,nrow=length(range_vals),ncol=1)
    }
    else
    {
      newvals=sum(Views(cur_coverage, ranges(bed[seqnames(bed)==j])))
      counts[last:(last+length(newvals)-1), i] <-newvals 
    }   
    last=last+length(range_vals) 
  }
  
}


counts <- round(counts / 36, 0)
row.names(counts) <- stringr::str_c(seqnames(peakAnno),start(peakAnno)-1,end(peakAnno),sep = "_")
write.csv(counts,file=outputfname,quote = F)
counts <- read.csv(outputfname,row.names = 1)



##############差异分析
library(DESeq2)
cData <- data.frame(sample = colnames(counts),
                    "group" = factor(ifelse(startsWith(colnames(counts),"N"),"Normal","Tumor"),
                                     levels = c("Tumor","Normal")))
rownames(cData) <- colnames(counts)
d.deseq <- DESeqDataSetFromMatrix(countData = counts, colData = cData,design = ~group)
d.deseq <- DESeq(d.deseq)
res <- results(d.deseq, contrast=c("group","Tumor","Normal"))
write.csv(res,"/home/mwshi/project/conservation/human_enhancer_DE.csv",quote = F)
####去掉header之后变成human_enhancer_DE.bed

###########################################################################
###mouse
bedfile = "/home/mwshi/project/conservation/NACC/enhancer.peak.unique.ID.gene_name.sorted.bed"
data<-read.table(bedfile,header=F)
data <- data[,c(1,2,3,4)]
bigwig_dir <- "/home/zhluo/Project/CRC/data_nazhang/step39_read_count/BamCompare/"
list.files(bigwig_dir,full.names = F, pattern = "*H3K27ac*")
pr_files = list.files(bigwig_dir,full.names = T, pattern = "*H3K27ac*")[c(1,2,3,10:15)]
# file_name = list.files(bigwig_dir,full.names = F)%>%stringr::str_remove(".bigwig")
names(pr_files)=c("T1","T2","T3","T4","T5","T6","N1","N2","N3")
outputfname="/home/mwshi/project/conservation/enhancer_count/mouse_enhancer_all_counts.csv"
colnames(data)<-c("chr","start","end","peakid")
bed<-with(data,GRanges(chr,IRanges(start,end)))

mm10 <-  makeTxDbFromGFF("/home/mwshi/project/conservation/gencode.vM20.basic.annotation.gtf", format = "gtf")
peakAnno <- annotatePeak(bed, tssRegion=c(-3000, 3000), TxDb=mm10)
peakAnno <- as.GRanges(peakAnno)
peakAnno$geneId = stringr::str_split_fixed(peakAnno$geneId,"\\.",2)[,1]
bw=BigWigFileList(pr_files) 
counts <- matrix(NA, nrow = length(bed), ncol = length(bw))
colnames(counts) <- names(bw)
chromnames=levels(seqnames(bed))

for(i in seq_len(length(bw))) 
{
  last=1 	 
  coverage <- import(bw[[i]], as = 'RleList')
  for(j in chromnames)
  {
    range_vals=ranges(bed[seqnames(bed)==j])
    cur_coverage=coverage[[j]] 
    if(is.null(cur_coverage))
    {
      counts[last:(last+length(range_vals)-1),i]=matrix(0,nrow=length(range_vals),ncol=1)
    }
    else
    {
      newvals=sum(Views(cur_coverage, ranges(bed[seqnames(bed)==j])))
      counts[last:(last+length(newvals)-1), i] <-newvals 
    }   
    last=last+length(range_vals) 
  }
  
}
.
counts <- round(counts / 36, 0)
row.names(counts) <- data$peakid
write.csv(counts,file=outputfname,quote = F)
counts <- read.csv("/home/mwshi/project/conservation/enhancer_count/mouse_enhancer_all_counts.csv",row.names = 1)
####差异
library(DESeq2)
cData <- data.frame(sample = colnames(counts),
                    "group" = factor(ifelse(startsWith(colnames(counts),"N"),"Normal","Tumor"),
                                     levels = c("Tumor","Normal")))
rownames(cData) <- colnames(counts)
d.deseq <- DESeqDataSetFromMatrix(countData = counts, colData = cData,design = ~group)
d.deseq <- DESeq(d.deseq)
res <- results(d.deseq, contrast=c("group","Tumor","Normal"))
write.csv(res,"/home/mwshi/project/conservation/mouse_enhancer_DE.csv",quote = F)
