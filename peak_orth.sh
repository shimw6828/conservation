bedtools intersect -sorted -wa -wb -a /home/zhluo/Project/CRC_conservation/step13_chip_seq_DEG/all_narrow_peak_human/human_H3K27ac_master.sort.merge.bed -b /home/zhluo/Project/CRC_conservation/step10_human_data/hs_mm_alignment/aggregate_base_output/hg19.mm10.50bp.h.gz | cut -f 1-3,7 > /home/mwshi/project/conservation/NACC/human_peak_id.bed


###enhancer.peak.unique.ID.gene_name.bed文件我自己手动去掉了第一行，vim
bedtools sort -i /home/mwshi/project/conservation/NACC/enhancer.peak.unique.ID.gene_name.bed > /home/mwshi/project/conservation/NACC/enhancer.peak.unique.ID.gene_name.sorted.bed
bedtools intersect -sorted -wa -wb -a /home/mwshi/project/conservation/NACC/enhancer.peak.unique.ID.gene_name.sorted.bed -b /home/zhluo/Project/CRC_conservation/step10_human_data/hs_mm_alignment/aggregate_base_output/hg19.mm10.50bp.m.gz | cut -f 1-6,10 > /home/mwshi/project/conservation/NACC/mouse_peak_id.bed



