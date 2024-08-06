/home/zhluo/miniconda3/envs/deeptools/bin/computeMatrix reference-point \
--referencePoint TSS \
-b 3000 -a 3000 \
--sortRegions keep \
--missingDataAsZero \
--skipZeros \
-R /home/mwshi/project/conservation/hg19_rank.bed \
-S /home/zhluo/Project/CRC_conservation/step12_intensity_comparison/human_bigwig/tumor_H3K4me3.bigwig \
	/home/zhluo/Project/CRC_conservation/step12_intensity_comparison/human_bigwig/tumor_H3K27ac.bigwig \
	/home/zhluo/Project/CRC_conservation/step12_intensity_comparison/human_bigwig/tumor_H3K4me1.bigwig \
	/home/zhluo/Project/CRC_conservation/step12_intensity_comparison/human_bigwig/tumor_H3K27me3.bigwig \
	/home/zhluo/Project/CRC_conservation/step12_intensity_comparison/human_bigwig/tumor_H3K9me3.bigwig \
-o /home/mwshi/project/conservation/human_tumor_all.mat.gz
	
	

/home/zhluo/miniconda3/envs/deeptools/bin/plotHeatmap -m /home/mwshi/project/conservation/human_tumor_all.mat.gz \
      -out /home/mwshi/project/conservation/human_tumor_all.png \
	  --sortRegions no \
	  --colorList 'white, red' 'white,red' 'white,red' 'white,red' 'white,red' \
	  --whatToShow "heatmap and colorbar" \
	  --zMax 4 4 3 2.5 2 \
	  --samplesLabel "H3K4me3" "H3K27ac" "H3K4me1" "H3K27me3" "H3K9me3" \
	  --xAxisLabel 