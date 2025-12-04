featureCounts -T 16 -p -B -C -s 2 \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Repeatmasker/T2T_CHM13v2_hs1_repeatmasker.gtf \
  -F GTF -t exon -g gene_id \
  -o TE_counts_gft_unique_12_3_2025.txt \
  case/hisat2_t2t_bam/*.clean.bam \
  control/hisat2_t2t_bam/*.clean.bam




featureCounts -T 16 -p -B -C -s 2 \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Repeatmasker/T2T_CHM13v2_hs1_repeatmasker.gtf \
  -F GTF -t exon -g gene_id \
  -M --fraction \
  -o TE_counts_gft_multi_fraction_12_3_2025.txt \
  case/hisat2_t2t_bam/*.clean.bam \
  control/hisat2_t2t_bam/*.clean.bam




  featureCounts -T 16 -p -B -C -s 2 \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Repeatmasker/T2T_CHM13v2_repeatmasker_family.saf \
  -F SAF \
  -M --fraction \
  -o T2T_CHM13v2_repeatmasker_TE_family_fraction_counts_all.txt \
  /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/case/hisat2_t2t_bam/*.clean.bam \
  /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/control/hisat2_t2t_bam/*.clean.bam
