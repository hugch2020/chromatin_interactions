# Identify significant chromatin interactions
  This method is aimed to identify significant chromatin loops among the specified chromatin elements (e.g. CTCF, promoters, enhancers). Files used in this code include loop file (file generated by HiCcompare, or any file with format: chr1, start1, end1, chr2, start2, end2, sample1 count, sample2 count), element regions (e.g. CTCF binding sites, promoter regions, enhancer regions) and segmented genome file at specified resolution (e.g. 10 kb).
  
  To identify the interactions mediated by CTCF, promoters and enhancers, the genome has been partitioned into element regions (CTCF, promoters, enhancers), element nearby regions and distal regions (remaining regions) at equal bin size (e.g. 10 kb). The read counts between distal regions and elements regions were  used as backgrounds. The distribution of background was evaluated and observed to be well consistent with negative binomial distribution. The background model is fitted for each genomic distance. When the distance between chromatin elements exceeded 1 Mb, all read counts between elements with the distance more than 1 Mb were used for fitting. The negative binomial distribution parameters were evaluated for each distance, then a P value was assigned to each interaction between two elements. Finally, the interactions from all distances were combined and Benjamini-Hochberg (BH) method was used for multiple testing correction.
