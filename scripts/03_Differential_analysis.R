library(DESeq2)

countsData_total = read.table("./data_analyses/01_Expression_level_screening/all_counts_1.0_expression_filter.txt", header=T, row.names=1)


info_file_name = './data_analyses/02_Dimension_Reduction_and_Clustering/Grouping_of_samples.csv'
info_total = read.csv(file = info_file_name)

## KO vs WT
dir.create('./data_analyses/03_Differential_analysis/KO_vs_WT')

for (tissue in unique(info_total[['group_1']])) { # tissue = 'RA'
    
  info = info_total[info_total[["group_1"]] == tissue, ]
  countsData = countsData_total[ , info[['sampleID']]]

  countsData = round(countsData)
  
  condition = as.factor(info[["group_3"]])

  infoData = data.frame(row.names = info[["sampleID"]], condition)

  dds = DESeqDataSetFromMatrix(countsData, infoData, design= ~ condition)

  dds = DESeq(dds)
  

  sample_combine = combn(levels(condition), 2, simplify=F)
  

  for (i in sample_combine) {
    message("start_", i[1], "_vs_", i[2])

    result01 = results(dds, contrast = c("condition", i[1], i[2]))

    result02 = result01[order(result01$pvalue),]

    filename1 = paste("./data_analyses/03_Differential_analysis/KO_vs_WT/", tissue, "_diff_", i[1], "_vs_", i[2], ".csv", sep = "")
    write.csv(result02, file = filename1)
    

    diff_gene_deseq2 = subset(result02, padj < 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1))
    

    filename2 = paste("./data_analyses/03_Differential_analysis/KO_vs_WT/", tissue, "_signification_padj0.05_logFold1_diff_", i[1], "_vs_", i[2], ".csv", sep = "")
    write.csv(diff_gene_deseq2, file = filename2)
    

    library(ggplot2)
    library(ggrepel)

    result02$threshold = factor(ifelse(result02$padj < 0.05 & abs(result02$log2FoldChange) >= 1, ifelse(result02$log2FoldChange >= 1, 'Up', 'Down'), 'NS'), levels=c('Up', 'Down', 'NS'))

    filename3 = paste("./data_analyses/03_Differential_analysis/KO_vs_WT/", tissue, "_signification_padj0.05_logFold1_diff_", i[1], "_vs_", i[2], "_marked.csv", sep = "")
    write.csv(result02, file = filename3)
    
    
    plot_data = na.omit(as.data.frame(result02))
    
    plot_data[['padj_modify']][plot_data[['log2FoldChange']] >= 0] = -log10(plot_data[['padj']][plot_data[['log2FoldChange']] >= 0])
    plot_data[['padj_modify']][plot_data[['log2FoldChange']] < 0] = log10(plot_data[['padj']][plot_data[['log2FoldChange']] < 0])
    

    diff_gene_VolcanoPlot = ggplot(plot_data, aes(x = log2FoldChange, y = padj_modify, color = threshold, alpha(0.5), ))+
      geom_point(size = 1.5)+
      scale_color_manual(values = c("Up" = "#EFC0B9", "Down" = "#385487", "NS" = "#d7d7d8"))+
      theme_test()+
      theme(legend.title = element_blank())+
      ylab('-log10(p-adj)')+
      xlab('log2(FoldChange)')+
      geom_vline(xintercept=c(-1, 1), lty=2, col="black", lwd=0.25) +
      geom_hline(yintercept = c(-log10(0.05), log10(0.05)), lty=2, col="black", lwd=0.25)+

      scale_y_continuous(labels = function(y) ifelse(y < 0, -y, y)) +
      theme_test()

    VolcanoPlot_name = paste("./data_analyses/03_Differential_analysis/KO_vs_WT/", tissue, "_diff_", i[1], "_vs_", i[2], "_padj0.05_logFD1_VolcanoPlot", ".pdf", sep = "")
    ggsave(diff_gene_VolcanoPlot, filename = VolcanoPlot_name, width = 4, height = 2.9)
    message("end_", i[1], "_vs_", i[2])
  }
}



rm(list=ls())

quit()