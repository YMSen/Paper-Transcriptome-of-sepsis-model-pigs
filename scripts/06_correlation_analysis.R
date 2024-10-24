library(pheatmap)
library(ggplot2)

dir.create('./data_analyses/06_correlation_analysis/')

tpmData_total = read.table("./data_analyses/01_Expression_level_screening/all_quant_TPM_1.0_expression_filter.txt", header=T, row.names=1)

info_file_name = './data_analyses/02_Dimension_Reduction_and_Clustering/Grouping_of_samples.csv'
info_total = read.csv(file = info_file_name)

# filter 'CLP' data
info_total = info_total[info_total[["group_2"]] == "CLP", ]
tpmData_total = tpmData_total[ , info_total[['sampleID']]]


for (type in unique(info_total[["group_3"]])) { # type = "WT"
  
  info = info_total[info_total[["group_3"]] == type, ]
  tpmData = tpmData_total[ , info[['sampleID']]]
  
  # top 5000 genes with variation
  genes = names(tail(sort(apply(tpmData,
                                MARGIN = 1,
                                FUN = function(x) var(x))), n = 3500))
  tpmData = tpmData[genes, ]
  
  spearman_cor_matrix = cor(tpmData, method = 'spearman')
  min = min(spearman_cor_matrix)
  max = max(spearman_cor_matrix)
  
  write.csv(spearman_cor_matrix, file = paste0("./data_analyses/06_correlation_analysis/", type, "_spearman_cor_matrix.csv"), row.names = T)
  
  
  annotation_col = data.frame(group = c(rep("BIP", 3),
                                        rep("EDL", 3),
                                        rep("GAS", 3),
                                        rep("LD", 3),
                                        rep("PMM", 3),
                                        rep("RA", 3),
                                        rep("SOL", 3),
                                        rep("TA", 3),
                                        rep("TIP", 3)))
  manual_order = c('TIP', 'BIP', 'RA', 'LD', 'PMM', 'SOL', 'GAS', 'EDL', 'TA')
  annotation_col[['group']] = factor(annotation_col[['group']], levels = manual_order)
  
  
  rownames(annotation_col) = rownames(spearman_cor_matrix)
  ann_colors = list(group = c("TIP" = '#466C1A',
                              "BIP" = '#166C1A',
                              "RA" = '#CD5C5C',
                              "LD" = '#FF6A6A',
                              "PMM" = '#CD3333',
                              "SOL" = '#8EE5EE',
                              "GAS" = '#63B8FF',
                              "EDL" = '#2DA2D6',
                              "TA" = '#1874CD'
                              ))
  
  save_pdf = paste0("./data_analyses/06_correlation_analysis/", type, "_multi_species_tpm_spearman_cor_heatmap.pdf")

  pheatmap(spearman_cor_matrix,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           cellwidth = 12,
           cellheight = 12,
           show_rownames = TRUE,
           show_colnames = FALSE,
           scale = "none",
           border_color = 'NA',
           breaks=seq(min, max, (max-min)/50),
           color = colorRampPalette(c("#414887", "#9692B9", "#DFDDE4", "#F6EAE0","#EFC18A","#E29C35"))(50),
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           fontsize = 9,
           fontsize_row = 9,
           fontsize_col = 9,
           angle_col = 45,
           main = "",
           filename = save_pdf)
}



rm(list=ls())
quit()