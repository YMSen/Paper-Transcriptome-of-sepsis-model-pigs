library(pheatmap)
library(ggplot2)
library(dplyr)
library(gprofiler2)

dir.create('./data_analyses/08_tissue-specific_gene_analysis/')

# read TAU data
tauData = read.table("./data_analyses/08_tissue-specific_gene_analysis/all_quant_TPM_1.0_expression_filter.TAU.txt", header = TRUE)
rownames(tauData) = tauData[['ID']]
# filter data
tauData = tauData[tauData[["Specific...0.75."]] == TRUE, ]

# perform homologous gene transformation for each gene set
for (part in unique(tauData[['Tissue']])) { # part = 'EDL'
  one_part = tauData[tauData[['Tissue']]==part, ]
  
  orth_info = gorth(
    query = one_part[['ID']],
    source_organism = "sscrofa",
    target_organism = "mmusculus",
    numeric_ns = "",
    mthreshold = Inf,
    filter_na = FALSE
  )
  
  # save information
  save_name = paste0('./data_analyses/08_tissue-specific_gene_analysis/', part, '_orth_genes_info.csv')
  write.csv(orth_info, file = save_name, row.names = FALSE, quote = FALSE)
}


# read gene expression data
gene_table = read.table('./data_analyses/08_tissue-specific_gene_analysis/all_quant_TPM_1.0_expression_filter.txt', header = TRUE)
rownames(gene_table) = gene_table[['ID']]

# add group info
merge_table = merge(x = tauData,
                    y = gene_table,
                    by.x = 'ID',
                    by.y = 'ID')
rownames(merge_table) = merge_table[['ID']]
# factor
merge_table[['Tissue']] = factor(merge_table[['Tissue']], levels = c("TIP","BIP","RA","LD","PMM",
                                                             "SOL","GAS","EDL","TA"))
# sort by Tissue
merge_table = merge_table %>% arrange(Tissue)

# filter data
plot_data = merge_table[ , (ncol(tauData)+1) : ncol(merge_table)]
plot_data = as.data.frame(t(plot_data))
plot_data[['sampleID']] = rownames(plot_data)

# read info data
info_file_name = './data_analyses/02_Dimension_Reduction_and_Clustering/Grouping_of_samples.csv'
info_total = read.csv(file = info_file_name)
# filter 'CLP' data
info_total = info_total[info_total[["group_2"]] == "CLP", ]

# merge
plot_data = merge(x = plot_data,
                  y = info_total,
                  by.x = 'sampleID',
                  by.y = 'sampleID')

# as factor
plot_data[['group_1']] = factor(plot_data[['group_1']], levels = c("TIP","BIP","RA","LD","PMM",
                                                                   "SOL","GAS","EDL","TA"))

# sort by group
plot_data = plot_data %>% arrange(group_1)
rownames(plot_data) = plot_data[['sampleID']]
# filter data
plot_data_2 = plot_data[ , 2:(nrow(merge_table)+1)]

for (gene_id in colnames(plot_data_2)) {
  plot_data_2[[gene_id]] = scale(as.numeric(plot_data_2[[gene_id]]), center=T, scale=T)
}


data = as.data.frame(t(plot_data_2))

annotation_col = data.frame(Group = factor(plot_data[["group_1"]]))
annotation_row = data.frame(Group = factor(merge_table[["Tissue"]]))

rownames(annotation_col) = rownames(plot_data_2)
rownames(annotation_row) = colnames(plot_data_2)

ann_colors = list(Group = c("TIP" = '#466C1A',
                            "BIP" = '#166C1A',
                            "RA" = '#CD5C5C',
                            "LD" = '#FF6A6A',
                            "PMM" = '#CD3333',
                            "SOL" = '#8EE5EE',
                            "GAS" = '#63B8FF',
                            "EDL" = '#2DA2D6',
                            "TA" = '#1874CD'))


pdf("./data_analyses/08_tissue-specific_gene_analysis/z-score_heatmap_of_tissue-specific_genes.pdf", height = 7, width = 3.5)

pheatmap(as.matrix(data),
         breaks=seq(-1,1,0.02),
         color = colorRampPalette(c("#295f32", "#a0c769", "#fbeeb0", "#d86e43", "#942320"))(100),
         cluster_cols = F,
         cluster_rows = F,
         border_color= NA,
         show_rownames = F,
         show_colnames = F,
         main = "",
         display_numbers=F,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors[1],
         angle_col = 90,
         fontsize_row = 12,
         fontsize_col = 12
)+theme_test()

dev.off()


rm(list=ls())
quit()