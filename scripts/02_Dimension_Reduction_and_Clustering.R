
library(Rtsne)
library(ggplot2)
library(ggforce)


data_1 = read.table(file = "./data_analyses/01_Expression_level_screening/all_quant_TPM_1.0_expression_filter.txt", header = T, check.names = F, row.names = 1)

info_file_name = './data_analyses/02_Dimension_Reduction_and_Clustering/Grouping_of_samples.csv'
info_total = read.csv(file = info_file_name)


info = info_total
used_data = data_1[ , info[['sampleID']]]

data = log2(used_data+1)

data_matrix = t(data)

# Set a seed if you want reproducible results
set.seed(123)
# Run TSNE
tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=11, theta=0.52)

tsne_res = as.data.frame(tsne_out$Y)
colnames(tsne_res) = c("tSNE1", "tSNE2")
tsne_res$sample = info$sampleID
tsne_res$condition = info$group
tsne_res$condition_1 = info$group_1
tsne_res$condition_2 = info$group_2
tsne_res$condition_3 = info$group_3

X_min = min(tsne_res[['tSNE1']]) -abs(min(tsne_res[['tSNE1']]))*0.15
X_max = max(tsne_res[['tSNE1']]) +abs(max(tsne_res[['tSNE1']]))*0.15
Y_min = min(tsne_res[['tSNE2']]) -abs(min(tsne_res[['tSNE2']]))*0.15
Y_max = max(tsne_res[['tSNE2']]) +abs(max(tsne_res[['tSNE2']]))*0.15

manual_order = c('TIP', 'BIP', 'RA', 'LD', 'PMM', 'SOL', 'GAS', 'EDL', 'TA')
tsne_res[['condition_1']] = factor(tsne_res[['condition_1']], levels = manual_order)

P = ggplot(tsne_res, aes(tSNE1, tSNE2))+ 
  geom_point(aes(shape = condition_1, colour = condition_1), size = 1.5)+
  geom_mark_ellipse(aes(color = condition_1), expand = unit(0, "mm"))+
  scale_color_manual(values=c("TIP" = '#466C1A', "BIP" = '#166C1A', "RA" = '#CD5C5C',
                              "LD" = '#FF6A6A', "PMM" = '#CD3333', "SOL" = '#8EE5EE',
                              "GAS" = '#63B8FF', "EDL" = '#2DA2D6', "TA" = '#1874CD'))+
  scale_shape_manual(values = c(1, 2, 3, 4, 5, 6, 7, 8, 9))+
  theme_test()+
  coord_cartesian(xlim = c(X_min, X_max), ylim = c(Y_min, Y_max))+
  geom_hline(yintercept = 0, lty=2, col="grey", lwd=0.3)+ 
  geom_vline(xintercept = 0, lty=2, col="grey", lwd=0.3)+
  xlab("t-SNE 1")+
  ylab("t-SNE 2")+
  labs(title = "")+
  guides(fill = guide_legend(title = 'group'),
         color = guide_legend(title = 'group'),
         shape = guide_legend(title = 'group'))

save_file_name_1 = paste0("./data_analyses/02_Dimension_Reduction_and_Clustering/", 'all_sample_t-SNE_plot.pdf')
ggsave(P, filename = save_file_name_1, width = 5, height = 4)


for (group in unique(info_total[['group_1']])) {
  
  dir.create(paste0('./data_analyses/02_Dimension_Reduction_and_Clustering/', group))
  
  info = info_total[info_total[["group_1"]] == group, ]
  used_data = data_1[ , info[['sampleID']]]
  
  data = log2(used_data+1)
  
  data_matrix = t(data)
  
  # Set a seed if you want reproducible results
  set.seed(123)
  # Run TSNE
  if (group == "BIP") {
    tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=0.95, theta=0.95)
  } else if (group == "EDL") {
    tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=0.5, theta=0.5)
  } else if (group == "GAS") {
    tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=0.87, theta=0.815)
  } else if (group == "LD") {
    tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=0.802, theta=0.75)
  } else if (group == "PMM") {
    tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=0.95, theta=0.65)
  } else if (group == "RA") {
    tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=0.90, theta=0.5)
  } else if (group == "SOL") {
    tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=0.985, theta=0.89)
  } else if (group == "TA") {
    tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=0.95, theta=0.70)
  } else if (group == "TIP") {
    tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=0.95, theta=0.75159)
  } 
  
  tsne_res = as.data.frame(tsne_out$Y)
  colnames(tsne_res) = c("tSNE1", "tSNE2")
  tsne_res$sample = info$sampleID
  tsne_res$condition = info$group
  tsne_res$condition_1 = info$group_1
  tsne_res$condition_2 = info$group_2
  tsne_res$condition_3 = info$group_3
  
  X_min = min(tsne_res[['tSNE1']]) -abs(min(tsne_res[['tSNE1']]))*0.30
  X_max = max(tsne_res[['tSNE1']]) +abs(max(tsne_res[['tSNE1']]))*0.30
  Y_min = min(tsne_res[['tSNE2']]) -abs(min(tsne_res[['tSNE2']]))*0.30
  Y_max = max(tsne_res[['tSNE2']]) +abs(max(tsne_res[['tSNE2']]))*0.30
  
  P = ggplot(tsne_res, aes(tSNE1, tSNE2))+ 
    geom_point(aes(shape = condition_3, colour = condition_3), size = 1.5)+
    geom_mark_ellipse(aes(color = condition_3), expand = unit(0, "mm"))+
    scale_color_manual(values=c(KO = '#CD3333', WT = '#1874CD'))+
    scale_fill_manual(values=c(KO = '#CD3333', WT = '#1874CD'))+
    scale_shape_manual(values = c(1, 2))+
    theme_test()+
    coord_cartesian(xlim = c(X_min, X_max), ylim = c(Y_min, Y_max))+
    xlab("t-SNE 1")+
    ylab("t-SNE 2")+
    labs(title = "")+
    guides(fill = guide_legend(title = 'group'),
           color = guide_legend(title = 'group'),
           shape = guide_legend(title = 'group'))
  

  save_file_name = paste0("./data_analyses/02_Dimension_Reduction_and_Clustering/", group, "/", group, '_t-SNE_plot.pdf')
  ggsave(P, filename = save_file_name, width = 5, height = 4)
}


rm(list=ls())

