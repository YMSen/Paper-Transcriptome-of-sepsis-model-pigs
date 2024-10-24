library(ggplot2)
library(dplyr)

dir.create('./data_analyses/05_mapping_info/')
# read mapping info
mappingData = read.table("./data_analyses/05_mapping_info/map_info.xls", header = TRUE)
rownames(mappingData) = mappingData[['samples']]
mappingData[['sampleID']] = mappingData[['samples']]

# read info data
info_file_name = './data_analyses/02_Dimension_Reduction_and_Clustering/Grouping_of_samples.csv'
info_total = read.csv(file = info_file_name)

# filter 'CLP' data
info_total = info_total[info_total[["group_2"]] == "CLP", ]
mappingData = mappingData[info_total[['sampleID']], ]

# add group info
merge_table = merge(x = mappingData,
                    y = info_total,
                    by.x = 'sampleID',
                    by.y = 'sampleID')

merge_table[['group_1']] = factor(merge_table[['group_1']], levels = c("TIP","BIP","RA","LD","PMM",
                                                                       "SOL","GAS","EDL","TA"))

# sort by group
merge_table = merge_table %>%
  arrange(group_1)
merge_table = merge_table %>%
  arrange(group_1, totalMapped...) %>%
  group_by (group_1) %>%
  mutate (rank = rank(-totalMapped...))
merge_table = merge_table %>%
  arrange(group_1, rank) %>%
  group_by (group_1) %>%
  mutate (rank = rank(rank))

merge_table[['sampleID']] = factor(merge_table[['sampleID']], levels = merge_table[['sampleID']])

# generate a series of gradient colors based on the main colors
# main colors
main_colors = c('#466C1A','#166C1A','#CD5C5C','#FF6A6A','#CD3333',
                '#8EE5EE','#63B8FF','#2DA2D6','#1874CD')
# generate six secondary colors
generate_secondary_colors = function(color) {
  
  library(colorspace)
  
  sapply(seq(0.1, 0.6, by = 0.1), function(amount) {
    lighten(color, amount = amount)
  })
}

# secondary colors
secondary_colors_list = lapply(main_colors, generate_secondary_colors)
secondary_colors = unlist(secondary_colors_list)

figure='color'

if (figure == 'grey'){
  P = ggplot(merge_table, aes(x=sampleID, y=totalMapped...))+
    geom_bar(fill = 'grey45', color = 'grey75', stat="identity", position=position_dodge(0.8), alpha = 0.7, width = 1.0)+
    geom_point(aes(x=sampleID, y=totalMapped..., color=group_1), pch=19, position=position_dodge(0.8), alpha = 0.7, size=1.5)+
    scale_color_manual(values=c("TIP" = '#466C1A',
                                "BIP" = '#166C1A',
                                "RA" = '#CD5C5C',
                                "LD" = '#FF6A6A',
                                "PMM" = '#CD3333',
                                "SOL" = '#8EE5EE',
                                "GAS" = '#63B8FF',
                                "EDL" = '#2DA2D6',
                                "TA" = '#1874CD'))+
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y=element_text(size=10),
          axis.title.y=element_text(size = 10),
          legend.text=element_text(colour="black",
                                   size=10),
          legend.title=element_text(colour="black",
                                    size=11),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
    geom_hline(yintercept = 95, color = "white", linetype = "dashed", size = 0.25)+
    ylab("total mapping ratio")+xlab("")
  
  save_file_name = paste0('./data_analyses/05_mapping_info/', 'mapping_rate_grey_col_plot.pdf')
  ggsave(P, filename = save_file_name, width = 6, height = 3.2)
  
} else if (figure == 'color') {
  P = ggplot(merge_table, aes(x=sampleID, y=totalMapped...))+
    geom_point(aes(x=sampleID, y=totalMapped..., color=group_1), pch=19, position=position_dodge(0.8), alpha = 0.7, size=1.5)+
    geom_bar(aes(fill=sampleID), stat="identity", position=position_dodge(0.8), alpha = 0.7, width = 1.0)+
    scale_fill_manual(values=secondary_colors)+
    scale_color_manual(values=c("TIP" = '#466C1A',
                                "BIP" = '#166C1A',
                                "RA" = '#CD5C5C',
                                "LD" = '#FF6A6A',
                                "PMM" = '#CD3333',
                                "SOL" = '#8EE5EE',
                                "GAS" = '#63B8FF',
                                "EDL" = '#2DA2D6',
                                "TA" = '#1874CD'))+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(size=10),
          axis.title.y=element_text(size = 10),
          legend.text=element_text(colour="black",
                                   size=10),
          legend.title=element_text(colour="black",
                                    size=11),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
    geom_hline(yintercept = 95, color = "grey35", linetype = "dashed", size = 0.25)+
    ylab("total mapping ratio")+xlab("")
  
  save_file_name = paste0('./data_analyses/05_mapping_info/', 'mapping_rate_color_col_plot.pdf')
  ggsave(P, filename = save_file_name, width = 6, height = 3.2)
}


rm(list=ls())
quit()