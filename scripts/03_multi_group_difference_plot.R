
library(ggplot2)
library(dplyr)
library(ggrepel)

dir.create('./data_analyses/03_Differential_analysis/')
# manual order
manual_order = c('TIP', 'BIP', 'RA', 'LD', 'PMM', 'SOL', 'GAS', 'EDL', 'TA')

diff_total = read.csv("./data_analyses/03_Differential_analysis/diff_total_table.csv", header=T)
# factor
diff_total[['group']] = factor(diff_total[["group"]], levels = manual_order)
diff_total[['threshold']] = factor(diff_total[['threshold']], levels = c('Up','Down'))


back_df = diff_total %>% group_by(group) %>% summarise(min = min(log2FoldChange),
                                                       max = max(log2FoldChange),
                                                       down_n = sum(log2FoldChange < 0),
                                                       up_n = sum(log2FoldChange > 0))
# factor
back_df[['group']] = factor(back_df[["group"]], levels = manual_order)


# label df
label_df = data.frame(group = back_df[['group']],
                      y = 0,
                      label = back_df[['group']])
# factor
label_df[['group']] = factor(label_df[["group"]], levels = manual_order)


colors = c("TIP" = '#466C1A',
           "BIP" = '#166C1A',
           "RA" = '#CD5C5C',
           "LD" = '#FF6A6A',
           "PMM" = '#CD3333',
           "SOL" = '#8EE5EE',
           "GAS" = '#63B8FF',
           "EDL" = '#2DA2D6',
           "TA" = '#1874CD')


logFC = 1
# plot
P = ggplot()+
  geom_col(data = back_df,
           mapping = aes(x = group, y = max),
           fill = "#dcdcdc",
           alpha = 0.55,
           show.legend = FALSE)+
  geom_col(data = back_df,
           mapping = aes(x = group, y = min),
           fill = "#dcdcdc",
           alpha = 0.55,
           show.legend = FALSE)+
  geom_jitter(data = diff_total,
              aes(x = group, y = log2FoldChange, color = threshold),
              size = 0.75,
              width =0.4)+
  geom_tile(data = label_df,
            aes(x = group, y = y, fill = group),
            height=1.5*logFC,
            color = "black",
            #fill = mycol,
            alpha = 1.0,
            show.legend = FALSE)+
  geom_text(data = label_df,
            aes(x = group, y = y, label = label),
            size = 3.5,
            color = "white")+
  geom_text(data = back_df,
            aes(x = group, y = max(max+1), label = up_n),
            size = 2,
            color = "black")+
  geom_text(data = back_df,
            aes(x = group, y = min(min-1), label = down_n),
            size = 2,
            color = "black")+
  scale_color_manual(name=NULL, values = c(Up = "#F39B7F", Down = "#4cbbd5"))+
  # Adjust the size of legend point
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  scale_fill_manual(name=NULL, values=colors)+
  labs(x = "", y = "log2(FoldChange)")+
  #theme_minimal()+
  theme_test()+
  theme(
    axis.title = element_text(size = 10, color = "black"),
    #axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.line.x = element_blank(),
    #axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1, 0),
    legend.text = element_text(color = "black", size = 10)
  )

# save name
save_name = paste0('./data_analyses/03_Differential_analysis/', 'diff_summary.pdf')
ggsave(P, filename = save_name, width = 5.5, height = 4.5)



rm(list=ls())

quit()