library(ggplot2)
library(tidyr)
library(ggpubr)

# define summarySE function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop, .fun = function(xx, col) {
    c(N = length2(xx[[col]], na.rm=na.rm),
      mean = mean(xx[[col]], na.rm=na.rm),
      sd = sd(xx[[col]], na.rm=na.rm)
    )
  },
  measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


ratio_file = './data_analyses/10_Inferring_muscle_fiber_ratios/Evaluation_results_of_sample_cell_proportion.csv'
cell_ratio = read.csv(file = ratio_file, header = TRUE, row.names = 1)


group_file = './data_analyses/02_Dimension_Reduction_and_Clustering/Grouping_of_samples.csv'
group_info = read.csv(file = group_file, header = TRUE, row.names = 1)


cell_ratio = cell_ratio[rownames(group_info), ]
merge_table = cbind(cell_ratio, group_info)
merge_table = merge_table[ , c("I", "IIA", "IIB", "group_1", "group_3")]

merge_table[['sample']] = rownames(merge_table)


plot_table = gather(merge_table, key = 'cell_type', value = 'ratio', 1:3)


data_summary = summarySE(plot_table,
                         measurevar="ratio",
                         groupvars=c("group_1", "group_3", "cell_type"))



P = ggplot(data_summary, aes(x=group_1, y=ratio, fill=group_3))+
  facet_grid(cell_type~.)+
  geom_point(aes(x=group_1, y=ratio, color=group_3), pch=19, position=position_dodge(0.85), size=1.5)+
  geom_bar(stat="identity", position=position_dodge(0.85), alpha = 0.7, width = 0.7)+
  geom_errorbar(aes(ymin = ratio-se, ymax = ratio+se),
                width=0.1,
                position=position_dodge(0.85), 
                color="black",
                alpha = 0.5,
                size=0.5)+
  stat_compare_means(data=plot_table, aes(x=group_1, y=ratio, group = group_3),
                     method = 'wilcox.test',
                     paired = FALSE,
                     label="p.format",
                     show.legend = F,
                     hjust=0.5,
                     vjust=1,
                     angle=20,
                     position = "identity")+
  scale_fill_manual(values=c(KO = '#CD3333', WT = '#1874CD'))+
  scale_color_manual(values=c(KO = '#CD3333', WT = '#1874CD'))+
  theme(axis.text.x=element_text(angle=40, hjust=1, colour="black", size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size = 8),
        panel.border = element_blank(),axis.line = element_line(colour = "black", linewidth=1),
        legend.text=element_text(colour="black",
                                 size=8),
        legend.title=element_text(colour="black",
                                  size=9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_test()+
  ylab("Proportion")+xlab("")

save_file_name2 = paste0('./data_analyses/10_Inferring_muscle_fiber_ratios/Comparison_of_muscle_fiber_ratios.pdf')
ggsave(P, filename = save_file_name2, width = 6, height = 6)


rm(list = ls())
quit()