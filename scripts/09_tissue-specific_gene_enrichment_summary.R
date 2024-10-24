# load R packages
library(readxl)
library(ggplot2)
library(dplyr)

dir.create('./data_analyses/09_tissue-specific_gene_function_enrichment/')
path_file = "./data_analyses/09_tissue-specific_gene_function_enrichment/enrichment_results_summary.xlsx"

table_1 = read_excel(path = path_file, sheet = 1, range = NULL, col_names = TRUE)
table_2 = read_excel(path = path_file, sheet = 2, range = NULL, col_names = TRUE)
table_3 = read_excel(path = path_file, sheet = 3, range = NULL, col_names = TRUE)


table_3$Description = factor(table_3$Description, levels =table_2$Description, ordered = TRUE)
table_3$Group = factor(table_3$Group, levels = rev(table_1$Group), ordered = TRUE)

P = ggplot(table_3, aes(x = Description, y = Group, size = Number, color = -LogP))+
  geom_point(stat ="identity")+
  scale_colour_gradientn(colors = c("#183f7f", "#4281b6", "#faefef", "#b1257a", "#572266"))+
  labs(
    color=expression(-LogP),
    size="Count Number",
    x=""
  )+
  scale_size_continuous(range = c(3, 10))+
  theme_test()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    )

# save as pdf
save_name = paste0("./data_analyses/09_tissue-specific_gene_function_enrichment/", "enrichment_results_summary.pdf")
ggsave(P, filename = save_name, width = 6.5, height = 5)

rm(list=ls())
quit()