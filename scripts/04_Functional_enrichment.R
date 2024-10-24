library("biomaRt")
library("clusterProfiler")
library("R.utils")
library("org.Hs.eg.db")
library("ggplot2")


ensembl = useMart("ensembl")

hsa = useMart('ensembl',dataset = "hsapiens_gene_ensembl", host = "https://oct2022.archive.ensembl.org/")

gga = useMart('ensembl',dataset = "ggallus_gene_ensembl", host = "https://oct2022.archive.ensembl.org/")

R.utils::setOption("clusterProfiler.download.method",'auto')
createKEGGdb::create_kegg_db("hsa")

install.packages("./KEGG.db_1.0.tar.gz", type="source", force=TRUE)

library(KEGG.db)

for (tissue in c('liver','stroma','swf','syf')) { # tissue = 'liver'
  

  diff_file = paste0("./data_analyses/03_Differential_analysis/", tissue, '_signification_padj0.05_logFold1_diff_high_production_vs_low_production_marked.csv')
  diff_table = na.omit(read.csv(diff_file, header = TRUE, row.names = 1))
  

  diff_table = diff_table[!diff_table[['threshold']] == 'Non-Sig', ]
  
  chicken_ensembl_id = rownames(diff_table)

  tissue_diff_gene = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id",
                                        values = chicken_ensembl_id, mart = gga,
                                        attributesL = c("external_gene_name", "entrezgene_id", "chromosome_name", 'ensembl_gene_id'),
                                        martL = gga, uniqueRows = T)
  colnames(tissue_diff_gene) = c('ensembl_gene_id', 'gene_name', 'entrezgene_id', 'chromosome', 'ensembl_gene_id')
  
  write.csv(tissue_diff_gene, file = paste0("./data_analyses/04_Functional_enrichment/", tissue, "_diff_gene_info.csv"), row.names = F)
  

  gene_name = tissue_diff_gene$gene_name
  

  gene_name = gene_name[!gene_name == '']
  

  go_all = enrichGO(gene       = gene_name,
                    OrgDb      = org.Hs.eg.db,
                    keyType    = "SYMBOL",
                    ont        = "ALL",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2,
                    minGSSize = 3,
                    maxGSSize = 500,
                    readable = TRUE)

  write.table(go_all, file = paste0("./data_analyses/04_Functional_enrichment/", tissue, "_GO_all_terms_", ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")

  length_sig_go = nrow(go_all@result[go_all@result[["p.adjust"]] < 0.05, ])
  if (length_sig_go >= 1){

    GO_dotplot = dotplot(go_all,
                         x = "GeneRatio",
                         color = "p.adjust",
                         font.size = 10,
                         label_format = 300,
                         
                         showCategory=20,
                         
                         split="ONTOLOGY")+

      facet_grid(ONTOLOGY~., scale='free')+

      scale_size_continuous(range=c(0.5, 4.5))+

      scale_colour_gradientn(values = seq(0, 1, 0.2),
                             colors = rev(c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25")))

    if (tissue == 'stroma') {
      ggsave(GO_dotplot, file = paste0("./data_analyses/04_Functional_enrichment/", tissue, "_GO_all_terms_", ".pdf"), width = 15, height = 10)
    } else {
      ggsave(GO_dotplot, file = paste0("./data_analyses/04_Functional_enrichment/", tissue, "_GO_all_terms_", ".pdf"), width = 10, height = 10)
    }
  }
  
  
  tissue_diff_gene_human_ncbi = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name",
                                       values = gene_name,
                                       mart = hsa,
                                       attributesL = c("external_gene_name", "entrezgene_id"),
                                       martL = hsa,
                                       uniqueRows = T)
  colnames(tissue_diff_gene_human_ncbi) = c('gene_name', 'gene_name', 'entrezgene_id')
  
  kegg = enrichKEGG(gene = unique(na.omit(tissue_diff_gene_human_ncbi[["entrezgene_id"]])),
                    organism = "hsa",
                    keyType = "kegg",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2,
                    minGSSize = 3,
                    maxGSSize = 500,
                    use_internal_data = TRUE)
  write.table(kegg, file = paste0("./data_analyses/04_Functional_enrichment/", tissue, "_KEGG_pathways_", ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  
  length_sig_kegg = nrow(kegg@result[kegg@result[["p.adjust"]] < 0.05, ])
  if (length_sig_kegg >= 1){
    KEGG_dotplot = dotplot(kegg, x = "GeneRatio",
                           color = "p.adjust",
                           font.size = 10,
                           label_format = 300,
                           showCategory=50,
    )+
      scale_size_continuous(range=c(0.5, 4.5))+
      scale_colour_gradientn(values = seq(0, 1, 0.2),
                             colors = rev(c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25")))
    ggsave(KEGG_dotplot, file = paste0("./data_analyses/04_Functional_enrichment/", tissue, "_KEGG_all_terms_", ".pdf"), width = 9, height = 2+0.3*length_sig_kegg)
  }
}


rm(list=ls())
quit()