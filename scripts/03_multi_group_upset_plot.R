library(ggplot2)
library(ggrepel)
library(sets)
library(UpSetR)

dir.create('./data_analyses/03_Differential_analysis/')

# manual order
manual_order = c('TIP', 'BIP', 'RA', 'LD', 'PMM', 'SOL', 'GAS', 'EDL', 'TA')

diff_total = read.csv("./data_analyses/03_Differential_analysis/diff_total_table.csv", header=T)
# factor
diff_total[['group']] = factor(diff_total[["group"]], levels = manual_order)
diff_total[['threshold']] = factor(diff_total[['threshold']], levels = c('Up','Down'))


for (type in c('Up', 'Down')) { # type = 'Up'
  # one type data
  one_type = diff_total[diff_total[['threshold']]==type, ]
  # IDs as list by each group
  result = split(one_type, one_type[['group']])
  result_list = lapply(result, function(x) x[['ID']])
  
  calculate_unique_elements = function(sets_list) {
    unique_elements = lapply(names(sets_list), function(set_name) {
      other_sets = sets_list[names(sets_list) != set_name]
      unique_elements = setdiff(sets_list[[set_name]], Reduce(union, other_sets))
      return(unique_elements)
    })
    names(unique_elements) = names(sets_list)
    return(unique_elements)
  }
  

  unique_elements = calculate_unique_elements(result_list)

  all_intersection = Reduce(intersect, result_list)
  

  stats = list()
  for (name in names(unique_elements)) {
    stats[[name]] = length(unique_elements[[name]])
  }
  stats[[paste(names(unique_elements), collapse = '&')]] = length(all_intersection)
  
  # unlist
  stats = unlist(stats)

  save_pdf = paste0("./data_analyses/03_Differential_analysis/", type, "_gene_upset.pdf")
  pdf(save_pdf, width = 7, height = 5.5)
  
  if (type == 'Up') {
    queries = list(
      list(query = elements, params = list("TIP"), color = '#466C1A', active = T),
      list(query = elements, params = list("BIP"), color = '#166C1A', active = T), 
      list(query = elements, params = list("RA"), color = '#CD5C5C', active = T),
      list(query = elements, params = list("LD"), color = '#FF6A6A', active = T),
      list(query = elements, params = list("PMM"), color = '#CD3333', active = T),
      list(query = elements, params = list("SOL"), color = '#8EE5EE', active = T),
      list(query = elements, params = list("GAS"), color = '#63B8FF', active = T),
      list(query = elements, params = list("EDL"), color = '#2DA2D6', active = T),
      list(query = elements, params = list("TA"), color = '#1874CD', active = T),
      list(query = intersects, params = list('TIP', 'BIP', 'RA', 'LD', 'PMM', 'SOL', 'GAS', 'EDL', 'TA'), color = "grey65", active = T)
    )
  } else if (type == 'Down') {
    queries = list(
      list(query = elements, params = list("TIP"), color = '#466C1A', active = T),
      list(query = elements, params = list("BIP"), color = '#166C1A', active = T), 
      list(query = elements, params = list("RA"), color = '#CD5C5C', active = T),
      list(query = elements, params = list("LD"), color = '#FF6A6A', active = T),
      list(query = elements, params = list("PMM"), color = '#CD3333', active = T),
      list(query = elements, params = list("SOL"), color = '#8EE5EE', active = T),
      list(query = elements, params = list("GAS"), color = '#63B8FF', active = T),
      list(query = elements, params = list("EDL"), color = '#2DA2D6', active = T),
      list(query = elements, params = list("TA"), color = '#1874CD', active = T)
    )
  }
  
  upset(fromExpression(stats),
        # fromList(result_list),
        sets = manual_order,
        keep.order = TRUE,
        nsets = 10,
        nintersects = NA,
        order.by = "freq",
        decreasing = c(TRUE),
        mb.ratio = c(0.6,0.4),
        text.scale = c(1.5, 1.5, 1.5, 1, 1.5, 1.5),
        mainbar.y.label = "Intersect Number",
        sets.x.label = "Sample Number",
        show.numbers = "yes",
        set_size.show = TRUE,
        set_size.numbers_size = 8.5,
        number.angles = 0,
        point.size = 5,
        line.size = 1.5,
        sets.bar.color = c("TIP" = '#466C1A',
                           "BIP" = '#166C1A',
                           "RA" = '#CD5C5C',
                           "LD" = '#FF6A6A',
                           "PMM" = '#CD3333',
                           "SOL" = '#8EE5EE',
                           "GAS" = '#63B8FF',
                           "EDL" = '#2DA2D6',
                           "TA" = '#1874CD'),
        matrix.color = "gray70",
        main.bar.color = "gray70",
        queries = queries
        )
  
  #
  dev.off()
}



rm(list=ls())

quit()