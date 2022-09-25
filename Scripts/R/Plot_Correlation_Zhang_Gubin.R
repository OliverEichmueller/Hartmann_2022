# Re Plot Corrplot
library(corrplot)
mod_zhang <- read.csv('Zhang_2020_processed/modules_Zhang_Myeloid.csv', header = T)
zhang <- openxlsx::read.xlsx('Zhang_2020_processed/Modules_Cor_Zhang.xlsx', sheet = 4, rowNames = T)
zhang_aggr <- openxlsx::read.xlsx('Zhang_2020_processed/Modules_Cor_Zhang.xlsx', sheet = 2, rowNames = F)
gubin <- openxlsx::read.xlsx('Macro_DC_2/Comparison_Cell_Paper/Modules_Cor_Gubin.xlsx', sheet = 3, rowNames = T)


col2 = colorRampPalette(c('#67001F', '#B2182B', '#D6604D', '#F4A582',
                          '#FDDBC7', '#FFFFFF', '#D1E5F0', '#92C5DE',
                          '#4393C3', '#2166AC', '#053061'))

for (i in 1:ncol(zhang)) {
  zhang[,i] <- unlist(as.numeric(zhang[,i]))
}

pdf('Zhang_2020_processed/Correlation_Zhang_Myeloid.pdf', width = 10, height = 10)
corrplot(as.matrix(zhang[c("Monocytes", "NC-Monocytes", "TAM-3", "TAM-4", "TAM-5",
                           "TAM-6", "TAM-7","TAM-1", "TAM-2", "TRM"), 
                         c(22,23,21,28,24,27,26,25)]), col = rev(col2(100)))
dev.off()

dc <- readRDS('processed.DS_2/DCs_all_monocle_211117.Rds')

cell_group_df <- dc@colData[,c("barcode", "DC_names")]

aggr_mod_dc <- aggregate_gene_expression(dc, gene_group_df = mod_zhang %>% select(-1), 
                                         cell_group_df = cell_group_df)
cor_dc <- cor(cbind(zhang_aggr, aggr_mod_dc))

pdf('Zhang_2020_processed/Correlation_Zhang_DCs.pdf', width = 5, height = 5)
corrplot::corrplot(as.matrix(cor_dc[16:19,3:7]), 
                   col = rev(col2(100)))
dev.off()


pdf('Macro_DC_2/Comparison_Cell_Paper/correlation_Gubin_name_new.pdf', width = 10, height = 10)
corrplot::corrplot(as.matrix(gubin[c("Monocytes", "NC-Monocytes", "TAM-3", "TAM-4", "TAM-5",
                               "TAM-6", "TAM-7","TAM-1", "TAM-2", "TRM"),
                             c("1:Macro", "4:Macro","5:Macro","3:Macro",   "2:Macro" )]), col = rev(col2(100)))
dev.off()

write.xlsx(as.data.frame(cor_dc[16:19,3:7]), 'Zhang_2020_processed/Cor_DCs.xlsx', rowNames = T, overwrite = T)
