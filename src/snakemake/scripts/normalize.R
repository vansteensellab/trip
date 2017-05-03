#!/usr/bin/R
args = commandArgs(trailingOnly=T)
cDNA_input = read.table(args[1], col.names=c('barcode','count'), stringsAsFactors=F)
gDNA_input = read.table(args[2], col.names=c('barcode','count'), stringsAsFactors=F)

if (length(args) == 4){
  spike_input = read.table(args[3], col.names=c('barcode','count'), stringsAsFactors=F)
  spike_sum = sum(spike_input$count)
  output = args[4]
} else if (length(args) == 3){
  output = args[3]
}

match_vec = match(cDNA_input$barcode, gDNA_input$barcode)
count_table = cbind.data.frame(cDNA_input$count, gDNA_input$count[match_vec])
cpm_table = t(t(count_table) / colSums(count_table) * 1000000)
normalized_by_gDNA = cpm_table[,1] / cpm_table[,2]
output_data = cbind.data.frame(cDNA_input$barcode, count_table, cpm_table, normalized_by_gDNA)
colnames(output_data) = c('barcode', 'cDNA_count', 'gDNA_count', 'cDNA_cpm', 'gDNA_cpm', 'normalized_by_gDNA')
if (length(args)==4){
  output_data$normalized_by_spike = output_data$normalized_by_gDNA / spike_sum * 1000000
}

write.table(output_data, file=output, row.names=F, quote=F, sep='\t')
