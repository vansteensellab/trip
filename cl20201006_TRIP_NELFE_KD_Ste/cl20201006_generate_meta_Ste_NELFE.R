library(stringr)

file_list = list.files('/DATA/usr/c.leemans/data/forge/s.manzo/6155/fastq_files',
                       full.names=T)


file_df = str_match(file_list, '.*/6155_[0-9]+_(.*?)_(.*_I+).*')
colnames(file_df) = c('file', 'PCR_type', 'ID')

write.table(file_df, "cl20201006_TRIP_NELFE_KD_meta.txt",
            sep='\t', quote=F, row.names=F)
