TTRIP K562 KRAB D2
==================

1. Read and process raw data
----------------------------

```{r}
require(stringr)
require(ggplot2)

counts = read.table('scratch/trip/norm_exp_cl20160702/bc_count.txt', stringsAsFactors=F, header=T, row.names=1)
spike_counts = read.table('scratch/trip/spikein_norm_exp_cl20160702/bc_count.txt', stringsAsFactors=F, header=T, row.names=1)
bc_repeats = read.table('scratch/trip/bc_repeats.txt', stringsAsFactors=F)
colnames(bc_repeats) = c('count', 'barcode', 'family', 'name')

bc_unique = lapply(unique(bc_repeats[,2]),
                   function(x,bc_repeats){
                      this_repeats = bc_repeats[bc_repeats[,2]==x,,drop=F]
                      return(this_repeats[which(this_repeats[,1]==max(this_repeats[,1],na.rm=T))[1],])
                   },bc_repeats)
bc_unique = do.call(rbind, bc_unique)

###IMPORTANT:
#for samples apart from KRAB the usual order of exp1 - POI-GaL4, exp2 - Gal4, exp3 - POI is changed
#now it is exp1 - Gal4, exp2 - Gal4POI, exp3 - POI
info = data.frame(POI=as.factor(rep(c('G9a','CBX5','KRAB'), each=12, 2)),
                  condition=as.factor(c(rep(rep(c('GAL4-POI', 'GAL4', 'POI'),4),2),rep(c('GAL4', 'GAL4-POI', 'POI'),4))),
                  type=as.factor(rep(c('norm', 'exp'),each=36)),
                  day=as.factor(rep(c(2,9), each=6, 6)))

krab_names = paste('KRAB', c('GAL4', 'GAL4-POI', 'POI'), c(rep('norm', 12), rep('exp', 12)), c(rep('D2',6),rep('D9',6)), c(rep('1',3),rep('2',3)), sep='_')
colnames(counts)[grep('KRAB', colnames(counts))] = krab_names
colnames(spike_counts)[grep('KRAB', colnames(counts))] = krab_names

cbx_names = paste('CBX5', c('GAL4-POI', 'GAL4', 'POI'), c(rep('norm', 12), rep('exp', 12)), c(rep('D2',6),rep('D9',6)), c(rep('1',3),rep('2',3)), sep='_')
colnames(counts)[grep('CBX5', colnames(counts))] = cbx_names
colnames(spike_counts)[grep('CBX5', colnames(counts))] = cbx_names

g9a_names = paste('G9a', c('GAL4-POI', 'GAL4', 'POI'), c(rep('norm', 12), rep('exp', 12)), c(rep('D2',6),rep('D9',6)), c(rep('1',3),rep('2',3)), sep='_')
colnames(counts)[grep('G9a', colnames(counts))] = g9a_names
colnames(spike_counts)[grep('G9a', colnames(counts))] = g9a_names

sum_spike_counts = colSums(spike_counts)
exp_norm_counts = t(t(counts[,info$type=='exp'])/sum_spike_counts[info$type=='exp'])

KRAB_norm = counts[,info$POI=='KRAB'&info$type=='norm']
KRAB_above_cut = rowSums(KRAB_norm>=1)==ncol(KRAB_norm)
KRAB_norm = KRAB_norm[KRAB_above_cut,]

G9a_norm = counts[,info$POI=='G9a'&info$type=='norm']
G9a_above_cut = rowSums(G9a_norm>=1)==ncol(G9a_norm)
G9a_norm = G9a_norm[G9a_above_cut,]

KRAB_norm_sum = colSums(KRAB_norm) / 1000000
KRAB_cpm = t(t(KRAB_norm)/KRAB_norm_sum)
KRAB_exp_gDNA = exp_norm_counts[KRAB_above_cut,info$POI[info$type=='exp']=='KRAB']/KRAB_cpm
KRAB_exp_info = info[info$POI=='KRAB' & info$type=='exp',]

G9a_norm_sum = colSums(G9a_norm) / 1000000
G9a_cpm = t(t(G9a_norm)/G9a_norm_sum)
G9a_exp_gDNA = exp_norm_counts[G9a_above_cut,info$POI[info$type=='exp']=='G9a']/G9a_cpm
G9a_exp_info = info[info$POI=='G9a' & info$type=='exp',]


foldch <- function(cond, contr) {
  if ((as.numeric(cond) > 0) & (as.numeric(contr) > 0))
    calc <- as.numeric(cond) / as.numeric(contr)
  if ((as.numeric(cond) == 0) | (as.numeric(contr) == 0))
    calc <- NA
  return(calc)
}

KRAB_GPvsG = apply(cbind(which(KRAB_exp_info$condition=='GAL4-POI'),
                         which(KRAB_exp_info$condition=='GAL4')), 1,
                  function(x){apply(KRAB_exp_gDNA[,x],1,function(x){foldch(x[1],x[2])})})
colnames(KRAB_GPvsG) = paste('GPvsG_day',KRAB_exp_info[KRAB_exp_info$condition=='GAL4-POI','day'],c('.1','.2'),sep='')
KRAB_GPvsP = apply(cbind(which(KRAB_exp_info$condition=='GAL4-POI'),
                         which(KRAB_exp_info$condition=='POI')), 1,
                   function(x){apply(KRAB_exp_gDNA[,x],1,function(x){foldch(x[1],x[2])})})
colnames(KRAB_GPvsP) = paste('GPvsP_day',KRAB_exp_info[KRAB_exp_info$condition=='GAL4-POI','day'],c('.1','.2'),sep='')
KRAB_PvsG = apply(cbind(which(KRAB_exp_info$condition=='POI'),
                        which(KRAB_exp_info$condition=='GAL4')), 1,
                  function(x){apply(KRAB_exp_gDNA[,x],1,function(x){foldch(x[1],x[2])})})
colnames(KRAB_PvsG) = paste('PvsG_day',KRAB_exp_info[KRAB_exp_info$condition=='GAL4-POI','day'],c('.1','.2'),sep='')

KRAB_fc = cbind(KRAB_GPvsG, KRAB_GPvsP, KRAB_PvsG)

G9a_GPvsG = apply(cbind(which(G9a_exp_info$condition=='GAL4-POI'),
                         which(G9a_exp_info$condition=='GAL4')), 1,
                   function(x){apply(G9a_exp_gDNA[,x],1,function(x){foldch(x[1],x[2])})})
colnames(G9a_GPvsG) = paste('GPvsG_day',G9a_exp_info[G9a_exp_info$condition=='GAL4-POI','day'],c('.1','.2'),sep='')
G9a_GPvsP = apply(cbind(which(G9a_exp_info$condition=='GAL4-POI'),
                         which(G9a_exp_info$condition=='POI')), 1,
                   function(x){apply(G9a_exp_gDNA[,x],1,function(x){foldch(x[1],x[2])})})
colnames(G9a_GPvsP) = paste('GPvsP_day',G9a_exp_info[G9a_exp_info$condition=='GAL4-POI','day'],c('.1','.2'),sep='')
G9a_PvsG = apply(cbind(which(G9a_exp_info$condition=='POI'),
                        which(G9a_exp_info$condition=='GAL4')), 1,
                  function(x){apply(G9a_exp_gDNA[,x],1,function(x){foldch(x[1],x[2])})})
colnames(G9a_PvsG) = paste('PvsG_day',G9a_exp_info[G9a_exp_info$condition=='GAL4-POI','day'],c('.1','.2'),sep='')

G9a_fc = cbind(G9a_GPvsG, G9a_GPvsP, G9a_PvsG)

KRAB_fc_mean = sapply(seq(1,ncol(KRAB_fc),by=2),function(x,y){rowMeans(cbind(y[,x],y[,x+1]))},KRAB_fc)
colnames(KRAB_fc_mean) = unlist(str_replace(colnames(KRAB_fc)[seq(1,ncol(KRAB_fc),by=2)],'.1',''))

G9a_fc_mean = sapply(seq(1,ncol(G9a_fc),by=2),function(x,y){rowMeans(cbind(y[,x],y[,x+1]))},G9a_fc)
colnames(G9a_fc_mean) = unlist(str_replace(colnames(G9a_fc)[seq(1,ncol(G9a_fc),by=2)],'.1',''))

KRAB_match_vec=match(bc_unique[,2],rownames(KRAB_fc_mean))
bc_table_KRAB = cbind(bc_unique[!is.na(KRAB_match_vec),],KRAB_fc_mean[KRAB_match_vec[!is.na(KRAB_match_vec)],])

G9a_match_vec=match(bc_unique[,2],rownames(G9a_fc_mean))
bc_table_G9a = cbind(bc_unique[!is.na(G9a_match_vec),],G9a_fc_mean[G9a_match_vec[!is.na(G9a_match_vec)],])

ggplot(bc_table_KRAB, aes(x = as.factor(family), y=log(GPvsG_day2), colour=as.factor(family))) + geom_boxplot() + geom_point(position=position_jitter(0.2))
ggplot(bc_table_KRAB[grepl('LINE', bc_table_KRAB$family),], aes(x = as.factor(name), y=log(GPvsG_day2), colour=as.factor(family))) + geom_boxplot() + geom_point(position=position_jitter(0.2))

ggplot(bc_table_G9a, aes(x = as.factor(family), y=log(GPvsG_day2), colour=as.factor(family))) + geom_boxplot() + geom_point(position=position_jitter(0.2))
ggplot(bc_table_G9a[grepl('LINE', bc_table_G9a$family),], aes(x = as.factor(name), y=log(GPvsG_day2), colour=as.factor(family))) + geom_boxplot() + geom_point(position=position_jitter(0.2))