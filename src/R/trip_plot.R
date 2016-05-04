# File: trip_plot.R
# Version: 0.1
# Author: Waseem ....
#
# Description: Generate plots on trip.pl output files

library(ggplot2)
library(plyr)
require(gridExtra)
require(reshape2)

args <- commandArgs(TRUE)
datapath = args[1]


# input files to be used
#datapath <- "~/perldr/TRIP_data/mPGK-A/data/replicate1/"
filename1 <- paste(datapath, "stats.txt", sep = "")
filename2 <- paste(datapath, "final_TRIP_data_table.txt", sep = "")

# output file
output_file = paste(datapath, "summary_results.pdf", sep = "")

# a list to store all the plots
plots <- list()	

################### *******  Basic stats of sequencing reads **********###################

df <- read.table (file = filename1, header = T, stringsAsFactors = F)
rows = df$Type
df <- melt(df, id.vars = "Type")  # reshape the data
df$Type <- factor(df$Type, levels = rows )

p1 <- ggplot (data = df , aes(x = Type, y = value, fill =variable)) + 
		geom_bar (stat = "identity", 
		position = position_dodge(), colour = F) + 
		scale_x_discrete (name = "Read Type") +
		scale_y_continuous (name = "No of Reads") +
		theme(panel.grid=element_blank()) +
		labs (title = "Summary of Reads", fill = "") +
		theme(plot.title = element_text (lineheight = 0.8, face = "bold", size = 12)) +
		theme (legend.position = "top")

plots [["p1"]] <- p1

################### *******  Density maps of normalized expression **********################

trip <- read.table (file = filename2, row.names=1, header = T, stringsAsFactors = F)
norm = as.vector(trip[,1])
norm = (norm + 1)/sum(norm)
exp_col_idx <- grep("exp", colnames (trip))
df1 <- trip[, exp_col_idx , drop = F]
cols1 = colnames(df1)

df1 <- as.data.frame (sapply(1:ncol(df1), function(m) { (df1[,m] + 1)/sum (df1[,m] + 1) }))
df1 <- as.data.frame (sapply(1:ncol(df1), function(m) { log2 ((df1[,m])/norm) }))
df1 <- as.data.frame(sapply(1:ncol(df1) , function(m) { scale (df1[,m], scale = F)}))
colnames(df1) = cols1
df1 <- stack(df1)




#cdf <- ddply(df1, .(ind), summarize, mean_counts = round(mean(values), 2))

p2 <- ggplot(data = df1, aes(x=values, fill = factor(ind, levels = cols1))) + geom_density(alpha = 0.3, colour = F) +
	# geom_vline(data = cdf, aes (xintercept = mean_counts, 
	# colour = factor(ind, levels = cols1)), linetype = "dashed", size = 0.5) +
	scale_fill_discrete(name = "Samples") +
	scale_x_continuous (name = "Expression") +
	scale_y_continuous (name = "Frequency") +
	labs (title = "Distribution of expression") +
	theme(plot.title = element_text (lineheight = 0.8, face = "bold", size = 12)) +
	theme(panel.grid=element_blank()) +
	theme (legend.position = "top")
	
plots [["p2"]] <- p2	
	


################### *******  Plots about mapping data **********###################

# first see if the mapping data is there

if (sum (grepl("chr", colnames(trip)))) {
	
	########## making a simple barplot about mapping ################
	
	lst <- list()
	total = nrow(trip)
	if (sum (grepl("_f$", colnames(trip)))) {
		mapped = total - sum(is.na (trip$chr_f))
		aligned = mapped - sum(trip$chr_f == "*", na.rm  = T)	
		map_stat_f <- data.frame (
			type = "forward",
			total = total,
			mapped = mapped,
			aligned = aligned,
			above.9 = aligned - sum(trip$freq1_f < 0.9, na.rm = T)
		)
		lst [["for"]] <- map_stat_f
	}
	if (sum (grepl("_r$", colnames(trip)))) {
		
		mapped = total - sum(is.na (trip$chr_r))
		aligned = mapped - sum(trip$chr_r == "*", na.rm  = T)	
		map_stat_r <- data.frame (
			type = "reverse",
			total = total,
			mapped = mapped,
			aligned = aligned,
			above.9 = aligned - sum(trip$freq1_r < 0.9, na.rm = T)
		)
		lst [["rev"]] <- map_stat_r
	}
	map_stat <- do.call(what = rbind, lapply(lst, "["))
	map_stat <- melt(map_stat, id.vars = "type")
	
	p_map_stat <- ggplot (data = map_stat , aes(x = variable, y = value, fill = type)) + 
		geom_bar (stat = "identity", 
		position = position_dodge(), colour = F) + 
		#scale_fill_discrete(name = "") +
		scale_x_discrete (name = "Mapping status") +
		scale_y_continuous (name = "Counts") +
		theme(panel.grid=element_blank()) +
		labs (title = "Summary of Mapping", fill = "") +
		theme(plot.title = element_text (lineheight = 0.8, face = "bold", size = 12)) +
		theme (legend.position = "top")
		
	plots [["p_map_stat"]] <- p_map_stat
		
	
	########## making plot about top location freq ##############
	
	df2 <- trip [ , !(grepl("exp|norm", colnames(trip))), drop = F]
	df2 <- df2 [ df2[,1] != "*" , , drop = F]
	df2 <- df2 [, grep("freq1", colnames(df2)), drop = F]
	df2 <- df2 [complete.cases(df2), , drop = F]
	freq_tbl <- data.frame()
	for (idx in 1:100) {
		for (jdx in 1:ncol(df2)) {
			freq_tbl[idx, jdx] <- (mean (df2[,jdx] >= idx/100)) * 100
		}
	}
	colnames(freq_tbl) <- colnames(df2)
	freq_tbl$thresh <- 1:100
	freq_tbl <- melt(freq_tbl, id.vars = "thresh")
	p_freq <- ggplot(data = freq_tbl, aes(x = thresh, y = value, colour = variable )) + geom_line(size = 1.2) +
		scale_x_continuous (name = "Threshold") +
		scale_y_continuous (name = "Percentage of barcodes") +
		#theme(panel.grid=element_blank()) +
		labs (title = "Top Location Freq Stats") +
		theme(plot.title = element_text (lineheight = 0.8, face = "bold", size = 12)) +
		theme (legend.position = "top") +
		labs (colour = "Mapping direction")
	plots [["p_freq"]] <- p_freq
	
}	
	
	
	
	
	
	
	
	# df2 <- trip [ , !(grepl("exp|norm", colnames(trip))), drop = F]
	# df2 <- df2 [ df2[,1] != "*" , ]
	# df2 <- df2 [, grep("freq", colnames(df2)), drop = F]
	# if (sum (grepl("_f$", colnames(df2)))) {
		# df_f <- df2 [, grep("_f$", colnames(df2)), drop = F]
		# df_f <- df_f[complete.cases(df_f), , drop = F]
		# cols_f = colnames(df_f)
		# df_f <- stack(df_f)
		# p3 <- ggplot(df_f, aes(x=values, fill = factor(ind, levels = cols_f))) + geom_histogram (binwidth = 0.02) +
				# scale_fill_discrete(name = "Samples") +
				# scale_x_continuous (name = "Expression") +
				# scale_y_continuous (name = "Counts") +
				# labs (title = "Distribution of frequencies of top locations") +
				# theme(plot.title = element_text (lineheight = 0.8, face = "bold", size = 12)) +
				# theme(panel.grid=element_blank()) +
				# theme (legend.position = "top")
		# plots[['p3']] <- p3
	# }
	# if (sum (grepl("_r$", colnames(df2)))) {
		# df_r <- df2 [, grep("_r$", colnames(df2)), drop = F]
		# df_r <- df_r[complete.cases(df_r), , drop = F]
		# cols_r = colnames(df_r)
		# df_r <- stack(df_r)
		# p4 <- ggplot(df_r, aes(x=values, fill = factor(ind, levels = cols_r))) + geom_histogram (binwidth = 0.02) +
				# scale_fill_discrete(name = "Samples") +
				# scale_x_continuous (name = "Expression") +
				# scale_y_continuous (name = "Counts") +
				# labs (title = "Distribution of frequencies of top locations") +
				# theme(plot.title = element_text (lineheight = 0.8, face = "bold", size = 12)) +
				# theme(panel.grid=element_blank()) +
				# theme (legend.position = "top")
		# plots[['p4']] <- p4
	# }	
# }

#printing all the plots to pdf
grobPlot = marrangeGrob(plots, nrow=2, ncol=2)
ggsave(output_file, width = 21, height = 21, units = 'cm',grobPlot )


