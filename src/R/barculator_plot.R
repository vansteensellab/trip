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
datapath <- "~/perldr/TRIP_data/mPGK-A/output/"
filename1 <- paste(datapath, "stats_table.txt", sep = "")
filename2 <- paste(datapath, "final_barcode_data_table.txt", sep = "")

# output file
output_file = paste(datapath, "summary_results.pdf", sep = "")

# a list to store all the plots
plots <- list()	

################### *******  Basic stats of sequencing reads **********###################

df <- read.table (file = filename1, header = T, stringsAsFactors = F)
rows = df$File
df <- melt(df, id.vars = "File")  # reshape the data
df$File <- factor(df$File, levels = rows )

p1 <- ggplot (data = df , aes(x = File, y = value, fill =variable)) + 
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


exp_df = trip[ , grepl("exp", colnames(trip)), drop = F]
cols_exp = colnames(exp_df)
exp_df = as.data.frame (sapply(1:ncol(exp_df), function(m) { (exp_df[,m] +1)/ sum(exp_df[,m] + 1)}))

norm = sum(grepl("norm", colnames(trip)))

if(norm) {
	norm_df = trip[ , grepl("norm", colnames(trip)), drop = F]
	cols_norm = colnames(norm_df)
	norm_df = as.data.frame (sapply(1:ncol(norm_df), function(m) { (norm_df[,m] + 1) / sum(norm_df[,m] + 1)}))
	colnames(norm_df) = cols_norm
	
	if (ncol(norm_df) == 1) {
		exp_df = as.data.frame (sapply(1:ncol(exp_df), function(m) { log2 (exp_df[,m]/norm_df[,1]) }))
	} else {
			exp_df = as.data.frame(sapply(1:ncol(exp_df), function(m) { log2 (exp_df[,m]/norm_df[,m]) }))
		}
	
} else {
	exp_df = as.data.frame (sapply(1:ncol(exp_df), function(m) { log2 (exp_df[,m])}))
}

exp_df <- as.data.frame(sapply(1:ncol(exp_df) , function(m) { scale (exp_df[,m], scale = F)}))

colnames(exp_df) = cols_exp	
df1 <- stack(exp_df)


p2 <- ggplot(data = df1, aes(x=values, fill = factor(ind, levels = cols_exp	))) + 
	geom_density(alpha = 0.3, colour = F) +
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
	
################### *******  Correlation plots **********################

if (ncol(exp_df) > 1) {
	mat <- cor(exp_df, method = "spearman")
	mat <- melt(mat)
	p3 <- ggplot(data = mat, aes(Var1, Var2, fill = value)) + geom_tile() + 
		scale_fill_gradient(name = "Correlation\n(Spearman)", low = "blue",  high = "yellow") +
		labs (title = "Correlation between different samples") +
		theme (plot.title = element_text (lineheight = 0.8, face = "bold", size = 12)) +
		#scale_fill_discrete(name = "Spearman Correlation") +
		scale_x_discrete (name = "") +
		scale_y_discrete (name = "")
	
}

plots [["p3"]] <- p3







#printing all the plots to pdf
ggsave(output_file, width = 21, height = 21, units = 'cm', do.call(marrangeGrob, c(plots, list(nrow=2, ncol=1))))


