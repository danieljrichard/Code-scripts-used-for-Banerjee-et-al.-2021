###########
##Analyzing Salmon-quantified viral gene expression.
##Banerjee et al 2021
##
######

metadata_file <- "exp_METADATA_viral.csv"
outfix <- "Calu3_SARS_COV2_UPDATED_VIRAL"

##move around folders a bit.

library(data.table)
metadata <- as.data.frame(fread(metadata_file))
metadata$File <- paste0("VIRAL_RESULTS/", metadata$File)

metadata$File <- gsub("SARS2_CALU3", "CALU3_COV2_TRANS", metadata$File)
metadata$File <- gsub("CALU3_COV2_TRANS", "CALU3_COV2_TRANS_UPDATED", metadata$File)

metadata_curate <- metadata

colnames(metadata_curate) <- c("FILE", "TYPE", "TIME")

metadata_curate$TIME <- as.numeric(metadata_curate$TIME)
metadata_curate$Name <-  unlist(lapply(unique(metadata_curate$TIME), function(x) paste0(as.character(x), "H_", 1:length(which(metadata_curate$TIME == x)))))
metadata_curate$Name <- sapply(1:dim(metadata_curate)[1], function(x) paste0(metadata_curate$Name[x], "_", metadata_curate$ORIG_TYPE[x]))

metadata_curate$ORIG_TYPE <- as.factor(metadata_curate$ORIG_TYPE)
metadata_curate$ORIG_TYPE <- relevel(metadata_curate$ORIG_TYPE, "MOCK")


library(tximport)
library(DESeq2)
metadata_curate <- metadata_curate[metadata_curate$ORIG_TYPE != "MOCK",]

tx_data_MAP <- tximport(as.character(metadata_curate$FILE), type="salmon", ignoreTxVersion=TRUE, txOut = T)

metadata_curate$TIME <- as.factor(metadata_curate$TIME)
de_data_MAP <- DESeqDataSetFromTximport(tx_data_MAP, colData = metadata_curate, design = ~ TIME)

##Applying filtering for lowly-expressed genes.
filter_de_data <- de_data_MAP[rowSums(counts(de_data_MAP) > 5) > 3,]
filter_de_data <- filter_de_data[order(-rowSums(counts(filter_de_data))),]

outfix <- "COV2_inf_CALU3_JUST_INFECTED_UPD_PCA"

test_vst <- vst(filter_de_data, nsub = 10)
pdf(paste0(outfix, "_sample_PCA_plots.pdf"), width=10, height=10)
plotPCA(test_vst, intgroup="TIME")
dev.off()

sampleDists <- dist(t(assay(test_vst)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- metadata_curate$Name
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

rownames(sampleDistMatrix) <- metadata_curate$Name

pdf(paste0(outfix, "_prettier_sample_PCA_plots2.pdf"), width=10, height=10)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

salmon_DE <- DESeq(filter_de_data)
all_by_all <- list()

dir.create("VIRAL_PROCESSING")
setwd("VIRAL_PROCESSING")

pdf(paste0(outfix, "_treament_comparisons.pdf"), width = 10, height = 12)
table_frame <- list()
table_frame_MAP <- list()

x <- 1
sig_hits <- results(salmon_DE, name = "TIME")
sig_hits <- sig_hits[complete.cases(sig_hits),]
plotMA(sig_hits, ylim=c(-5, 5))
fdr_pass_MAP <- sig_hits[which(sig_hits$padj < 0.05),]
table_frame_MAP[[x]] <- fdr_pass_MAP

dev.off()

library(rtracklayer)

sig_genes <- rownames(fdr_pass_MAP)

##read in viral gene data.
test <- readGFF("GCF_009858895.2_ASM985889v3_genomic.gff")
test_clean <- test[which(!is.na(test$locus_tag)),]
test_clean <- as.data.frame(test_clean)
test_names <- test_clean$ID
test_names_clean <- sapply(test_names, function(x) unlist(strsplit(x, "YP_"))[2])
test_names_clean[is.na(test_names_clean)] <- "NA"
test_grep <- lapply(test_names_clean, function(x) grep(x, sig_genes))
good_bool <- sapply(test_grep, function(x) length(x) == 1)
subset_test <- test_clean[which(good_bool),]
subset_test$sig_genes <- sig_genes[unlist(test_grep)]

subset_test <- subset_test[!duplicated(subset_test$sig_genes),]
subset_rebool <- lapply(sig_genes, function(x) which(subset_test$sig_genes == x))

subset_test <- subset_test[unlist(subset_rebool),]

write.csv(subset_test, "SIGNIFICANT_GENES_ORF_MAPPING.csv")
write.csv(subset_test, "TRULY_SIGNIFICANT_GENES_ORF_MAPPING.csv")

names(table_frame_MAP) <- names(all_by_all)

saveRDS(salmon_DE, paste0(outfix, "_salmon_DESEQ_obj_BIOMART.rds"))

try(dir.create("finished_salmon_proc"))
setwd("finished_salmon_proc")

out_mat_count <- as.data.frame(counts(filter_de_data)) ##Output the quants only

colnames(out_mat_count) <- filter_de_data$Name

out_mat_count$transcript <- rownames(out_mat_count)

library(data.table)
fwrite(out_mat_count, paste0(outfix, "_salmon_RAW_counts_matrix_biomart_genenames.csv"), row.names=F)

##Now cutting out 'raw' quants for just the significant genes...
sub_thing <- lapply(subset_test$sig_genes, function(x) which(out_mat_count$transcript == x))
out_mat_subset <- out_mat_count[unlist(sub_thing),]
out_mat_subset <- data.frame(out_mat_subset, subset_test)
fwrite(out_mat_subset, paste0(outfix, "_significant_transcripts_RAW_COUNTS_genenames.csv"), row.names=F)

fpkm_vals <- as.data.frame(fpkm(filter_de_data, robust = T))
colnames(fpkm_vals) <- filter_de_data$Name

fpkm_vals$transcript <- rownames(fpkm_vals)

fwrite(fpkm_vals, paste0(outfix, "_salmon_FPKM_counts_matrix_genenames.csv"), row.names=F)

sub_thing_fpkm <- lapply(subset_test$sig_genes, function(x) which(fpkm_vals$transcript == x))
fpkm_out_mat_subset <- fpkm_vals[unlist(sub_thing_fpkm),]
fpkm_out_mat_subset <- data.frame(fpkm_out_mat_subset, subset_test)
fwrite(fpkm_out_mat_subset, paste0(outfix, "_significant_transcripts_FPKM_counts_matrix_genenames.csv"), row.names=F)

salmon_DE_MAP_sizes <- estimateSizeFactors(salmon_DE) ##Correcting for library size

salmon_DE_MAP_counts <- as.data.frame(counts(salmon_DE_MAP_sizes, normalized=TRUE))
colnames(salmon_DE_MAP_counts) <-  as.character(salmon_DE$Name)

salmon_DE_MAP_counts$gene <- rownames(salmon_DE_MAP_counts) ##Assuming that order of DE_map doesn't change from data_MAP
salmon_DE_MAP_counts$transcript <- rownames(salmon_DE_MAP_counts)

salmon_DE_mat <- do.call("cbind", lapply(1:(dim(salmon_DE_MAP_counts)[2] - 2), function(x) as.numeric(salmon_DE_MAP_counts[,x])))
library(ComplexHeatmap)


z_score_normed_counts <- t(scale(t(salmon_DE_mat), center=TRUE, scale=TRUE))
rownames(z_score_normed_counts) <- salmon_DE_MAP_counts$transcript

colnames(z_score_normed_counts) <- salmon_DE$Name

check_na <- apply(z_score_normed_counts, 1, function(x) length(which(is.na(x))) != 0)
z_score_normed_counts <- z_score_normed_counts[!check_na,]

thresh_map = Heatmap(z_score_normed_counts, name= paste0(outfix, "_all_samples_all_genes_salmon_heatmap"),
 cluster_columns=FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=F)

pdf(paste0(outfix, "_all_samples_all_genes_salmon_heatmap.pdf"), width=10, height=12)
print(thresh_map)
dev.off()

startdir <- getwd()

x <- 1
curr_dir <- "VIRUS"
try(dir.create(curr_dir))
try(setwd(curr_dir))

curr_rel_cols <- 1:dim(metadata_curate)[1]
curr_norm_counts <- salmon_DE_MAP_counts[, unlist(curr_rel_cols)]
	
curr_z_score_normed_counts <- t(scale(t(curr_norm_counts), center=TRUE, scale=TRUE))
colnames(curr_z_score_normed_counts) <- colnames(curr_norm_counts)

curr_check_na <- apply(curr_z_score_normed_counts, 1, function(x) length(which(is.na(x))) != 0)
curr_z_score_normed_counts <- curr_z_score_normed_counts[!curr_check_na,]

curr_norm_counts_clean <- curr_norm_counts[!curr_check_na,]

curr_thresh_map = Heatmap(curr_z_score_normed_counts, name = "POST_INFECTION",
 cluster_columns=TRUE, row_names_gp = gpar(fontsize = 6), show_row_names=F)

pdf("UPDATED_CALU3_COV2_POST_INFECTION_allgenes_salmon_heatmap.pdf", width = 10, height = 12)
print(curr_thresh_map)
dev.off()

curr_sig_results <- table_frame_MAP[[x]]

sig_matchback <- lapply(rownames(curr_sig_results), function(x) which(rownames(curr_z_score_normed_counts) == x))
over_check <- sapply(sig_matchback, function(x) length(x) > 1)
if(sum(over_check) != 0) {
	print("overmatch")
	browser()
}	
under_check <- sapply(sig_matchback, function(x) length(x) == 0)
if(sum(under_check) != 0) {
	print("undermatch?")
	browser()
}

z_score_sig_set <- curr_z_score_normed_counts[unlist(sig_matchback),]
sig_curr_norm_counts_clean <- curr_norm_counts_clean[unlist(sig_matchback),]

rownames(z_score_sig_set) <- subset_test$product

z_score_sig_set <- z_score_sig_set[, order(colnames(z_score_sig_set))]

sig_curr_thresh_map = Heatmap(z_score_sig_set, name= "COV2_SIG005_genes_salmon_heatmap",
 cluster_columns=FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=T)

pdf("UPD_CALU3_COV2_POSTINFECTION_SIGNIFICANT_VIRAL_GENES_HEATMAP.pdf", width = 10, height = 12)
print(sig_curr_thresh_map)
dev.off()

fwrite(sig_curr_norm_counts_clean, paste0(outfix, "_SIG005_salmon_normalized_counts_matrix.csv"), row.names=F)

setwd(startdir)

filter_normed_salmon_MAP_counts <- salmon_DE_MAP_counts


gene_sets <- list()


	setwd("VIRUS")
	x <- 1
	names(table_frame_MAP)[1] <- "POST_COVID_INFECTION"
	
	write.csv(table_frame_MAP[[x]], paste0("significant_comparisons_BIOMART_", names(table_frame_MAP)[x], "_FDR005.csv"), row.names=F)
	writeLines(as.character(table_frame_MAP[[x]]$gene_rep), paste0("significant_comparisons_BIOMART_", names(table_frame_MAP)[x], "_FDR005_GENES.txt"))
	
	writeLines(as.character(rownames(table_frame_MAP[[x]])[which(table_frame_MAP[[x]]$log2FoldChange >= 0)]), paste0("significant_comparisons_BIOMART_", "UP_IN_", "COVID_infection", "_FDR005_GENES.txt"))
	writeLines(as.character(rownames(table_frame_MAP[[x]])[which(table_frame_MAP[[x]]$log2FoldChange < 0)]), paste0("tissue_fragment_significant_comparisons_BIOMART_", "DOWN_IN_", "COVID_infection", "_FDR005_GENES.txt"))

	curr_list <- list()
	curr_list[[1]] <- as.character(table_frame_MAP[[x]]$gene_rep[which(table_frame_MAP[[x]]$log2FoldChange >= 0)])
	curr_list[[2]] <- as.character(table_frame_MAP[[x]]$gene_rep[which(table_frame_MAP[[x]]$log2FoldChange < 0)])
	setwd("..")

big_sig_hits <- results(salmon_DE, alpha=0.05, name = "TIME")

if(dim(big_sig_hits)[1] != dim(filter_normed_salmon_MAP_counts)[1]) {
	print("error to be fixed!")
	big_sig_hits$gene <- mcols(salmon_DE)$gene_sym
	big_sig_hits <- big_sig_hits[!duplicated(big_sig_hits$gene),]
	filter_normed_salmon_MAP_counts <- filter_normed_salmon_MAP_counts[!duplicated(filter_normed_salmon_MAP_counts$gene),]
	filter_normed_salmon_MAP_counts <- data.table(filter_normed_salmon_MAP_counts)
	subset <- intersect(big_sig_hits$gene, filter_normed_salmon_MAP_counts$gene)
	big_sig_hits <- data.table(as.data.frame(big_sig_hits))
	setkey(filter_normed_salmon_MAP_counts, gene)
	setkey(big_sig_hits, gene)
	filter_normed_salmon_MAP_counts <- filter_normed_salmon_MAP_counts[.(subset)]
	big_sig_hits <- big_sig_hits[.(subset)]
	filter_normed_salmon_MAP_counts <- as.data.frame(filter_normed_salmon_MAP_counts)
	big_sig_hits <- as.data.frame(big_sig_hits)
	browser()
}

big_curr_output <- data.frame(filter_normed_salmon_MAP_counts, big_sig_hits)

##similar to the host expression analyses, outputting a large CSV table for timecourse plotting purposes:
curr_name <- paste0(outfix, "_salmon_counts_FILTERED_NORMALIZED_PERSAMPLE_CROSS_BATCH_MATRIX.csv")

fwrite(big_curr_output, curr_name, row.names=F, sep=",")

##Cutdown.

big_sig_output <- big_curr_output[big_curr_output$padj < 0.05,]
big_sig_output <- big_sig_output[!is.na(big_sig_output$padj),]

big_sig_output$trans <- rownames(big_sig_output)
big_sig_output <- data.frame(big_sig_output, subset_test)

fwrite(big_sig_output, gsub(".csv", "_GENE_NAMES_significant.csv", curr_name), row.names=F, sep=",")


}



gen_timecourse_plot_virus_ADV <- function(csv, geneset = NULL, outfix, ordered = F, fix_axis = F, do_log = F) {
	rna <- read.csv(csv)
	
	if(!is.null(geneset)) {
	gene_list <- readLines(geneset) 
	
	relevant_lines <- lapply(gene_list, function(x) which(rna$gene == x))
	relevant_dat <- rna[unlist(relevant_lines),]
}else{
	relevant_dat <- rna
}
	
	if(ordered) {
		print("sorting by start")
		relevant_dat <- relevant_dat[order(relevant_dat$start),]	
	}
	
		if(fix_axis) {
			
		grab_thing <- function(line) {
				colnames(line) <- gsub("X", "", colnames(line))
		rel_line <- line[, 1:(which(colnames(line) == "gene") - 1)]
		if(identical(grep("gene", colnames(line)), integer(0))) {
		rel_line <- line[, 1:(which(colnames(line) == "transcript") - 1)]	
		}
		time_points <- sapply(colnames(rel_line), function(x) unlist(strsplit(x, "H"))[1])
		type <- sapply(colnames(rel_line), function(x) unlist(strsplit(x, "_"))[length(unlist(strsplit(x, "_")))])
		
		super_set_mean <- list()
		for (x in unique(time_points)) {
			rel_data <- rel_line[, which(time_points == x)]
			rel_data_frame <- data.frame(val = t(rel_data), type = type[which(time_points == x)])
			colnames(rel_data_frame)[1] <- "val"
			rel_data_frame$type <- as.factor(rel_data_frame$type)
			rel_data_frame$TIME <- rep(x, dim(rel_data_frame)[1])
			
			if(do_log) {
				rel_data_frame$val[rel_data_frame$val == 0] <- 1
				rel_data_frame$val <- log10(rel_data_frame$val)
				
			}
			
			rel_data_frame_mean <- sapply(unique(rel_data_frame$type), function(e) mean(rel_data_frame$val[which(rel_data_frame$type == e)]))
			rel_data_frame_sd <- sapply(unique(rel_data_frame$type), function(e) sd(rel_data_frame$val[which(rel_data_frame$type == e)]))
			
			rel_mean <- data.frame(mean = rel_data_frame_mean, sd = rel_data_frame_sd, type = unique(rel_data_frame$type), TIME = x)
			
			
			rel_mean$TIME <- gsub("4", "6", rel_mean$TIME)
			rel_mean$TIME <- gsub("5", "12", rel_mean$TIME)
			rel_mean$TIME <- paste0(rel_mean$TIME, "H")
			
			colnames(rel_mean)[1:2] <- paste0(c("Mean INF", "SD INF"), " ", rel_mean$TIME)
			
			super_set_mean[[x]] <- rel_mean
		}
		
		super_duper_set <- do.call("cbind", lapply(super_set_mean, function(x) x[, 1:2]))
		colnames(super_duper_set) <- sapply(colnames(super_duper_set), function(x) unlist(strsplit(x, ".", fixed=T))[2])
		
	return(super_duper_set)
	}
	
		super_stuff <- lapply(1:dim(relevant_dat)[1], function(x) grab_thing(relevant_dat[x,]))
		
		super_squash <- do.call("rbind", super_stuff)
				
		super_squash$gene = relevant_dat$gene
		super_squash$transcript = relevant_dat$transcript
		
		write.csv(super_squash, paste0(outfix, "_summarized_read_count_data.csv"), row.names=F, quote = F)
		
		ylim_range <- c((min(super_squash$mean) - super_squash$sd[which.min(super_squash$mean)]), (max(super_squash$mean) + super_squash$sd[which.max(super_squash$mean)]))
				
}
	
	plot_thing <- function(line,  END = FALSE) {
		colnames(line) <- gsub("X", "", colnames(line))
		rel_line <- line[, 1:(which(colnames(line) == "gene") - 1)]
	
		if(identical(grep("gene", colnames(line)), integer(0))) {
		rel_line <- line[, 1:(which(colnames(line) == "transcript") - 1)]	
		}
		time_points <- sapply(colnames(rel_line), function(x) unlist(strsplit(x, "H"))[1])
		type <- sapply(colnames(rel_line), function(x) unlist(strsplit(x, "_"))[length(unlist(strsplit(x, "_")))])
		
		##UPD: FIXED AXES:
		
		library(ggplot2)
		super_set <- list()
		super_set_mean <- list()
		for (x in unique(time_points)) {
			rel_data <- rel_line[, which(time_points == x)]
			rel_data_frame <- data.frame(val = t(rel_data), type = type[which(time_points == x)])
			colnames(rel_data_frame)[1] <- "val"
			rel_data_frame$type <- as.factor(rel_data_frame$type)
			rel_data_frame$TIME <- rep(x, dim(rel_data_frame)[1])
			
			if(do_log) {
				rel_data_frame$val[rel_data_frame$val == 0] <- 1
				rel_data_frame$val <- log2(rel_data_frame$val)			
			}
			
			rel_data_frame_mean <- sapply(unique(rel_data_frame$type), function(e) mean(rel_data_frame$val[which(rel_data_frame$type == e)]))
			rel_data_frame_sd <- sapply(unique(rel_data_frame$type), function(e) sd(rel_data_frame$val[which(rel_data_frame$type == e)]))
			
			rel_mean <- data.frame(mean = rel_data_frame_mean, sd = rel_data_frame_sd, type = unique(rel_data_frame$type), TIME = x)
			super_set[[x]] <- rel_data_frame
			super_set_mean[[x]] <- rel_mean
		}
		super_frame <- do.call("rbind", super_set)
		
		super_frame_mean <- do.call("rbind", super_set_mean)

		super_frame_mean$TIME <- as.numeric(as.character(super_frame_mean$TIME))
		super_frame_mean <- super_frame_mean[order(super_frame_mean$TIME),]

		super_frame_mean$TIME <- paste0(super_frame_mean$TIME, "H")
		super_frame_mean$TIME[super_frame_mean$TIME == "4H"] <- "6H"
		super_frame_mean$TIME[super_frame_mean$TIME == "5H"] <- "12H"
						
		library(forcats)
		p <- ggplot(super_frame_mean, aes(x = fct_inorder(TIME), y = mean, group = type, color = type)) + geom_line() + geom_point() 
		p <- p + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                 position=position_dodge(0.05)) + ggtitle(as.character(line$product)) + labs(y = "Normalized Counts ") + theme_classic()
 
 p <- p + theme(legend.position = "none")
 
 if(END) {
	 p <- p + theme(axis.title.x = element_blank())
 }else{
	p <- p + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
 }
  
	if(fix_axis) {
		print("fixing axis")
		p <- p + coord_cartesian(ylim = ylim_range)
		
	}
 print(as.character(line$product))
 
     p <- p + theme(axis.title.y = element_blank())
	p <- p + theme(text = element_text(size = 12))
	p <- p + theme(plot.title = element_text(size = 12))
        return(p)
	}
	
	all_genes <- lapply(1:dim(relevant_dat)[1], function(x) plot_thing(relevant_dat[x,]))
	
	all_genes <- lapply(1:dim(relevant_dat)[1], function(x) if(x == 11 || x == 12) {plot_thing(relevant_dat[x,], T)} else {plot_thing(relevant_dat[x,]) })
	
	library(gridExtra)
	
	pdf(paste0(outfix, "_timecourse.pdf"), width = 11, height = 8)
	library(ggplot2)
	library(grid)
	if(do_log) {
	title <- textGrob("Log2 Transcript Counts", rot = 90, vjust = 1, gp = gpar(fontsize = 12))#, gp = gpar(fontfamily = "Arial"))
	bot = textGrob("Time Post-Infection", gp = gpar(fontsize = 12))#, gp = gpar(fontfamily = "Arial")) 
}else{
	title <- textGrob("Normalized Read Counts", rot = 90, vjust = 1)
}
	print(grid.arrange(grobs = all_genes, left = title, ncol = 2, bottom = bot))

	dev.off()
	
}


##takes the output from above, along with an optional geneset.
gen_timecourse_plot_virus_ADV("COV2_inf_CALU3_JUST_INFECTED_UPD_salmon_counts_FILTERED_NORMALIZED_PERSAMPLE_CROSS_BATCH_MATRIX_GENE_NAMES_significant.csv",
	outfix = "UPDATED_COV2_significant_transcripts_normalized")

