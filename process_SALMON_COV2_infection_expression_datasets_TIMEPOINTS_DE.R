###########
##Processing Salmon Output for Timepoints DE Analysis
##Banerjee et al 2021
##
##


gtf_file <- "Homo_sapiens.GRCh37.67.gtf"

gen_tx2gene <- function(gtf_file) {
	
	library(GenomicFeatures)
	
	txdb <- makeTxDbFromGFF(gtf_file)
	
#https://support.bioconductor.org/p/81954/
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]
a<- gsub("\\..*","",tx2gene[,1])
b<- gsub("\\..*","",tx2gene[,2])
c<-cbind(a,b)
colnames(c)=colnames(tx2gene)
tx2gene <- as.data.frame(c)
#head(tx2gene)

saveRDS(tx2gene, "hg19_gene_annotations_mapped_genesymbols.rds")

}


gtf_rds <- "hg19_gene_annotations_mapped_genesymbols.rds"

metadata_file <- "UPDATED_METADATA.csv"
outfix <- "Calu3_SARS_COV2"

library(data.table)
metadata <- as.data.frame(fread(metadata_file))
##Adjust paths slightly:
metadata$File <- paste0(getwd(), "/RESULTS/", metadata$File)
metadata$File <- gsub("_CALU3", "_CALU3_UPDATED", metadata$File)

metadata_curate <- metadata

colnames(metadata_curate) <- c("FILE", "TYPE", "TIME")

times <- metadata_curate$TIME
range <- sort(unique(times))
time_scale <- sapply(times, function(x) which(range == x))
time_scale <- time_scale - 1 ##Just to set things to 0

metadata_curate$ORIG_TYPE <- metadata_curate$TYPE

metadata_curate$TIME <- time_scale
metadata_curate$Name <- unlist(lapply(unique(metadata_curate$orig_time), function(x) paste0(as.character(x), 1:length(which(metadata_curate$orig_time == x)))))

metadata_curate$TIME <- as.numeric(metadata_curate$TIME)
metadata_curate$Name <-  unlist(lapply(unique(metadata_curate$TIME), function(x) paste0(as.character(x), "H_", 1:length(which(metadata_curate$TIME == x)))))
metadata_curate$Name <- sapply(1:dim(metadata_curate)[1], function(x) paste0(metadata_curate$Name[x], "_", metadata_curate$ORIG_TYPE[x]))

metadata_curate$ORIG_TYPE <- as.factor(metadata_curate$ORIG_TYPE)
metadata_curate$ORIG_TYPE <- relevel(metadata_curate$ORIG_TYPE, "MOCK")

library("DESeq2")

if(is.null(gtf_rds)) {
library('biomaRt')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest")

map_to_symbol <- getBM(filters = "ensembl_gene_id", attributes = c("hgnc_symbol", "ensembl_gene_id"), values = tx2gene$GENEID, mart=mart)

map_df <- map_to_symbol[which(map_to_symbol[,1] != ""),]
map_back <- lapply(map_df[,2], function(x) which(tx2gene$GENEID == x))

super_check <- sapply(map_back, length)

gene_rep <- lapply(1:length(map_back), function(x) rep(map_df[x,1], super_check[x]))
ens_gene_rep <- lapply(1:length(map_back), function(x) rep(map_df[x,2], super_check[x]))

final_out_frame <- data.frame(tx2gene[unlist(map_back),], gene_rep = unlist(gene_rep), ens_gene_rep = unlist(ens_gene_rep))
final_check <- sapply(1:dim(final_out_frame)[1], function(x) final_out_frame$GENEID[x] == final_out_frame$ens_gene_rep[x])

if(sum(final_check) != dim(final_out_frame)[1]) {
print("messed up")
browser()
}

saveRDS(final_out_frame, "hg19_gene_annotations_mapped_genesymbols.rds")
}else{
	library(biomaRt)
	final_out_frame <- readRDS(gtf_rds)
}


start_super <- getwd()
super_curate <- metadata_curate

for (x in unique(super_curate$TIME)) {
	metadata_curate <- super_curate[super_curate$TIME == x,]
	do_diff <- T
	SUPER_DIFF <- T

# BiocManager::install("tximport")
library(tximport)

tx_data_MAP <- tximport(as.character(metadata_curate$FILE), type="salmon", ignoreTxVersion=TRUE, tx2gene = final_out_frame)

de_data_MAP <- DESeqDataSetFromTximport(tx_data_MAP, colData = metadata_curate, design = ~ORIG_TYPE)

try(dir.create(paste0("CLEAN_TIME_POINT_", as.character(x))))
setwd(paste0("CLEAN_TIME_POINT_", as.character(x)))

map_back_genes_dedata <- final_out_frame[unlist(sapply(rownames(de_data_MAP), function(x) which(final_out_frame$ens_gene_rep == x)[1])), c("GENEID", "gene_rep")]
##Trying to add gene names into de object to save map-backs later:

mcols(de_data_MAP) <- cbind(mcols(de_data_MAP), gene_sym = map_back_genes_dedata$gene_rep)

##Applying filtering for lowly-expressed genes. At least counts of 5 in more than three samples.
filter_de_data <- de_data_MAP[rowSums(counts(de_data_MAP) > 5) > 3,]
filter_de_data <- filter_de_data[order(-rowSums(counts(filter_de_data))),]
filter_de_data <- filter_de_data[!duplicated(mcols(filter_de_data)$gene_sym),]

##General PCA:

test_vst <- vst(filter_de_data)
pdf(paste0(outfix, "_sample_PCA_plots.pdf"), width=10, height=10)
plotPCA(test_vst, intgroup="ORIG_TYPE")
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

treats <- as.character(unique(metadata_curate$ORIG_TYPE))
for (x in 1:(length(treats)-1)) {
	for (y in (x+1):length(treats)) {
		curr_cond <- c("ORIG_TYPE", treats[x], treats[y])
		all_by_all[[paste0(treats[x], "_", treats[y])]] <- curr_cond
	}
}

filter_biomart <- T ##Otherwise not implemented.

sig_hits <- results(salmon_DE, alpha=0.05, contrast = all_by_all[[1]])
sig_hits <- sig_hits[complete.cases(sig_hits),]
fdr_pass <- sig_hits[which(sig_hits$padj < 0.05),]
map_back_genes <- final_out_frame[unlist(sapply(rownames(fdr_pass), function(x) which(final_out_frame$ens_gene_rep == x)[1])), c("GENEID", "gene_rep")]
out_frame <- data.frame(fdr_pass, map_back_genes)
table_frame_MAP <- list()

table_frame_MAP[[1]] <- out_frame

##Save this out
if(filter_biomart) {
saveRDS(salmon_DE, paste0(outfix, "_salmon_DESEQ_obj_BIOMART.rds"))
}else{
saveRDS(salmon_DE, paste0(outfix, "_salmon_DESEQ_obj.rds"))
}

try(dir.create("finished_salmon_proc"))
setwd("finished_salmon_proc")

out_mat_count <- as.data.frame(counts(filter_de_data)) 

colnames(out_mat_count) <- filter_de_data$Name

out_mat_count$gene <- mcols(filter_de_data)$gene_sym
out_mat_count$ensembl <- rownames(out_mat_count)

library(data.table)
fwrite(out_mat_count, paste0(outfix, "_salmon_RAW_counts_matrix_biomart_genenames.csv"), row.names=F)

fpkm_vals <- as.data.frame(fpkm(filter_de_data, robust = T))
colnames(fpkm_vals) <- filter_de_data$Name

fpkm_vals$gene <- mcols(filter_de_data)$gene_sym
fpkm_vals$ensembl <- rownames(fpkm_vals)
fwrite(fpkm_vals, paste0(outfix, "_salmon_FPKM_counts_matrix_genenames.csv"), row.names=F)

salmon_DE_MAP_sizes <- estimateSizeFactors(salmon_DE) ##Correcting for library size

salmon_DE_MAP_counts <- as.data.frame(counts(salmon_DE_MAP_sizes, normalized=TRUE))
colnames(salmon_DE_MAP_counts) <-  as.character(salmon_DE$Name)

salmon_DE_MAP_counts$gene <- mcols(filter_de_data)$gene_sym
salmon_DE_MAP_counts$ensembl <- rownames(salmon_DE_MAP_counts)

salmon_DE_mat <- do.call("cbind", lapply(1:(dim(salmon_DE_MAP_counts)[2] - 2), function(x) as.numeric(salmon_DE_MAP_counts[,x])))

##output heatmap.
library(ComplexHeatmap)

z_score_normed_counts <- t(scale(t(salmon_DE_mat), center=TRUE, scale=TRUE))
rownames(z_score_normed_counts) <- salmon_DE_MAP_counts$ensembl

colnames(z_score_normed_counts) <- salmon_DE$Name

check_na <- apply(z_score_normed_counts, 1, function(x) length(which(is.na(x))) != 0)
z_score_normed_counts <- z_score_normed_counts[!check_na,]

thresh_map = Heatmap(z_score_normed_counts, name= paste0(outfix, "_all_samples_all_genes_salmon_heatmap"),
 cluster_columns=TRUE, row_names_gp = gpar(fontsize = 6), show_row_names=F)

pdf(paste0(outfix, "_all_samples_all_genes_salmon_heatmap.pdf"), width=10, height=12)
print(thresh_map)
dev.off()


##Let's do it all-by-all;
startdir <- getwd()

for (x in 1:length(all_by_all)) {
	
	curr_dir <- paste0(all_by_all[[x]][2], "-", all_by_all[[x]][3])
	try(dir.create(curr_dir))
	try(setwd(curr_dir))
	
	curr_rel_cols <- sapply(all_by_all[[x]][2:3], function(y) grep(y, colnames(salmon_DE_MAP_counts)))
	curr_rel_cols <- 1:dim(metadata_curate)[1]
	curr_norm_counts <- salmon_DE_MAP_counts[, unlist(curr_rel_cols)]
	
	curr_z_score_normed_counts <- t(scale(t(curr_norm_counts), center=TRUE, scale=TRUE))
	colnames(curr_z_score_normed_counts) <- colnames(curr_norm_counts)

curr_check_na <- apply(curr_z_score_normed_counts, 1, function(x) length(which(is.na(x))) != 0)
curr_z_score_normed_counts <- curr_z_score_normed_counts[!curr_check_na,]

curr_norm_counts_clean <- curr_norm_counts[!curr_check_na,]

curr_thresh_map = Heatmap(curr_z_score_normed_counts, name= paste0(outfix, "_", all_by_all[[x]][2], "-vs-", all_by_all[[x]][3], "_all_genes_salmon_heatmap"),
 cluster_columns=TRUE, row_names_gp = gpar(fontsize = 6), show_row_names=F)

pdf(paste0(outfix, "_", all_by_all[[x]][2], "-vs-", all_by_all[[x]][3], "_all_genes_salmon_heatmap.pdf"), width=10, height=12)
print(curr_thresh_map)
dev.off()

	curr_sig_results <- table_frame_MAP[[x]]
	##We'll take the cleaned frame. If the significant genes didn't pass NA filtering above, they're not going to pass it below.
	sig_matchback <- lapply(curr_sig_results$GENEID, function(x) which(rownames(curr_z_score_normed_counts) == x))
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

	rownames(z_score_sig_set) <- as.character(curr_sig_results$gene_rep)
	sig_curr_norm_counts_clean$gene <- as.character(curr_sig_results$gene_rep)
	sig_curr_norm_counts_clean$ensembl <- rownames(sig_curr_norm_counts_clean)
	
sig_curr_thresh_map = Heatmap(z_score_sig_set, name= paste0(outfix, "_", all_by_all[[x]][2], "-vs-", all_by_all[[x]][3], "_SIG005_genes_salmon_heatmap"),
 cluster_columns=TRUE, row_names_gp = gpar(fontsize = 6), show_row_names=F)
	
pdf(paste0(outfix, "_", all_by_all[[x]][2], "-vs-", all_by_all[[x]][3], "_SIG005_salmon_heatmap2.pdf"), width=10, height=12)
print(sig_curr_thresh_map)
dev.off()

fwrite(sig_curr_norm_counts_clean, paste0(outfix, "_", all_by_all[[x]][2], "-vs-", all_by_all[[x]][3], "_SIG005_salmon_normalized_counts_matrix.csv"), row.names=F)

setwd(startdir)
print(paste0(all_by_all[[x]][2], "-vs-", all_by_all[[x]][3]))
}

##Applying an additional filter to remove genes for which adjusted counts are quite low.

filter_bool <- rep(TRUE, dim(salmon_DE_MAP_counts)[1])
big_check <-  apply(salmon_DE_MAP_counts[, 1:(which(grepl("gene", colnames(salmon_DE_MAP_counts))) - 1)], 1, function(x) length(which(x >= 1)))
big_avg <-  apply(salmon_DE_MAP_counts[, 1:(which(grepl("gene", colnames(salmon_DE_MAP_counts))) - 1)], 1, function(x) mean(x))
##These should be overlapping, but still doing:
filter_bool[which(big_check < 3)] <- FALSE ##Not expressed in at least 3 of 12 samples.
filter_bool[which(big_avg < 1)] <- FALSE ##Average expression is less than one adjusted count, should largely overlap above check.

filter_normed_salmon_MAP_counts <- salmon_DE_MAP_counts[filter_bool,]
colnames(filter_normed_salmon_MAP_counts)[1:length(as.character(salmon_DE$Name))] <- as.character(salmon_DE$Name)

fwrite(filter_normed_salmon_MAP_counts, paste0(outfix, "_salmon_counts_matrix_biomart_genenames_FILTERED_NORMALIZED.csv"), row.names=F, sep=",")

writeLines(as.character(filter_normed_salmon_MAP_counts$gene), paste0(outfix, "_salmon_counts_matrix_biomart_genenames_FILTERED_NORMALIZED_genes.txt"))

gene_sets <- list()

for (x in 1:length(table_frame_MAP)) {
	curr_dir <- paste0(all_by_all[[x]][2], "-", all_by_all[[x]][3])
	setwd(curr_dir)	
	
	write.csv(table_frame_MAP[[x]], paste0("significant_comparisons_BIOMART_", names(table_frame_MAP)[x], "_FDR005.csv"), row.names=F)
	writeLines(as.character(table_frame_MAP[[x]]$gene_rep), paste0("significant_comparisons_BIOMART_", names(table_frame_MAP)[x], "_FDR005_GENES.txt"))
	writeLines(as.character(table_frame_MAP[[x]]$gene_rep[which(table_frame_MAP[[x]]$log2FoldChange >= 0)]), paste0("significant_comparisons_BIOMART_", "UP_IN_", all_by_all[[x]][2], "_vs_", all_by_all[[x]][3], "_FDR005_GENES.txt"))
	writeLines(as.character(table_frame_MAP[[x]]$gene_rep[which(table_frame_MAP[[x]]$log2FoldChange < 0)]), paste0("significant_comparisons_BIOMART_", "UP_IN_", all_by_all[[x]][3], "_vs_", all_by_all[[x]][2], "_FDR005_GENES.txt"))

	curr_list <- list()
	curr_list[[1]] <- as.character(table_frame_MAP[[x]]$gene_rep[which(table_frame_MAP[[x]]$log2FoldChange >= 0)])
	curr_list[[2]] <- as.character(table_frame_MAP[[x]]$gene_rep[which(table_frame_MAP[[x]]$log2FoldChange < 0)])
	gene_sets[[curr_dir]] <- curr_list
	
print(x)
	setwd(startdir)

big_sig_hits <- results(salmon_DE, alpha=0.05, contrast = all_by_all[[1]])

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

big_averages_frame <- do.call("cbind", lapply(levels(metadata_curate$ORIG_TYPE), function(x) apply(filter_normed_salmon_MAP_counts[, grep(x, colnames(filter_normed_salmon_MAP_counts))], 1, mean)))

colnames(big_averages_frame) <- paste0(levels(metadata_curate$ORIG_TYPE), "_AVG")

big_curr_output <- data.frame(filter_normed_salmon_MAP_counts, big_averages_frame, big_sig_hits)

curr_name <- paste0(outfix, "_salmon_counts_FILTERED_NORMALIZED_PERSAMPLE_CROSS_BATCH_MATRIX.csv")

fwrite(big_curr_output, curr_name, row.names=F, sep=",")
writeLines(as.character(big_curr_output$gene), gsub(".csv", "_geneset.txt", curr_name))

##padj significant output.
big_sig_output <- big_curr_output[big_curr_output$padj < 0.05,]
big_sig_output <- big_sig_output[!is.na(big_sig_output$padj),]
fwrite(big_sig_output, gsub(".csv", "_significant.csv", curr_name), row.names=F, sep=",")

print("finished!")
print(x)
setwd(start_super) 
}


}


