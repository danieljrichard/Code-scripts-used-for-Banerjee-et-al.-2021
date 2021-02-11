###########
##Running Salmon Pseudo-alignment for gene expression quantification
##Banerjee et al 2021
##
##

#module load R/3.5.1-fasrc01


human_reference <- "/n/capellini_lab/users/drichard/FRAG_RNA_SEQ/hg19_reference_salmon014/hg19_transcriptome_salmon014"

salmon_quant <- function(fastq_list, reference, cores = 12, boots=100, special_fix="", local = TRUE, outfolder) {
	
	salmon_command <- "/n/capellini_lab/users/drichard/FRAG_RNA_SEQ/salmon_014/salmon-latest_linux_x86_64/bin/./salmon"
	
	fileset <- readLines(fastq_list)
	fileset <- fileset[grep(".fastq", fileset)]
		
	fileset <- fileset[!grepl("_map", fileset)]
	
	##Using file paths.
	##Lining up forward/reverse reads.
	orig_files <- sapply(fileset, function(x) unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))])
	
	fores <- fileset[grep("_1", fileset)]
	hinds <- fileset[grep("_2", fileset)]
	orig_files <- unique(unlist(sapply(orig_files, function(x) unlist(strsplit(x, "_[1-2]"))[1])))

if(length(fores) == 0) {
	
	fores <-  fileset[grep(".R1.", fileset, fixed=T)]
	hinds <- fileset[grep(".R2.", fileset, fixed=T)]
	orig_files <- unique(unlist(sapply(orig_files, function(x) unlist(strsplit(x, ".R[1-2]"))[1])))

	}
	
	if(length(fores) == 0) {
		
	fores <-  fileset[grep("_R1", fileset, fixed=T)]
	hinds <- fileset[grep("_R2", fileset, fixed=T)]
	orig_files <- sapply(fileset, function(x) unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))])
	orig_files <- unique(unlist(sapply(orig_files, function(x) unlist(strsplit(x, "_R[1-2]"))[1])))

	}
	
	
	for(orig in orig_files) {
		print("sample")
		print(orig)
		
		curr_fore <- fores[grep(orig, fores)]
		
		if(length(curr_fore) != 1) {
			print("Error, expected R1.fastq and R2.fastq files")
			browser()
			
		}
		
		curr_hind <- hinds[grep(orig, hinds)]
		print(curr_fore)
		print(curr_hind)
		
		if(length(curr_hind) == 0) {
			print("single")
			system(paste0("echo '", curr_fore, "' >> singles_set.txt"))
			next()
		}
		
		sbatch_lines <- list()
		sbatch_lines [1] <- "#!/bin/bash -l"
		sbatch_lines[2] <- "#SBATCH -p shared"
		sbatch_lines[3] <- paste0("#SBATCH -n ", cores)
		sbatch_lines[4] <- "#SBATCH -N 1"
		sbatch_lines[5] <- "#SBATCH -t 0-06:00"
		sbatch_lines[6] <- "#SBATCH --mem=20G"
		sbatch_lines[7] <- paste0("#SBATCH -J ", orig, "_salmon_map")
		sbatch_lines[8] <- paste0("#SBATCH -o ", orig, "_salmon_map.out")
		sbatch_lines[9] <- paste0("#SBATCH -e ", orig, "_salmon_map.err")
		
if(local == FALSE) {
	curr_out <- unlist(strsplit(curr_fore, "R1"))[1]
	curr_out <- paste0(curr_out, "_salmon_results")
}else{
	if(special_fix != "") {
		
	curr_out <- paste0(outfolder, "/", orig, "_salmon_results", "_", special_fix)
}else{
	curr_out <- paste0(outfolder, "/", orig, "_salmon_results")

}
}
	
	sbatch_lines[10] <- paste0(salmon_command, " quant -i ", reference, " -l A -1 <(zcat ", curr_fore, ")", " -2 <(zcat ", curr_hind, ") -o ", curr_out, " --numBootstraps ", boots, " -p ", cores, " --gcBias --validateMappings")
			
		sbatch_lines[11] <- "echo 'finished'"
		if(special_fix == "") {
		outfile <- file(paste0(orig, "_salmon_quant.sbatch"))
	}else{
		outfile <- file(paste0(orig, "_salmon_quant_", special_fix, ".sbatch"))
	}
		writeLines(unlist(sbatch_lines), outfile)
		close(outfile)
	

}

	
}


salmon_quant("gz_set.txt", human_reference, special_fix = "SARS2_CALU3", outfolder = "/n/holyscratch01/capellini_lab/drichard/SECRET/RESULTS")
