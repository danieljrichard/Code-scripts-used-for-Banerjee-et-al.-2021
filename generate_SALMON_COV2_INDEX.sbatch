#!/bin/bash
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 10                  # Number of cores
#SBATCH -t 180                 # Runtime in minutes (0~10080)
#SBATCH -p general             # Partition
#SBATCH --mem=50000            # Total memory (varies across nodes)
#SBATCH -o salmon_ind_%j.out         # Standard out goes to this file
#SBATCH -e salmon_ind_%j.err         # Standard err goes to this file
#SBATCH --mail-type=END        # Email
#SBATCH --mail-user=YOUR_EMAIL

#module load salmon

/n/capellini_lab/users/drichard/FRAG_RNA_SEQ/salmon_014/salmon-latest_linux_x86_64/bin/./salmon index -t GCF_009858895.2_ASM985889v3_cds_from_genomic.fna -i COVID_SALMON_INDEX -k 31
