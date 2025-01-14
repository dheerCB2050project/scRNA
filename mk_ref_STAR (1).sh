#!/bin/bash
# naiss2023-5-259
#SBATCH -A naiss2023-5-259
#SBATCH -p node
#SBATCH -t 0-1:00:00
#SBATCH -n 1
#SBATCH -J mk_ref_STAR.sh
#SBATCH -e mk_ref_STAR.err
#SBATCH -o mk_ref_STAR.o
#SBATCH --mail-user emilio.skarwan@scilifelab.se
#SBATCH --mail-type=ALL

## Script to setup the reference to be used with STARsolo

module load bioinfo-tools
module load star/2.7.11a


VERSION=109
ref_name=Macaca_fascicularis_6.0."$VERSION".index
gtf_file=Macaca_fascicularis.Macaca_fascicularis_6.0."$VERSION".gtf
fasta_file=Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa

cd /proj/rnaatlas/nobackup/private/EmilioTemp/macaque_sc/data/alignment

STAR --runMode genomeGenerate \
 --runThreadN 18 \
 --genomeDir $ref_name \
 --genomeFastaFiles $fasta_file \
 --sjdbGTFfile $gtf_file \
 --sjdbOverhang 99
