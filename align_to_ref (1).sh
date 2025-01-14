#!/bin/bash
# naiss2023-5-259
#SBATCH -A naiss2023-5-259
#SBATCH -p node
#SBATCH -t 0-6:00:00
#SBATCH -n 1
#SBATCH -J align_to_ref.sh
#SBATCH -e align_to_ref.err
#SBATCH -o align_to_ref.o
#SBATCH --mail-user emilio.skarwan@scilifelab.se
#SBATCH --mail-type=ALL


module load bioinfo-tools
module load star/2.7.11a

cd /proj/rnaatlas/private/macaque/macaque_sc

output_dir=/proj/rnaatlas/nobackup/private/EmilioTemp/macaque_sc/data/raw_data
reference_directory=/proj/rnaatlas/nobackup/private/EmilioTemp/macaque_sc/data/alignment/Macaca_fascicularis_6.0.109.index


for organ in *; do
  if [ ! -d "$output_dir"/"$organ"/ ]; then
    echo ********************************
    echo "$output_dir"/"$organ"/
  mkdir "$output_dir"/"$organ"/
fi
  for sample in "$organ"/*; do
    if [[ "$sample" == *_Oligo ]]; then
      continue
    fi
    if [ ! -d "$output_dir"/"$sample"/ ]; then
      echo "$output_dir"/"$sample"/

      for read in "$sample"/*.fq.gz; do
      if [[ "$read" == *1.fq.gz ]]; then
        #echo "Storing barcode reads as read_1: $read"
        read_1=$read
      elif [[ "$read" == *2.fq.gz ]]; then
        #echo "Storing biological reads as read_2: $read"
        read_2=$read
      fi
    done
    mkdir "$output_dir"/"$sample"
    ## Perform alignment for the sample
    STAR --genomeDir "$reference_directory" \
                   --soloType CB_UMI_Simple \
                   --readFilesCommand zcat \
                   --readFilesIn "$read_2" "$read_1" \
                   --soloCBwhitelist /proj/rnaatlas/nobackup/private/EmilioTemp/macaque_sc/data/alignment/BC_whitelist.txt \
                   --soloCBstart 1 \
                   --soloCBlen 20 \
                   --soloUMIstart 21 \
                   --soloUMIlen 10 \
                   --outFileNamePrefix "$output_dir"/"$sample"/ \
             --limitOutSJcollapsed 2000000 \
             --runThreadN 18 \
             --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
             --soloCellFilter EmptyDrops_CR \
             --soloUMIfiltering MultiGeneUMI_CR \
             --soloUMIdedup 1MM_CR \
             --soloFeatures Gene Velocyto \
             --soloMultiMappers EM \
             --outSAMtype None  # just as a test, then decide which one to take

    fi
  done
done



#

#for organ in *;do
#  if [[ "$organ" != Heart ]]; then #
#

#    mkdir "$output_dir"/"$organ"/#

#    for sample in "$organ"/*; do#

#    # Skip olig samples, that have no biological reads
#    if [[ "$sample" == *_Oligo ]]; then
#      continue
#    fi#

#     for read in "$sample"/*.fq.gz; do
#     if [[ "$read" == *1.fq.gz ]]; then
#       #echo "Storing barcode reads as read_1: $read"
#       read_1=$read
#     elif [[ "$read" == *2.fq.gz ]]; then
#       #echo "Storing biological reads as read_2: $read"
#       read_2=$read
#     fi#

#   done#

#   mkdir "$output_dir"/"$sample"
#   ## Perform alignment for the sample
#   STAR --genomeDir "$reference_directory" \
#                  --soloType CB_UMI_Simple \
#                  --readFilesCommand zcat \
#                  --readFilesIn "$read_2" "$read_1" \
#                  --soloCBwhitelist /proj/rnaatlas/nobackup/private/EmilioTemp/macaque_sc/data/alignment/BC_whitelist.txt \
#                  --soloCBstart 1 \
#                  --soloCBlen 20 \
#                  --soloUMIstart 21 \
#                  --soloUMIlen 10 \
#                  --outFileNamePrefix "$output_dir"/"$sample"/ \
#            --limitOutSJcollapsed 2000000 \
#            --runThreadN 18 \
#            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
#            --soloCellFilter EmptyDrops_CR \
#            --soloUMIfiltering MultiGeneUMI_CR \
#            --soloUMIdedup 1MM_CR \
#            --soloFeatures Gene Velocyto \
#            --soloMultiMappers EM \
#            --outSAMtype None  # just as a test, then decide which one to take#

#   # eventually test out --soloMultiMappers PropUnique or Rescue
#   # –outFilterScoreMin 30.#

#  done#

#fi
###·####
#done



#for organ in *;do
#  if [[ "$organ" == Heart ]]; then
#    mkdir "$output_dir"/"$organ"/#

#    for sample in "$organ"/*; do#

#       for read in "$sample"/*.fq.gz; do
#       if [[ "$read" == *1.fq.gz ]]; then
#         #echo "Storing barcode reads as read_1: $read"
#         read_1=$read
#       elif [[ "$read" == *2.fq.gz ]]; then
#         #echo "Storing biological reads as read_2: $read"
#         read_2=$read
#       fi#

#     done
#     mkdir "$output_dir"/"$organ"_PropUnique/"$sample"
#     ## Perform alignment for the sample
#     STAR --genomeDir "$reference_directory" \
#                    --soloType CB_UMI_Simple \
#                    --readFilesCommand zcat \
#                    --readFilesIn "$read_2" "$read_1" \
#                    --soloCBwhitelist /proj/rnaatlas/nobackup/private/EmilioTemp/macaque_sc/data/alignment/BC_whitelist.txt \
#                    --soloCBstart 1 \
#                    --soloCBlen 20 \
#                    --soloUMIstart 21 \
#                    --soloUMIlen 10 \
#                    --outFileNamePrefix "$output_dir"/"$organ"_PropUnique/"$sample"/ \
#              --limitOutSJcollapsed 2000000 \
#              --runThreadN 18 \
#              --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
#              --soloCellFilter EmptyDrops_CR \
#              --soloUMIfiltering MultiGeneUMI_CR \
#              --soloUMIdedup 1MM_CR \
#              --soloFeatures Gene GeneFull Velocyto \
#              --soloMultiMappers PropUnique \
#              --outSAMtype None  # just as a test, then decide which one to take#

#     # eventually test out --soloMultiMappers PropUnique or Rescue
#     # –outFilterScoreMin 30.#

#    done
#  fi#

#done


#
#

# STAR --genomeDir "$reference_directory" \
#                --soloType CB_UMI_Simple \
#                --readFilesCommand zcat \
#                --readFilesIn "$read_2" "$read_1" \
#                --soloCBwhitelist /proj/rnaatlas/nobackup/private/EmilioTemp/whitelist_10x.txt \
#                --soloCBstart 1 \
#                --soloCBlen 16 \
#                --soloUMIstart 17 \
#                --soloUMIlen 10 \
#                --outFileNamePrefix "$output_dir" \
#          --limitOutSJcollapsed 2000000 \
#          --runThreadN 20 \
#          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
#          --soloCellFilter EmptyDrops_CR \
#          --soloUMIfiltering MultiGeneUMI_CR \
#          --soloUMIdedup 1MM_CR \
#          --soloFeatures Gene GeneFull Velocyto # just as a test, then decide which one to take
# # eventually test out --soloMultiMappers PropUnique or Rescue
# # –outFilterScoreMin 30.#
#
#

# –soloFeatures Gene Velocyto
# –soloBarcodeReadLength 0;
# –soloType CB_UMI_Simple;
# –soloCellFilter Empty_Drops_CR %s 0.99 10 45000 90000 500 0.01 20000 0.01 10000;
# –soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts;
# –soloUMIfiltering MultiGeneUMI_CR;
# –soloUMIdedup 1MM_CR
# –clipAdapterType CellRanger4;
# –outFilterScoreMin 30.#
#
#

# STAR --genomeDir $reference_directory \
#                --soloType CB_UMI_Simple \
#                --readFilesCommand zcat \
#                --readFilesIn "$2" "$1" \
#                --soloCBwhitelist /proj/rnaatlas/nobackup/private/EmilioTemp/macaque_sc/data/alignment/BC_whitelist.txt \
#                --soloCBstart 1 \
#                --soloCBlen 16 \
#                --soloUMIstart 17 \
#                --soloUMIlen 10 \
#                --outFileNamePrefix $homedir"/single_nuclei_10x_star/"$acc"/" \
#          --limitOutSJcollapsed 2000000 \
#          --runThreadN 20 \
#          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
#          --soloUMIfiltering MultiGeneUMI_CR \
#          --soloUMIdedup 1MM_CR
