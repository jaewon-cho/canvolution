#!/bin/bash
#usage: ./1_variant_calling_single.sh fastq_file_name number_of_core
#$1="fastq_file_name" (ex: pat1.fastq --> pat1)
#$2="number of core"

#output usage
#VCF: $cwd/$1/$1.filtered.vcf.gz
#BAM file: $cwd/$1/$1.star.Aligned.sortedByCoord.out.bam


module load singularity


cwd=$(pwd)
str1=".fastq"
file1="$1$str1"

singularity exec -e -B ${CTAT_GENOME_LIB}:/data/hemberg/jaewon/ctat/new5/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ ~/ctat_mutations.v3.2.0.simg ~/ctat-mutations-CTAT-Mutations-v3.2.0/ctat_mutations --left $file1 --genome_lib_dir ~/sfa_jaewon/ctat/new5/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ -O $1 --cpu $2 --sample_id $1 --boosting_method=none



