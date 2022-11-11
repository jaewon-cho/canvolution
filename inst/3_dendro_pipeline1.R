#!/PHShome/jc1062/.conda/envs/R4/bin/Rscript

#Program: dendro_pipeline.R
#Author: Jae-Won,Cho
#Date: 20220322(Tue)
#Description: pipeline for one patient in  dendro
#Usage: ./dendro_pipeline.R [merged.vcf.name]
#ex: ./dendro_pipeline.R patient_merged
#vcf_name: AZ003.vcf.gz -> AZ003

library(DENDRO)
library(data.table)
library(tidyr)
library(dplyr)

source("inst/dendro_function.R")
name <- commandArgs(TRUE)

dendro_process1(name)

