#adjusting K value for clonotyping and merging into seurat object

library(DENDRO)
source("inst/dendro_function.R")
pat_data <- read.table("data/pat_optk", sep = "\t")
#load("data/lung_seurat")
#not provided
#seurat object; name: d2

d2$cancer_clone <- NA
b1<-names(d2$orig.ident)

for (j in c(1: dim(pat_data)[1])){
	print(rownames(pat_data)[j])
	tmp_label <- dendro_process2(rownames(pat_data)[j],pat_data$optk[j])
	a1<-names(tmp_label)
	a2<-match(a1,b1)
	d2$cancer_clone[a2] <- tmp_label
}


save(d2, file = "lung_seurat_mut")

