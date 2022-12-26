#' path_gmt
#'
#' Evaulate jaccard coefficient between input gene and geneset from a given gmt (pathway)
#'
#' @param gmt gmt object from GSEABase package
#' @param gn input genelist
#'
#' @return jaccard coefficient between input gene and geneset from a given gmt (pathway)
#'
#' @examples
#' tmp_gmt_res <- path_gmt(gmt, gn)
#' @export
path_gmt <- function(gmt, gn){
	tmp_jacc  <- c()
	for(i in c(1:length(gmt@.Data))){
		path_gn <- gmt@.Data[[i]]@geneIds
                path <- gmt@.Data[[i]]@setName
		jacc <- length(intersect(gn, path_gn)) / length(union(gn, path_gn))
		tmp_jacc  <- c(tmp_jacc, jacc )
		names(tmp_jacc)[i] <- path
	}

	return (tmp_jacc)
}            

#' path_gmt_lr
#'
#' Evaulate jaccard coefficient between ligand-receptor gene and a given geneset
#'
#' @param res list object contains receptor (pre_name: target) or ligand (pre_name: sourcr) list of each cancer clone between each celltype
#' @param gn input gene
#' @param pre_name "target" (receptor gene from cancer) or "source" (ligand gene from cancer)
#'
#' @return jaccard coefficient for each interaction
#' 
#' @examples
#' load("data/cellchat/lung_target_gn")
#' #object name: res
#' tmp_res <- path_gmt_lr(res, gn, "target")
#'
#' load("data/cellchat/lung_source_gn")
#' #object name: res
#' tmp_res <- path_gmt_lr(res, gn, "source")
#' 
#' @export
path_gmt_lr<-function(res, gn, pre_name){
	tmp_jacc<- c()

	for(i in c(1:length(res))){

		counter<-res[[i]]
		counter_name <- names(res)[i]
		for (j in c(1:length(counter))){
			path_gn <-counter[[j]]
	                cluster <- names(counter)[j]
			jacc <- length(intersect(gn, path_gn)) / length(union(gn, path_gn))
			tmp_jacc  <- c(tmp_jacc, jacc )
			names(tmp_jacc)[length(tmp_jacc)] <- paste(pre_name, cluster,counter_name, sep = "_")
		}
	}
	return (tmp_jacc)
}            


#' mut_gmt
#'
#' measuring jaccard coefficient between mutational profile and multiple pathways for each clone in the sample 
#'
#' @param sample_list list of samples
#' @param gmt_list multiple pathways in GSEABase gmt object 
#' @param clone_info_prefix prefix of clone_info
#' @param clone_info_suffix suffix of clone_info
#'
#' @return jaccard coefficient between mutational profile and multiple pathways for each clone in the sample
#'
#' @examples
#' gmt_list <- list()
#' gmt <- getGmt("cancersea_gmt")
#' gmt_list$cancersea<- gmt
#' gmt <- getGmt("core_gmt")
#' gmt_list$core_fitness <- gmt
#' gmt <- getGmt("lung_cm_gmt")
#' gmt_list$cancermine <- gmt
#' mut_gmt(sample_list, sample_clone_mut_list, gmt_list, "data/clone_mutation_info/", "_clone_mut_list")
#' @export
mut_gmt <- function(sample_list, gmt_list, clone_info_prefix,clone_info_suffix){
	for (j in c(1:length(sample_list))){
		pat_tmp <- sample_list[j]
		print(pat_tmp)
		load(paste(clone_info_prefix,pat_tmp, clone_info_suffix, sep = ""))
		res <- c()
		clone_list <- names(sample_clone_mut_list)
		for (k in c(1:length(clone_list))){
			clone_name <- clone_list[k]
			gn <- sample_clone_mut_list[[k]]
			gmt_res <- c()
			for (m in c(1:length(gmt_list))){
				gmt <- gmt_list[[m]]
				tmp_gmt_res <- path_gmt(gmt, gn)
				gmt_res <- c(gmt_res, tmp_gmt_res)
			}
###
			res<-cbind(res, gmt_res)

			colnames(res)[k] <- clone_name
		}
	write.table(res, paste(pat_tmp, "_mut_gmt_jacc",sep = ""), quote = F, sep = "\t")
#	write.csv(res, paste(pat_tmp, "_mut_gmt_jacc.csv",sep = ""), quote = F)
#
	}
##
} 



#' path_module
#'
#' measure signature score for each clone by gene expression in the sample with a given pathway
#'
#' @param gmt gmt object from GSEABase package
#' @param d seurat object
#' @param total_gn total gene list in the data
#'
#' @return signature score for each clone by gene expression in the sample
#'
#' @examples
#' gmt <- gmt_list[[z]]
#' tmp_res <- path_module(gmt, d3, total_gn)
#' @export
path_module<-function(gmt, d, total_gn){
	res <- c()
	for (i in c(1:length(gmt@.Data))){
                tmp_path_gn <- gmt@.Data[[i]]@geneIds
                path <- gmt@.Data[[i]]@setName
		print(path)
		path_gn<-intersect(tmp_path_gn, total_gn)
		d1<-AddModuleScore(d, features = path_gn)
		#d1<-AddModuleScore(d, features = path_gn, ctrl = length(path_gn))
		pre_size <- dim(d@meta.data)[2]
		pre_size <- pre_size + 1
		d2<-apply(d1@meta.data[,pre_size:dim(d1@meta.data)[2]], 1, mean)

		res<-rbind(res, d2)
		rownames(res)[dim(res)[1]] <- path
	}
	return(res)
}


#' addmodulescore_multiple_pathway
#'
#' measure signature score for each clone in the sample with multiple pathways
#'
#' @param d3 seurat_object
#' @param gmt_list gmt list
#' @param name title of the data
#'
#' @return signature score of pathways by gene expression for each clone in the sample (multiple pathways)
#'
#' @examples
#' gmt_list <- list()
#' gmt <- getGmt("cancersea_gmt")
#' gmt_list$cancersea<- gmt
#' gmt <- getGmt("core_gmt")
#' gmt_list$core_fitness <- gmt
#' gmt <- getGmt("lung_cm_gmt")
#' gmt_list$cancermine <- gmt
#' addmodulescore_multiple_pathway(seurat_object, gmt_list, "lung")
#' #cluster name: seurat_clusters
#' @export
addmodulescore_multiple_pathway <- function(d3, gmt_list, name){
	total_gn <- rownames(d3)

	res <- c()
	for(z in c(1:length(gmt_list))){
		gmt <- gmt_list[[z]]
		gmt_name <- names(gmt_list)[z]

		tmp_res <- path_module(gmt, d3, total_gn)
		res <- rbind(res, tmp_res)
	}
	write.table(res, paste(name, "_cell_addmodule",sep = ""), quote = F, sep = "\t")

}

#' mut_gmt_lr
#'
#' jaccard coefficient between LR gene and counter celltypes for each clone in the sample
#'
#' @param sample_list list of samples
#' @param file_prefix prefix of lr file
#' @param name title of the data
#' @param clone_info_prefix prefix of clone_info
#' @param clone_info_suffix suffix of clone_info
#'
#' @return jaccard coefficient between LR gene and counter celltypes for each clone in the sample
#'
#' @examples
#' mut_gmt_lr(sample_list, "data/cellchat/", "lung", "data/clone_mutation_info/", "_clone_mut_list")
#' @export
mut_gmt_lr <- function(sample_list, file_prefix, name, clone_info_prefix,clone_info_suffix){

##
	for (j in c(1:length(sample_list))){
		pat_tmp <- sample_list[j]
		print(pat_tmp)
                load(paste(clone_info_prefix,pat_tmp, clone_info_suffix, sep = ""))
		res2 <- c()
		clone_list <- names(sample_clone_mut_list)
		for (k in c(1:length(clone_list))){
			clone_name <- clone_list[k]
			gn <- sample_clone_mut_list[[k]]

			res1 <- c()
			load(paste(file_prefix,name, "_target_gn", sep = ""))
			tmp_res <- path_gmt_lr(res, gn, "target")
			res1<-c(res1, tmp_res)

			load(paste(file_prefix,name, "_source_gn", sep = ""))
			tmp_res <- path_gmt_lr(res, gn, "source")
			res1 <- c(res1, tmp_res)

			res2<-cbind(res2, res1)

			colnames(res2)[dim(res2)[2]] <- clone_name
		}
		write.table(res2, paste(pat_tmp, "_mut_LR_jacc",sep = ""), quote = F, sep = "\t")

#
	}
##
} 


#' mut_gmt_module
#'
#' mean expression of mutated gene of a given pathway in each clone
#'
#' @param exp expression matrix
#' @param meta_info seurat_obj@meta.data
#' @param sample_info variable name for sample in the "meta_info"
#' @param gmt_list gmt list
#' @param clone_info_prefix prefix of clone_info
#' @param clone_info_suffix suffix of clone_info
#' @param name title of the data
#'
#' @return mean expression of pathways by (mutated) gene expression for each clone in the sample (multiple pathways)
#'
#' @examples
#' gmt_list <- list()
#' gmt <- getGmt("cancersea_gmt")
#' gmt_list$cancersea<- gmt
#' gmt <- getGmt("core_gmt")
#' gmt_list$core_fitness <- gmt
#' gmt <- getGmt("lung_cm_gmt")
#' gmt_list$cancermine <- gmt
#' load("your_seurat_object")
#' #object name: d3
#' sample_info <- "sample_name"
#' exp <- GetAssayData(d3)
#' meta_info <- d3@meta.data
#' name <- "lung"
#' clone_info_prefix <- "data/clone_mutation_info/"
#' clone_info_suffix <- "_clone_mut_list"
#' mut_gmt_module(exp, meta_info, sample_info, gmt_list, clone_info_prefix,clone_info_suffix, name)
#' @export
mut_gmt_module <- function(exp, meta_info, sample_info, gmt_list, clone_info_prefix,clone_info_suffix, name){
sample_list <- as.character(meta_info[,sample_info])
sample_list_uniq <- unique(sample_list)
clone_list <- as.character(meta_info[,"cancer_clone"])
	clone_list <- paste("clone", clone_list, sep = "_")

total_gn <- rownames(exp)
res <- c()
for (m in c(1:length(gmt_list))){
	gmt <- gmt_list[[m]]
	print(names(gmt_list)[m])
	for(i in c(1:length(gmt@.Data))){
		path_gn <- gmt@.Data[[i]]@geneIds
		path <- gmt@.Data[[i]]@setName
	
	tmp_res <- rep(0, dim(exp)[2])
		for (j in c(1:length(sample_list_uniq))){
			pat_tmp <- sample_list_uniq[j]
load(paste(clone_info_prefix,pat_tmp, clone_info_suffix, sep = ""))



clone_list2 <- names(sample_clone_mut_list)

for (k in c(1:length(clone_list2))){
clone_name <- clone_list2[k]
					tmp_clone_list <- paste("clone", clone_list, sep = "_")
clone_loc <- (which(sample_list == pat_tmp & clone_list == clone_name))
gn <- sample_clone_mut_list[[k]]
tmp_path_gn <- intersect(gn, path_gn)
path_gn1 <- intersect(tmp_path_gn, total_gn)

if (length(path_gn1) < 1){
next
}

tmp_exp<- exp[path_gn1, clone_loc, drop = F]
tmp_exp2 <- apply(tmp_exp, 2, mean)
tmp_res[clone_loc] <- tmp_exp2
}
}
res <- rbind(res, tmp_res)
rownames(res)[dim(res)[1]] <- path
}
	}
colnames(res) <- colnames(exp)
	write.table(res, paste(name, "_mut_gmt_addmodule",sep = ""), quote = F, sep = "\t")

}

#' lr_gmt
#'
#' make cellchat result into list format
#'
#' @param source list object of (ligand gene from cancer)
#' @param target list object of (receptor gene from cancer)
#' @param name title of the data
#'
#' @return list variable composed of geneset for each source/target of counter_celltype
#'
#' @examples
#' name <- "lung"
#' 
#' load(paste("../cellchat/",name, "_target_gn", sep = ""))
#' target<- res
#' load(paste("../cellchat/",name, "_source_gn", sep = ""))
#' source<- res
#' res <- lr_gmt(source, target, name)
#' save(res, file = paste(name, "_cellchat_gmt", sep = "")) 
#'
#' @export
lr_gmt <- function(source, target, name){
	
	res <- list()
	
	for (m in c(1:length(source))){
		cell <- names(source)[m]
		cell1 <- paste("source",cell, sep = "__")
		cell1 <- gsub("-", ".",cell1)
		gn <- c()
		for (k in c(1:length(source[[m]]))){
			gn <- c(gn, source[[m]][[k]])
		}
		gn <- unique(gn)
		res[[m]] <- gn
		names(res)[m] <- cell1	
		
	}
	for (m in c(1:length(target))){
				
		cell <- names(target)[m]
		cell1 <- paste("target",cell, sep = "__")
		cell1 <- gsub("-", ".",cell1)
		gn <- c()
		for (k in c(1:length(source[[m]]))){
			gn <- c(gn, source[[m]][[k]])
		}
		gn <- unique(gn)
		res[[length(res)+1]] <- gn
		names(res)[length(res)] <- cell1	

	}	
	return (res)
	
} 

#' lr_module
#'
#' addmodulescore of each source/target gene in each cancer cell
#'
#' @param res object from lr_gmt
#' @param d seurat object
#' @param total_gn total gene list in the seurat object
#'
#' @return signature score of LR by gene expression for each cancer cell in the sample 
#'
#' @examples
#' name <- "lung"
#' print(name)
#' load(paste(name, "_cellchat_gmt", sep = ""))
#' #cellchat gmt object name: res (from lr_gmt)
#' load(paste(name, "_cancer_seurat_mut", sep = ""))
#' #seurat object name: d3
#' total_gn <- rownames(d3)
#' res2 <- lr_module(res, d3,total_gn)
#' write.table(res2, paste(name,"_cell_addmodule", sep = ""), quote = F, sep = "\t")
#'
#' @export
lr_module<-function(res, d, total_gn){
	res1 <- c()
	for (i in c(1:length(res))){
                tmp_path_gn <- res[[i]]
                path <- names(res)[i]
		print(path)
		path_gn<-intersect(tmp_path_gn, total_gn)
		if (length(path_gn) < 3){
			next
		}
		d1<-AddModuleScore(d, features = path_gn)
		#d1<-AddModuleScore(d, features = path_gn, ctrl = length(path_gn))
		pre_size <- dim(d@meta.data)[2]
		pre_size <- pre_size + 1
		d2<-apply(d1@meta.data[,pre_size:dim(d1@meta.data)[2]], 1, mean)

		res1<-rbind(res1, d2)
		rownames(res1)[dim(res1)[1]] <- path
	}
	return(res1)
}

