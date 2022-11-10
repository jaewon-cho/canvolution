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
#' gmt <- getGmt("lung_gmt")
#' gmt_list$specific_fitness <- gmt
#' gmt <- getGmt("lung_gmt")
#' gmt_list$ctg <- gmt
#' mut_gmt(sample_list, gmt_list, "data/clone_mutation_info/", "_clone_mut_list")
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
			for (m in c(1:length(gmt_list))){
				gmt <- gmt_list[[m]]
				tmp_gmt_res <- path_gmt(gmt, gn)
				gmt_res <- c(gmt_res, tmp_gmt_res)
			}
###
			res<-cbind(res, tmp_res)

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
#' gmt <- getGmt("lung_gmt")
#' gmt_list$specific_fitness <- gmt
#' gmt <- getGmt("lung_gmt")
#' gmt_list$ctg <- gmt
#' mut_gmt(seurat_object, gmt_list, "lung")
#'
#' @export
#cluster name: seurat_clusters
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
#'
#' @return jaccard coefficient between LR gene and counter celltypes for each clone in the sample
#'
#' @examples
#' mut_gmt_lr(sample_list, "data/cellchat/", "lung")
#' @export
mut_gmt_lr <- function(sample_list, file_prefix, name){

##
	for (j in c(1:length(sample_list))){
		pat_tmp <- sample_list[j]
		print(pat_tmp)
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
		write.table(res, paste(pat_tmp, "_mut_LR_jacc",sep = ""), quote = F, sep = "\t")

#
	}
##
} 




