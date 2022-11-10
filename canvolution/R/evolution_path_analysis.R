#' clone_path_list
#'
#' split each evolution path from total evolution path (igraph) object
#'
#' @param sample_clone_path total evolution path (igraph) object
#'
#' @return each evolution path in a list object
#'
#' @examples
#' clone_path <- read.table("data/robustclone/TH155_clone_path")
#' clone_path <- as.matrix(clone_path)
#' clone_path[1,] <- as.character(clone_path[1,])
#' clone_path[2,] <- as.character(clone_path[2,])
#' clone_path_graph <-graph_from_edgelist(clone_path)
#' sample_clone_path <- subgraph(clone_path_graph, sample_clone_list)
#' #sample_clone_list: clone list in that sample
#' total_path <- clone_path_list(sample_clone_path)
#'
#' @export
clone_path_list <- function(sample_clone_path){
leaves= names(which(degree(sample_clone_path, v = V(sample_clone_path), mode = "out")==0, useNames = T))
root= names(which(degree(sample_clone_path, v = V(sample_clone_path), mode = "in")==0, useNames = T))

total_path <- list()
for(j in c(1: length(root))){
tmp_root <- root[j]
tmp_path <- all_simple_paths(sample_clone_path, from = tmp_root, to = leaves)
total_path <- append(total_path, tmp_path)
}

return (total_path)
}

#' evolution_mut_corr_sample_wise
#'
#' Correlation analysis between each evolution path and mutational features
#'
#' @param meta_info seurat_obj@meta.data
#' @param name title for your data
#' @param pat_info variable name for patient in the "meta_info"
#' @param sample_info variable name for sample in the "meta_info"
#' @param clone_path_file_prefix prefix of clone path file
#' @param feature_file_prefix prefix of feature file
#' @param output_file_prefix prefix of outputfile
#' @param cell_cnt_thr threshold for excluding clones less than the given cell count
#'
#' @return Correlation analysis between each evolution path and mutational features
#'
#' @examples
#' #meta_info <- seurat_obj@meta.data
#' meta_info <- read.table("data/meta_info_lung", sep = "\t", check.names = F)
#' name <- "lung"
#' pat_info <- "pat_collapsed"
#' sample_info <- "sample_name"
#' clone_path_file_prefix <- "data/robustclone/"
#' feature_file_prefix <- "data/feature/"
#' output_file_prefix <- "./"
#' evolution_mut_corr_samplewise(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, output_file_prefix)
#'
#' @export
evolution_mut_corr_samplewise <- function(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, output_file_prefix, cell_cnt_thr=10){
	sample_list <- as.character(meta_info[,sample_info])
	sample_list_uniq <- unique(sample_list)

	pat_list <- as.character(meta_info[,pat_info])
	sample_pat <- table(sample_list, pat_list)
	pat_list_uniq <- colnames(sample_pat)

	sample_pat_map <- apply(sample_pat, 1, function(x) which(x!=0))

	clone_list <- as.character(meta_info[,"cancer_clone"])
	sample_clone <- table(sample_list, clone_list)
	clone_list_uniq <- colnames(sample_clone)

	sample_clone_map <- apply(sample_clone, 1, function(x) which(x!=0))


	for(j in c(1:length(sample_pat_map))){
		sample <- names(sample_pat_map)[j]
		print(sample)
		sample_size <- sum(sample_clone[j,])
			if (sample_size < cell_cnt_thr){
			next
		}


		pat <- pat_list_uniq[sample_pat_map[j]]
		sample_clone_list <- as.character(colnames(sample_clone)[sample_clone_map[sample][[1]]])

		clone_path <- read.table(paste(clone_path_file_prefix, pat, "_clone_path", sep =""))
		clone_path <- as.matrix(clone_path)
		clone_path[1,] <- as.character(clone_path[1,])
		clone_path[2,] <- as.character(clone_path[2,])

		clone_path_graph <-graph_from_edgelist(clone_path)

		sample_clone_path <- subgraph(clone_path_graph, sample_clone_list)
			if(length(E(sample_clone_path))==0){
			next
		}

		total_path <- clone_path_list(sample_clone_path)
		

		if (length(total_path) <1){

			next
		}

		
		a2 <- read.delim(paste(feature_file_prefix, pat, "_mut_gmt_jacc", sep = ""),sep = "\t", check.names = F)
		#a2<-na.omit(a)

		total_path_order <- c()
		total_path_scc <- c()
		pp <- 1:length(total_path)
		for (k in c(1:length(total_path))){
			path_order <- total_path[[k]]
			path_order1<- as.character(names(path_order))

			tmp <- paste(path_order1, collapse = "_")
			path_order_name <- paste(sample, tmp, sep = "__")
			path_order <- paste(pat, path_order1, sep ="_")


			a3 <- a2[,path_order]

			single_path_scc <- c()
			for(m in c(1: dim(a3)[1])){
				spearman <- cor.test(1:dim(a3)[2], as.numeric(a3[m,]), method = "spearman")
				val <- spearman$estimate

				single_path_scc <- c(single_path_scc, val)
				names(single_path_scc)[m] <- rownames(a3)[m]

			}
			total_path_scc  <- cbind(total_path_scc, single_path_scc)
			colnames(total_path_scc)[k] <- path_order_name

		}



		write.table(total_path_scc, paste(output_file_prefix, sample, "_clone_evolution_mut_scc_sample", sep = ""), quote = F, sep = "\t") 

	}
}


#' evolution_signature_cellwise_samplewise
#'
#' Correlation analysis between each evolution path and transcriptional features
#'
#' @param meta_info seurat_obj@meta.data
#' @param name title for your data
#' @param pat_info variable name for patient in the "meta_info"
#' @param sample_info variable name for sample in the "meta_info"
#' @param clone_path_file_prefix prefix of clone path file
#' @param feature_file_prefix prefix of feature file
#' @param output_file_prefix prefix of outputfile
#' @param cell_cnt_thr threshold for excluding clones less than the given cell count
#'
#' @return Correlation analysis between each evolution path and mutational features
#'
#' @examples
#' #meta_info <- seurat_obj@meta.data
#' meta_info <- read.table("data/meta_info_lung", sep = "\t", check.names = F)
#' name <- "lung"
#' pat_info <- "pat_collapsed"
#' sample_info <- "sample_name"
#' clone_path_file_prefix <- "data/robustclone/"
#' feature_file_prefix <- "data/feature/"
#' output_file_prefix <- "./"
#' evolution_signature_cellwise_samplewise(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, output_file_prefix)
#'
#' @export
evolution_signature_cellwise_samplewise <- function(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, output_file_prefix, cell_cnt_thr=10){
	sample_list <- as.character(meta_info[,sample_info])
	sample_list_uniq <- unique(sample_list)

	pat_list <- as.character(meta_info[,pat_info])
	sample_pat <- table(sample_list, pat_list)
	pat_list_uniq <- colnames(sample_pat)

	sample_pat_map <- apply(sample_pat, 1, function(x) which(x!=0))

	clone_list <- as.character(meta_info[,"cancer_clone"])
	sample_clone <- table(sample_list, clone_list)
	clone_list_uniq <- colnames(sample_clone)

	sample_clone_map <- apply(sample_clone, 1, function(x) which(x!=0))


	for(j in c(1:length(sample_pat_map))){
		sample <- names(sample_pat_map)[j]
		print(sample)
		sample_size <- sum(sample_clone[j,])
			if (sample_size < cell_cnt_thr){
			next
		}


		pat <- pat_list_uniq[sample_pat_map[j]]
		sample_clone_list <- as.character(colnames(sample_clone)[sample_clone_map[sample][[1]]])

		clone_path <- read.table(paste(clone_path_file_prefix, pat, "_clone_path", sep =""))
		clone_path <- as.matrix(clone_path)
		clone_path[1,] <- as.character(clone_path[1,])
		clone_path[2,] <- as.character(clone_path[2,])

		clone_path_graph <-graph_from_edgelist(clone_path)

		sample_clone_path <- subgraph(clone_path_graph, sample_clone_list)
			if(length(E(sample_clone_path))==0){
			next
		}

		total_path <- clone_path_list(sample_clone_path)
		

		if (length(total_path) <1){

			next
		}

		
		
		a2 <- read.delim(paste(feature_file_prefix, name, "_cell_addmodule", sep = ""),sep = "\t", check.names = F)
		

		total_path_order <- c()
		total_path_scc <- c()

		pp <- 1:length(total_path)
		for (k in c(1:length(total_path))){
			path_order <- total_path[[k]]
			path_order1<- as.character(names(path_order))

			tmp <- paste(path_order1, collapse = "_")
			path_order_name <- paste(sample, tmp, sep = "__")
			path_order <- paste(pat, path_order1, sep ="_")

			clone_a3 <- c()
			a3 <- c()
			for (z in c(1:length(path_order1))){
				tmp_path <- path_order1[z]
				clone_loc <- which(sample_list == sample & clone_list == tmp_path)
				tmp_a3 <- a2[,clone_loc]
				if (length(clone_loc) == 1){
					tmp_clone_a3 <- path_order[z]

				} else{

					tmp_clone_a3<-rep(path_order[z], dim(tmp_a3)[2])
				}
				clone_a3 <- c(clone_a3, tmp_clone_a3)
				a3<-cbind2(a3, as.matrix(tmp_a3))
			}
			rownames(a3) <- rownames(a2)

			

			single_path_scc <- c()
			
			for(m in c(1: dim(a3)[1])){

				b <- cbind(a3[m,],clone_a3)
				colnames(b) <- c("value", "clone")
				b<-data.frame(b)
				b$value <- as.numeric(b$value)
				b$clone <- factor(b$clone, levels=path_order, ordered=TRUE)
				bbb<- as.numeric(b$clone)
				spearman <- cor.test(bbb, b$value, method = "spearman")
				val <- spearman$estimate

			

				single_path_scc <- c(single_path_scc, val)
				names(single_path_scc)[m] <- rownames(a3)[m]
			

			}
			total_path_scc  <- cbind(total_path_scc, single_path_scc)
			colnames(total_path_scc)[k] <- path_order_name
			

		}



		write.table(total_path_scc, paste(output_file_prefix, name, "_", sample, "_clone_evolution_signature_scc_sample", sep = ""), quote = F, sep = "\t") 
		

	}
}


#' evolution_mut_signature_cellwise_corr_samplewise
#'
#' Correlation analysis between the same mutational features and transcriptional features in each evolution path
#'
#' @param meta_info seurat_obj@meta.data
#' @param name title for your data
#' @param pat_info variable name for patient in the "meta_info"
#' @param sample_info variable name for sample in the "meta_info"
#' @param clone_path_file_prefix prefix of clone path file
#' @param feature_file_prefix prefix of feature file
#' @param output_file_prefix prefix of outputfile
#' @param cell_cnt_thr threshold for excluding clones less than the given cell count
#'
#' @return Correlation analysis between the same mutational features and transcriptional features in each evolution path
#'
#' @examples
#' #meta_info <- seurat_obj@meta.data
#' meta_info <- read.table("data/meta_info_lung", sep = "\t", check.names = F)
#' name <- "lung"
#' pat_info <- "pat_collapsed"
#' sample_info <- "sample_name"
#' clone_path_file_prefix <- "data/robustclone/"
#' feature_file_prefix <- "data/feature/"
#' output_file_prefix <- "./"
#' evolution_mut_signature_cellwise_corr_samplewise(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, output_file_prefix)
#'
#' @export
evolution_mut_signature_cellwise_corr_samplewise <- function(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, output_file_prefix, cell_cnt_thr=10){
	sample_list <- as.character(meta_info[,sample_info])
	sample_list_uniq <- unique(sample_list)

	pat_list <- as.character(meta_info[,pat_info])
	sample_pat <- table(sample_list, pat_list)
	pat_list_uniq <- colnames(sample_pat)

	sample_pat_map <- apply(sample_pat, 1, function(x) which(x!=0))

	clone_list <- as.character(meta_info[,"cancer_clone"])
	sample_clone <- table(sample_list, clone_list)
	clone_list_uniq <- colnames(sample_clone)

	sample_clone_map <- apply(sample_clone, 1, function(x) which(x!=0))


	for(j in c(1:length(sample_pat_map))){
		sample <- names(sample_pat_map)[j]
		print(sample)
		sample_size <- sum(sample_clone[j,])
			if (sample_size < cell_cnt_thr){
			next
		}


		pat <- pat_list_uniq[sample_pat_map[j]]
		sample_clone_list <- as.character(colnames(sample_clone)[sample_clone_map[sample][[1]]])

		clone_path <- read.table(paste(clone_path_file_prefix, pat, "_clone_path", sep =""))
		clone_path <- as.matrix(clone_path)
		clone_path[1,] <- as.character(clone_path[1,])
		clone_path[2,] <- as.character(clone_path[2,])

		clone_path_graph <-graph_from_edgelist(clone_path)

		sample_clone_path <- subgraph(clone_path_graph, sample_clone_list)
			if(length(E(sample_clone_path))==0){
			next
		}

		total_path <- clone_path_list(sample_clone_path)
		

		if (length(total_path) <1){

			next
		}

		
		
		a2 <- read.delim(paste(feature_file_prefix, name, "_cell_addmodule", sep = ""),sep = "\t", check.names = F)
		a2_1 <- read.delim(paste(feature_file_prefix, pat, "_mut_kegg_jacc", sep = ""),sep = "\t", check.names = F)

		total_path_order <- c()
		total_path_scc <- c()

		pp <- 1:length(total_path)
		for (k in c(1:length(total_path))){
			path_order <- total_path[[k]]
			path_order1<- as.character(names(path_order))

			tmp <- paste(path_order1, collapse = "_")
			path_order_name <- paste(sample, tmp, sep = "__")
			path_order <- paste(pat, path_order1, sep ="_")

			a3_1 <- a2_1[,path_order]


			clone_a3 <- c()
			a3 <- c()
			for (z in c(1:length(path_order1))){
				tmp_path <- path_order1[z]
				clone_loc <- which(sample_list == sample & clone_list == tmp_path)
				tmp_a3 <- a2[,clone_loc]
				if (length(clone_loc) == 1){
					tmp_clone_a3 <- path_order[z]

				} else{

					tmp_clone_a3<-rep(path_order[z], dim(tmp_a3)[2])
				}
				clone_a3 <- c(clone_a3, tmp_clone_a3)
				a3<-cbind2(a3, as.matrix(tmp_a3))
			}
			rownames(a3) <- rownames(a2)
			tmp_total_path_result_scc<- c()
			tmp_total_path_result_pval<- c()

			for (w in c(1: dim(a3_1)[1])){

				mut_feature <- rownames(a3_1)[w]
	
				check_mut <- which(a3_1[w,] == a3_1[w,1])
				if (length(check_mut) == dim(a3_1)[2]){
					next
				}

				mut_feature_order <- order(as.numeric(a3_1[w,]))

				new_path_order <- path_order[mut_feature_order]


			

				single_path_coeff <- c()
			
				for(m in c(1: dim(a3)[1])){
#print(m)
					b <- cbind(a3[m,],clone_a3)
					colnames(b) <- c("value", "clone")
					b<-data.frame(b)
					b$value <- as.numeric(b$value)


					bbb<- melt(a3_1[w,], b$clone)
					b$value2 <- as.numeric(bbb)
					spearman <- cor.test(b$value2, b$value, method = "spearman")
					val <- spearman$estimate


			

					single_path_coeff <- c(single_path_coeff, val)
					names(single_path_coeff)[length(single_path_coeff)] <- paste(mut_feature, rownames(a3)[m], sep = "____")
			

				}
				tmp_total_path_result_scc  <- c(tmp_total_path_result_scc, single_path_coeff)
			}
			total_path_result_scc  <- cbind(total_path_result_scc, tmp_total_path_result_scc)
			colnames(total_path_result_scc)[k] <- path_order_name
		}



		write.table(total_path_scc, paste(output_file_prefix, name, "_", sample, "_clone_evolution_mut_signature_scc_sample", sep = ""), quote = F, sep = "\t") 
		

	}
}

#' evolution_clone_abundance_mut_corr_samplewise
#'
#' Correlation analysis between abundance of each clone and mutational features in a given evolution path
#'
#' @param meta_info seurat_obj@meta.data
#' @param name title for your data
#' @param pat_info variable name for patient in the "meta_info"
#' @param sample_info variable name for sample in the "meta_info"
#' @param clone_path_file_prefix prefix of clone path file
#' @param feature_file_prefix prefix of feature file
#' @param abundance_file_prefix prefix of abundance file
#' @param output_file_prefix prefix of outputfile
#' @param cell_cnt_thr threshold for excluding clones less than the given cell count
#'
#' @return Correlation analysis between abundance of each clone and mutational features in a given evolution path
#'
#' @examples
#' #meta_info <- seurat_obj@meta.data
#' meta_info <- read.table("data/meta_info_lung", sep = "\t", check.names = F)
#' name <- "lung"
#' pat_info <- "pat_collapsed"
#' sample_info <- "sample_name"
#' clone_path_file_prefix <- "data/robustclone/"
#' feature_file_prefix <- "data/feature/"
#' abundance_file_prefix <- "data/"
#' output_file_prefix <- "./"
#' evolution_clone_abundance_mut_corr_samplewise(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, output_file_prefix)
#'
#' @export
evolution_clone_abundance_mut_corr_samplewise <- function(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, abundance_file_prefix, output_file_prefix, cell_cnt_thr=10){
	sample_list <- as.character(meta_info[,sample_info])
	sample_list_uniq <- unique(sample_list)

	pat_list <- as.character(meta_info[,pat_info])
	sample_pat <- table(sample_list, pat_list)
	pat_list_uniq <- colnames(sample_pat)

	sample_pat_map <- apply(sample_pat, 1, function(x) which(x!=0))

	clone_list <- as.character(meta_info[,"cancer_clone"])
	sample_clone <- table(sample_list, clone_list)
	clone_list_uniq <- colnames(sample_clone)

	sample_clone_map <- apply(sample_clone, 1, function(x) which(x!=0))


	for(j in c(1:length(sample_pat_map))){
		sample <- names(sample_pat_map)[j]
		print(sample)
		sample_size <- sum(sample_clone[j,])
			if (sample_size < cell_cnt_thr){
			next
		}


		pat <- pat_list_uniq[sample_pat_map[j]]
		sample_clone_list <- as.character(colnames(sample_clone)[sample_clone_map[sample][[1]]])

		clone_path <- read.table(paste(clone_path_file_prefix, pat, "_clone_path", sep =""))
		clone_path <- as.matrix(clone_path)
		clone_path[1,] <- as.character(clone_path[1,])
		clone_path[2,] <- as.character(clone_path[2,])

		clone_path_graph <-graph_from_edgelist(clone_path)

		sample_clone_path <- subgraph(clone_path_graph, sample_clone_list)
			if(length(E(sample_clone_path))==0){
			next
		}

		total_path <- clone_path_list(sample_clone_path)
		

		if (length(total_path) <1){

			next
		}

		
		a <- read.csv(paste(abundance_file_prefix, name, "_clone_ratio.csv", sep = ""))
		a1<-a[, 2:dim(a)[2]]
		a2<-as.matrix(a1)
		rownames(a2)<-a$X
		a2<-data.frame(a2)

		
		a2_1 <- read.delim(paste(feature_file_prefix, pat, "_mut_kegg_jacc", sep = ""),sep = "\t", check.names = F)

		total_path_order <- c()
		total_path_scc <- c()

		pp <- 1:length(total_path)
		for (k in c(1:length(total_path))){
			path_order <- total_path[[k]]
			path_order1<- as.character(names(path_order))

			tmp <- paste(path_order1, collapse = "_")
			path_order_name <- paste(sample, tmp, sep = "__")
			
			path_order_abundance <- paste("clone", path_order1, sep ="_")
			path_order_abundance2  <- paste(pat, path_order1, sep ="_")
			a3_1 <- a2_1[,path_order_abundance2]
			a3<-a2[sample, path_order_abundance]


			single_path_scc <- c()
	

			for (w in c(1: dim(a3_1)[1])){
				mut_feature <- rownames(a3_1)[w]
	
				check_mut <- which(a3_1[w,] == a3_1[w,1])
				if (length(check_mut) == dim(a3_1)[2]){
					next
				}

				spearman <- cor.test(as.numeric(a3[1,]), as.numeric(a3_1[w,]), method = "spearman")
				val <- spearman$estimate
	
				single_path_scc <- c(single_path_scc, val)
				names(single_path_scc)[length(single_path_scc)] <- rownames(a3_1)[w]
	

			}

		total_path_scc  <- cbind(total_path_scc, single_path_scc)
		colnames(total_path_scc)[k] <- path_order_name

		}

	write.table(total_path_scc, paste(output_file_prefix, name, "_", sample, "_clone_abundance_mut_scc_sample", sep = ""), quote = F, sep = "\t") 
	}
}

#' evolution_clone_abundance_transcriptomic_corr_samplewise
#'
#' Correlation analysis between abundance of each clone and transcriptional features in a given evolution path
#'
#' @param meta_info seurat_obj@meta.data
#' @param name title for your data
#' @param pat_info variable name for patient in the "meta_info"
#' @param sample_info variable name for sample in the "meta_info"
#' @param clone_path_file_prefix prefix of clone path file
#' @param feature_file_prefix prefix of feature file
#' @param abundance_file_prefix prefix of abundance file
#' @param output_file_prefix prefix of outputfile
#' @param cell_cnt_thr threshold for excluding clones less than the given cell count
#'
#' @return Correlation analysis between abundance of each clone and transcriptional features in a given evolution path
#' @examples
#' #meta_info <- seurat_obj@meta.data
#' meta_info <- read.table("data/meta_info_lung", sep = "\t", check.names = F)
#' name <- "lung"
#' pat_info <- "pat_collapsed"
#' sample_info <- "sample_name"
#' clone_path_file_prefix <- "data/robustclone/"
#' feature_file_prefix <- "data/feature/"
#' abundance_file_prefix <- "data/"
#' output_file_prefix <- "./"
#' evolution_clone_abundance_transcriptomic_corr_samplewise(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, output_file_prefix)
#'
#' @export
evolution_clone_abundance_transcriptomic_corr_samplewise <- function(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, abundance_file_prefix, output_file_prefix, cell_cnt_thr=10){
	sample_list <- as.character(meta_info[,sample_info])
	sample_list_uniq <- unique(sample_list)

	pat_list <- as.character(meta_info[,pat_info])
	sample_pat <- table(sample_list, pat_list)
	pat_list_uniq <- colnames(sample_pat)

	sample_pat_map <- apply(sample_pat, 1, function(x) which(x!=0))

	clone_list <- as.character(meta_info[,"cancer_clone"])
	sample_clone <- table(sample_list, clone_list)
	clone_list_uniq <- colnames(sample_clone)

	sample_clone_map <- apply(sample_clone, 1, function(x) which(x!=0))


	for(j in c(1:length(sample_pat_map))){
		sample <- names(sample_pat_map)[j]
		print(sample)
		sample_size <- sum(sample_clone[j,])
			if (sample_size < cell_cnt_thr){
			next
		}


		pat <- pat_list_uniq[sample_pat_map[j]]
		sample_clone_list <- as.character(colnames(sample_clone)[sample_clone_map[sample][[1]]])

		clone_path <- read.table(paste(clone_path_file_prefix, pat, "_clone_path", sep =""))
		clone_path <- as.matrix(clone_path)
		clone_path[1,] <- as.character(clone_path[1,])
		clone_path[2,] <- as.character(clone_path[2,])

		clone_path_graph <-graph_from_edgelist(clone_path)

		sample_clone_path <- subgraph(clone_path_graph, sample_clone_list)
			if(length(E(sample_clone_path))==0){
			next
		}

		total_path <- clone_path_list(sample_clone_path)
		

		if (length(total_path) <1){

			next
		}

		
		a <- read.csv(paste(abundance_file_prefix, name, "_clone_ratio.csv", sep = ""))
		a1<-a[, 2:dim(a)[2]]
		a2<-as.matrix(a1)
		rownames(a2)<-a$X
		a2_1<-data.frame(a2)

		
		a2 <- read.delim(paste(feature_file_prefix, name, "_cell_addmodule", sep = ""),sep = "\t", check.names = F)

		total_path_order <- c()
		total_path_scc <- c()

		pp <- 1:length(total_path)
		for (k in c(1:length(total_path))){
			path_order <- total_path[[k]]
			path_order1<- as.character(names(path_order))

			tmp <- paste(path_order1, collapse = "_")
			path_order_name <- paste(sample, tmp, sep = "__")
			
			path_order <- paste(pat, path_order1, sep ="_")
			path_order_abundance <- paste("clone", path_order1, sep ="_")
			path_order_abundance2  <- paste(pat, path_order1, sep ="_")
			a3_1 <- a2_1[sample,path_order_abundance]



			colnames(a3_1) <- path_order_abundance 

			clone_a3 <- c()
			a3 <- c()
			for (z in c(1:length(path_order1))){
				tmp_path <- path_order1[z]
				clone_loc <- which(sample_list == sample & clone_list == tmp_path)
				tmp_a3 <- a2[,clone_loc]
				if (length(clone_loc) == 1){
					tmp_clone_a3 <- path_order_abundance[z]

				} else{

					tmp_clone_a3<-rep(path_order_abundance[z], dim(tmp_a3)[2])
				}
				clone_a3 <- c(clone_a3, tmp_clone_a3)
				a3<-cbind2(a3, as.matrix(tmp_a3))
			}
			rownames(a3) <- rownames(a2)
	
	
			mut_feature_order <-order(as.numeric(a3_1))

			new_path_order <- names(a3_1)[mut_feature_order]

			single_path_coeff <- c()
		
			for(m in c(1: dim(a3)[1])){
#print(m)
				b <- cbind(a3[m,],clone_a3)
				colnames(b) <- c("value", "clone")
				b<-data.frame(b)
				b$value <- as.numeric(b$value)
				bbb<- melt(a3_1, b$clone)
				b$abundance <- as.numeric(bbb)
				spearman <- cor.test(b$abundance, b$value, method = "spearman")
				coeff<- spearman$estimate
		
				single_path_coeff <- c(single_path_coeff, coeff)
				names(single_path_coeff)[length(single_path_coeff)] <- rownames(a3)[m]
	
			}


			total_path_scc  <- cbind(total_path_scc, single_path_coeff)
			colnames(total_path_scc)[k] <- path_order_name

		}

		write.table(total_path_scc, paste(output_file_prefix, name, "_", sample, "_clone_abundance_transcriptomic_scc_sample", sep = ""), quote = F, sep = "\t") 

	}
}



