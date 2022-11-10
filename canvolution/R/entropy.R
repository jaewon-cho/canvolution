#' shannon_entropy
#'
#' normalized shannon entropy
#'
#' @param x given clone or cluster
#' @param batch_vector cluster vector or clone vector (batch)
#' @param N_batches number of cluster or clone
#' @param cluster clone vector or cluster vector
#'
#' @return entropy
#'
#' @examples 
#' ent <- shannon_entropy(clone, cluster_vec, cluster_size, clone_vec)
#' @export
shannon_entropy <- function(x, batch_vector, N_batches, cluster) {
	freq_batch = table(batch_vector[which(cluster == x)])/length(batch_vector[which(cluster == x)])
	freq_batch_positive = freq_batch[freq_batch > 0]
	if (N_batches == 1){
		return (0)
	} else {
		return(-sum(freq_batch_positive * log(freq_batch_positive))/log(N_batches))
	}
}

#' random_entropy
#'
#' shannon entropy for random expectation
#'
#' @param cluster_vec cluster(clone) vector for a given sample
#'
#' @return entropy
#'
#' @examples
#' ent <- random_entropy(as.numeric(obj$seurat_clusters))
#' @export
random_entropy <- function(cluster_vec){
	#clone: NA exclude
	a<- table(cluster_vec)
	a1 <- a / sum(a)
	a2<- log(a1)
	a3 <- a1 * a2
	if(length(unique(cluster_vec))==1){
		return (0)
	} else{
		entropy <- -sum(a3) / log(length(unique(cluster_vec)))
		return (entropy)
	}
}


#' clone_shannon_process
#'
#' Evaluating entropy for each clone in a given sample
#'
#' @param tmp_d seurat object with "seurat_clusters" and "cancer_clone"
#' @param filter_cnt threshold of cell count for each clone
#' @param sample sample name
#'
#' @return entropy vector
#'
#' @examples
#' d3 <- seurat_obj
#' d3 <- FindClusters(d3, resolution = resolution, verbose = FALSE)
#' tmp_d <- subset(d3, subset = cancer_clone != "NA")
#' tmp_res1 <- shannon_process(tmp_d, 0, "sample1")
#'
#' @export
clone_shannon_process <- function(tmp_d, filter_cnt, sample){
	clone_list <- unique(tmp_d$cancer_clone)
	tmp_d$seurat_clusters <- as.character(tmp_d$seurat_clusters)

	cluster_vec <- tmp_d$seurat_clusters
	cluster_size <- length(unique(cluster_vec))
	clone_vec <- tmp_d$cancer_clone

	tmp_res <- c()
	tmp_res1 <- c()
	tmp_res2 <- c()

	for(k in c(1:length(clone_list))){
		clone <- clone_list[k]
		if (is.na(clone)){
			next
		}


		clone_size <- length(which(clone_vec == clone))
		if (clone_size < filter_cnt) {
			next
		}

		ent <- shannon_entropy(clone, cluster_vec, cluster_size, clone_vec)
		if (is.na(ent)){
			next
		}

		tmp_res <- c(tmp_res, ent)
		names(tmp_res)[length(tmp_res)] <- paste(sample, clone, sep ="__")
	}

	return (tmp_res)
}

#' cluester_shannon_process
#'
#' Evaluating entropy for each cluster in a given sample
#'
#' @param tmp_d seurat object with "seurat_clusters" and "cancer_clone"
#' @param filter_cnt threshold of cell count for each cluester
#' @param sample sample name
#'
#' @return entropy vector
#'
#' @examples
#' d3 <- seurat_obj
#' d3 <- FindClusters(d3, resolution = resolution, verbose = FALSE)
#' tmp_d <- subset(d3, subset = cancer_clone != "NA")
#' tmp_res1 <- shannon_process(tmp_d, 0, "sample1")
#'
#' @export
cluster_shannon_process <- function(tmp_d, filter_cnt, sample){
	tmp_d$seurat_clusters <- as.character(tmp_d$seurat_clusters)
	cluster_list <- unique(tmp_d$seurat_clusters)
	clone_vec <- tmp_d$cancer_clone
	clone_size <- length(unique(clone_vec))
	cluster_vec <- tmp_d$seurat_clusters


	tmp_res <- c()
	tmp_res1 <- c()
	tmp_res2 <- c()

	for(k in c(1:length(cluster_list))){
		cluster<- cluster_list[k]
		if (is.na(cluster)){
			next
		}


		cluster_size <- length(which(cluster_vec == cluster))
		if (cluster_size  < filter_cnt) {
			next
		}

		ent <- shannon_entropy(cluster, clone_vec, clone_size, cluster_vec)
		if (is.na(ent)){
			next
		}


		tmp_res <- c(tmp_res, ent)
		names(tmp_res)[length(tmp_res)] <- paste(sample, cluster, sep ="__")
	}

	return (tmp_res)
}


