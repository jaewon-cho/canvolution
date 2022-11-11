#' select_path_by_thr
#'
#' select evolution paths by a given threshold of a given feature
#'
#' @param total_table output of evolution path analysis
#' @param sample sample of interest
#' @param feature feature of interest 
#' @param thr threshold for feature score
#' @param sign selecting higher or lower paths than the threshold
#' \itemize{
#' \item high (>thr)
#' \item low (<= thr)
#' }
#' @return selected evolution path vector
#'
#' @examples
#' path_list <- select_path_by_thr(total_table,  "LT_S21", "Angiogenesis",1, "high")
#' @export
select_path_by_thr <- function(total_table, sample, feature, thr, sign){
	p<-total_table[which(total_table$sample == sample),,drop=F]
	p1<-p[,feature,drop=F]
	if (sign == "high"){
		p2<-p1>thr
	} else{
		p2 <- p1 <= the
	}
	res <- rownames(p1)[which(p2)]
	return(res)
}



###
#file_name → obj name: sample_clone_mut_list

#' mutation_from_path
#'
#' mutation profile for each clone in a given evolution path 
#'
#' @param path evolution path 
#' @param sample sample name
#' @param sample_clone_mut_list list object which contains mutational profile as vector for each clone
#' @param feature_gn gene list for a given feature
#' \itemize{
#' \item "all": calling all the mutation from a given clone
#' \item gene_vector: calling only the mutation within a given gene list
#' }
#' @param root_exclue whether user wants to exclude mutation profile of root state or not
#' \itemize{
#' \item True 
#' \item False 
#' }
#'
#' @return list object of mutation profile for each clone in a given evolution path
#'
#' @examples
#' load("data/clone_mutation_info/LT_S21_clone_mut_list")
#' path_mut_gn_list <- mutation_from_path("LT_S21__29_23_27_25", "LT_S21", sample_clone_mut_list, "all", F)
#' path_mut_gn_list <- mutation_from_path("LT_S21__29_23_27_25", "LT_S21", sample_clone_mut_list, angiogenesis_gn_list, F)
#' @export
mutation_from_path <- function(path, sample, sample_clone_mut_list, feature_gn, root_exclude){
	clone_list<- strsplit(path, "__")
	clone_list<-clone_list[[1]][length(clone_list[[1]])]
	clone_list<-strsplit(clone_list, "_")
	clone_list<-clone_list[[1]]
	clone_list <- paste("clone", clone_list, sep = "_")

	res<-list()
	for (j in c(1:length(clone_list))){
		clone <- clone_list[j]
		mut_gn <- sample_clone_mut_list[[clone]]
		if (feature_gn != "all"){
			mut_gn <- intersect(mut_gn, feature_gn)
		}

		if (root_exclude){
			if (j==1){
				root_gn <- mut_gn
				next
			} else{
				mut_gn <- setdiff(mut_gn, root_gn)
			}
		}	
			
		res[[length(res)+1]] <- mut_gn
		names(res)[length(res)] <- clone
	}
	return (res)
}


#' common_mutation_from_path
#'
#' extracting overlapping mutational profile in multiple evolution path
#'
#' @param path_list_info list object with sample, path, sample_clone_mut_list for each sample
#' @param feature_gn gene list for a given feature
#' \itemize{
#' \item "all": calling all the mutation from a given clone
#' \item gene_vector: calling only the mutation within a given gene list
#' }
#' @param root_exclue whether user wants to exclude mutation profile of root state or not
#' \itemize{
#' \item True
#' \item False
#' }
#'
#' @return return overlapping genelist
#'
#' @examples
#' path_list_info<-list()
#' a<- c("LT_S21__29_23_27_25", "LT_S53__3_1_2")
#'
#' for (j in c(1:length(a))){
#' 	path <- a[j]
#' 	sample<-strsplit(a[1], "__")[[1]][1]
#' 	file_name <- paste(sample, "_clone_mut_list",sep = "")
#' 	load(file_name)
#' 	path_list_info[[j]] <- list(sample, sample_clone_mut_list)
#' 	names(path_list_info)[j] <- path
#' } 
#'
#' #names(path_list_info): path
#' #path_list_info[[1]][[1]]: sample
#' #path_list_info[[1]][[2]]: sample_clone_mut_list
#'
#' genes <- common_mutation_from_path(path_list_info, "all", F)
#' @export
common_mutation_from_path <- function(path_list_info, feature_gn, root_exclude){
	res <- c()
	for (j in c(1:length(path_list_info))){
		path <- names(path_list_info)[j]
		sample <- path_list_info[[j]][1]
		sample_clone_mut_list <- path_list_info[[j]][2]
		mut_list <- mutation_from_path(path, sample, sample_clone_mut_list, feature_gn, root_exclude)
		mut_tmp <- c()
		for (k in c(1:length(mut_list))){
			tmp <- mut_list[[k]]
			mut_tmp <- c(mut_tmp, tmp)
		}
		mut_tmp <- unique(mut_tmp)
		if(j == 1){
			res <- mut_tmp
		} else{
			res <- intersect(res, mut_tmp)
		}

	}
	return (res)

}




