#' sample_table_making
#'
#' making table for the statistic analysis in evolution analysis
#'
#' @param tmp_file_name file name of evolution analysis
#' @param sample_name sample name
#' @param sample_position index of given sample from total single-cell data
#' @param meta_table table for meta information
#'
#' @return table for the statistic analysis in evolution analysis
#'
#' @examples
#' sample_list <- as.character(meta_info[,sample_info])
#' sample_list_uniq <- unique(sample_list)
#' meta_table<- c()
#' for (j in c(1:length(meta_info_list))){
#'  meta_info_tmp <- meta_info_list[j]
#'  meta_list <- as.character(meta_info[,meta_info_tmp])
#'  meta_table <- rbind(meta_table, meta_list)
#'  rownames(meta_table)[j] <- meta_info_tmp
#' }
#' sample_name <- sample_list_uniq[1]
#' sample_position <- which(sample_list == sample_name)
#' tmp_file_name <- paste(file_prefix, sample_name,file_suffix, sep = "")
#' a2<-sample_table_making(tmp_file_name, sample_name, sample_position, meta_table)
#'
#' @export
sample_table_making <- function(tmp_file_name, sample_name, sample_position, meta_table){
	a <- read.delim(tmp_file_name, check.names = F, sep = "\t")
	a2 <- t(a)
	a2<-cbind(a2, sample_name)
	colnames(a2)[dim(a2)[2]] <- "sample"

	tmp_meta_info <- c()

	for(k in c(1:dim(meta_table)[1])){
		tmp_meta_list <- meta_table[k,]
		tmp_meta_name <- rownames(meta_table)[k]
		sample_meta <- tmp_meta_list[sample_position[1]]
		tmp_meta_info <- cbind(tmp_meta_info, rep(sample_meta, dim(a2)[1]))
		colnames(tmp_meta_info)[k] <- tmp_meta_name
	
	}
	a2<-cbind(a2, tmp_meta_info)
	return (a2)
}





#' sample_table_making_mut_trans
#'
#' making table for the statistic analysis in evolution analysis (mut_trans analysis)
#'
#' @param tmp_file_name file name of evolution analysis
#' @param sample_name sample name
#' @param sample_position index of given sample from total single-cell data
#' @param meta_table table for meta information
#'
#' @return table for the statistic analysis in evolution analysis (mut_trans analysis)
#'
#' @examples
#' sample_list <- as.character(meta_info[,sample_info])
#' sample_list_uniq <- unique(sample_list)
#' meta_table<- c()
#' for (j in c(1:length(meta_info_list))){
#'  meta_info_tmp <- meta_info_list[j]
#'  meta_list <- as.character(meta_info[,meta_info_tmp])
#'  meta_table <- rbind(meta_table, meta_list)
#'  rownames(meta_table)[j] <- meta_info_tmp
#' }
#' sample_name <- sample_list_uniq[1]
#' sample_position <- which(sample_list == sample_name)
#' tmp_file_name <- paste(file_prefix, sample_name,file_suffix, sep = "")
#' a2<-sample_table_making_mut_trans(tmp_file_name, sample_name, sample_position, meta_table)
#'
#' @export
sample_table_making_mut_trans <- function(tmp_file_name, sample_name, sample_position, meta_table){
	a <- read.delim(tmp_file_name, check.names = F, sep = "\t")
	b<-rownames(a)
	b1<-strsplit(b, "____")

	b_loc <- c()
	new_row <- c()
	for (k in c(1:length(b))){
		if(b1[[k]][1] == b1[[k]][2]){
			b_loc <- c(b_loc, k)
			new_row <- c(new_row, b1[[k]][1])
		}

	}
	a1<-a[b_loc,, drop = F]
	
	rownames(a1) <- new_row
	a1<-data.frame(a1)
	#a1<-a1[gmt_list_input,,drop = F]
	
	#rownames(a1)<-gmt_list_input

	a2 <- t(a1)
	a2<-cbind(a2, sample_name)
	colnames(a2)[dim(a2)[2]] <- "sample"

	tmp_meta_info <- c()

	for(k in c(1:dim(meta_table)[1])){
		tmp_meta_list <- meta_table[k,]
		tmp_meta_name <- rownames(meta_table)[k]
		sample_meta <- tmp_meta_list[sample_position[1]]
		tmp_meta_info <- cbind(tmp_meta_info, rep(sample_meta, dim(a2)[1]))
		colnames(tmp_meta_info)[k] <- tmp_meta_name
	
	}
	a2<-cbind(a2, tmp_meta_info)
	
	return (a2)
}




#' generate_table_for_analysis
#'
#' making a table for statistic analysis
#ult <- generate_table_for_analysis(name, meta_info, pat_info, sample_info, meta_info_list,file_prefix, file_suffix)
#' @param name title for the data
#' @param meta_info meta information
#' @param pat_info variable name for patient
#' @param sample_info variable name for sample
#' @param meta_info_list list of meta information to analyze
#' @param file_prefix prefix of the file
#' @param file_suffix suffix of the file
#' @param mut_trans_index whether feature is mut-trans coassociation 
#' @return total_table for statistical analysis
#'
#' @examples
#' #meta_info <- seurat_obj@meta.data
#' meta_info <- read.table("data/meta_info_lung", sep = "\t", check.names = F)
#' name <- "lung"
#' pat_info <- "pat_collapsed"
#' sample_info <- "sample_name"
#' meta_info_list <- c("primary_or_metastaic","biopsy_site","smokingHx","histolgy","best_rxn_status")
#' file_prefix <- "./"
#' file_suffix <- "_clone_abundance_mut_scc_sample"
#' result_suffix <- "abundance_mut_scc_analysis"
#' result <- generate_table_for_analysis(name, meta_info, pat_info, sample_info, meta_info_list,file_prefix, file_suffix)
#' total_table <- result[[1]]
#' feature_size <- result[[2]]
#' @export
generate_table_for_analysis <- function(name, meta_info, pat_info, sample_info, meta_info_list,file_prefix,file_suffix, mut_trans_index = FALSE){
	sample_list <- as.character(meta_info[,sample_info])
	sample_list_uniq <- unique(sample_list)

	meta_table<- c()
	for (j in c(1:length(meta_info_list))){
		meta_info_tmp <- meta_info_list[j]
		meta_list <- as.character(meta_info[,meta_info_tmp])
		meta_table <- rbind(meta_table, meta_list)
		rownames(meta_table)[j] <- meta_info_tmp
	}

	meta_length <- length(meta_info_list)
	total_table <- c()

	for (j in c(1:length(sample_list_uniq))){
		sample_name <- sample_list_uniq[j]
		sample_position <- which(sample_list == sample_name)


		tmp_file_name <- paste(file_prefix, name, "_", sample_name,file_suffix, sep = "")
		if(file.exists(tmp_file_name)){
			if (mut_trans_index){
				a2<-sample_table_making_mut_trans(tmp_file_name, sample_name, sample_position, meta_table)
			} else{
				a2<-sample_table_making(tmp_file_name, sample_name, sample_position, meta_table)
			}
			a2<-data.frame(a2)


			if (j != 1){
				p1<-a2[,1:(dim(a2)[2] - meta_length-1)]
				p1_2 <- a2[(dim(a2)[2] - meta_length): dim(a2)[2]]
	                        p2<-total_table[,1:(dim(total_table)[2] - meta_length-1)]
        	                p2_2 <- total_table[(dim(total_table)[2] - meta_length): dim(total_table)[2]]

                	        new_col <- union(colnames(p1), colnames(p2))


				if(length(new_col) != dim(p1)[2]){
					p1<-t(p1)
					p1<-data.frame(p1)
					p1_1<-p1[new_col,]
					rownames(p1_1)<-new_col
					p1<-t(p1_1)
					a2 <- cbind(p1, p1_2)

				}
				if (length(new_col) != dim(p2)[2]){
					p2<-t(p2)
                                        p2<-data.frame(p2)
                                        p2_1<-p2[new_col,]
                                        rownames(p2_1)<-new_col
                                        p2<-t(p2_1)
                                        total_table <- cbind(p2, p2_2)
						
				}
			}
			total_table <- rbind(total_table, a2)

		}
	}

	tmp_n <- 1+ length(meta_info_list)
	feature_size <- dim(total_table)[2] - tmp_n
	total_table<-data.frame(total_table)

	for(m in c(1:(dim(total_table)[2] -tmp_n))){
		total_table[,m] <- as.numeric(total_table[,m])
	}
	for(m in c((dim(total_table)[2] -tmp_n +1):(dim(total_table)[2]))){
		total_table[,m] <- as.factor(as.character(total_table[,m]))
	}

	return (list(total_table, feature_size))
}


#' statistic_test
#'
#' statistical test with a given meta information in the dataset (+ plotting)
#'
#' @param total_table result of generate_table_for_analysis
#' @param feature_size number of the features
#' @param meta_info_input name of the meta information of interest
#' @param method wilcoxon-rank sum test or Kruskal-wallis test
#' @param name tile of the data
#' @param result_suffix suffix of the result file
#' @param plot_index plot (TRUE or FALSE)
#' @param format png or pdf
#' @param width width of the plot
#' @param height height of the plot
#' @param resolution resolution of png plot
#' @return statistical result and image (if applied)
#'
#' @examples
#' result_suffix <- "evolution_mut_scc_pathwise_analysis"
#' result <- generate_table_for_analysis(name, meta_info, pat_info, sample_info, meta_info_list,file_prefix, file_suffix)
#' total_table <- result[[1]]
#' feature_size <- result[[2]]
#' statistic_res<-statistic_test(total_table, feature_size, "biopsy_site", "kruskal", "lung", result_suffix, T)
#' 
#' statistic_res<-statistic_test(total_table, feature_size, "histolgy", "wilcox", "lung", result_suffix, T)
#'
#' @export
statistic_test<-function(total_table, feature_size, meta_info_input, method, name, result_suffix, plot_index, format = "png", width =640, height = 1200, resolution = 80){
	tmp_meta_list <- total_table %>% dplyr::select(starts_with(meta_info_input))
	meta_uniq <- as.character(unique(tmp_meta_list[,1]))
	colnames(tmp_meta_list) <-  "Meta"
	if (length(unique(tmp_meta_list$Meta))<2){
		print("The data has only one group")
		return(NA)
	}
 
	statistic_res<-data.frame()
	for (m in c(1:feature_size)){
		tmp_feature <-colnames(total_table)[m]
		tmp_table <- total_table[,m]

	
		tmp_table1 <- cbind(tmp_table, tmp_meta_list)
		colnames(tmp_table1)<-c("Value", "Meta")
		tmp_table1<-data.frame(tmp_table1)
		tmp_table1$Value<-as.numeric(tmp_table1$Value)
		tmp_table1$Meta<-as.factor(tmp_table1$Meta)

		tmp_table1 <- na.omit(tmp_table1)
		if (length(unique(tmp_table1$Meta))<2){
			next
		}
		if(method == "wilcox"){
			test_res <- wilcox.test(Value~Meta, tmp_table1, alternative = "two.sided")
			pval_res <- test_res$p.value
		} else if (method == "kruskal"){
			test_res <- kruskal.test(Value~Meta, tmp_table1)
			pval_res <- test_res$p.value
		}
		a<-aggregate(tmp_table1$Value, list(tmp_table1$Meta), mean)
		bbb <- match(meta_uniq, a$Group.1)
		a<-a[bbb,]
		a$Group.1 <- meta_uniq


		max_meta<-a$Group.1[which.max(a$x)]

		tmp_res <- c(pval_res, a$x, as.character(max_meta))
		statistic_res<-rbind(statistic_res, tmp_res)

		rownames(statistic_res)[dim(statistic_res)[1]] <- tmp_feature

	}
	colnames(statistic_res) <- c("pval", as.character(a$Group.1), "max_meta")

	for (j in c(1:(dim(statistic_res)[2]-1))){
		statistic_res[,j] <- as.numeric(statistic_res[,j])
	}

	statistic_res[,dim(statistic_res)[2]]<-as.factor(statistic_res[,dim(statistic_res)[2]])

	qval <- p.adjust(statistic_res$pval, method = "BH")
	statistic_res <- cbind(statistic_res, -log10(qval))
	statistic_res$pval <- -log10(statistic_res$pval)
	colnames(statistic_res)[dim(statistic_res)[2]] <- "qval"
	statistic_res <- cbind(statistic_res, rownames(statistic_res))
	colnames(statistic_res)[dim(statistic_res)[2]] <- "feature"

	statistic_res<-statistic_res[order(statistic_res$pval, decreasing = T),]

	write.csv(statistic_res, file = paste(name, "_", meta_info_input, "_", result_suffix, "_comparison.csv", sep =""), quote = F)
	if (plot_index){
		if (format == "png"){
			png(paste(name, "_", meta_info_input, "_", result_suffix, "_comparison.png", sep =""), width = width, height= height, res = resolution)
		} else if (format == "pdf"){
			pdf(paste(name, "_", meta_info_input, "_", result_suffix, "_comparison.pdf", sep =""), width = width, height= height)
		}
		plt <- ggplot(statistic_res,aes(y=qval, x=reorder(feature, qval), fill = max_meta)) + geom_bar(stat = "identity") + coord_flip() +labs(y = "-log10(qval)", x = "Feature")


		print(plt)
		dev.off()
	}
	return (statistic_res)
}

