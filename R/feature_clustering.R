#' avg_funct
#'
#' averaging multiple feature scores
#'
#' @param a1 path(row)-feature(column) table
#'
#' @return mean feature score for each path
#'
#' @examples
#' a<-total_table[,1:feature_size,drop=F]
#' a1<-a[,cluster_member, drop = F]
#' a2 <- avg_funct(a1)
#'
#' @export
avg_funct <- function(a1){
	new_mean <- c()
	for (k in c(1:dim(a1)[1])){
		val_list <- as.numeric(a1[k,])
		val_list <- na.omit(val_list)
		val_mean <- mean(val_list)
	if (is.nan(val_mean) ){
		val_mean <- NA
	}
		new_mean <- c(new_mean, val_mean)
	}
	names(new_mean) <- rownames(a1)
	return (new_mean)
}

#' feature_graph
#'
#' make feature table into igraph object
#'
#' @param total_table result of generate_table_for_analysis
#' @param feature_size number of features
#'
#' @return igraph object of features
#'
#' @examples
#' result <- generate_table_for_analysis(name, meta_info, pat_info, sample_info, meta_info_list,file_prefix,file_suffix)
#' total_table_pre <- result[[1]]
#' feature_size_pre <- result[[2]]
#' g3 <-feature_graph(total_table_pre,  feature_size_pre)
#'  
#' @export
feature_graph <- function(total_table, feature_size){

	a<-total_table[,c(1:feature_size)]
	b<-total_table[,c((1+feature_size):dim(total_table)[2])]

	g<-dist(t(a))
	g1<-as.matrix(g)
	g2<-graph_from_adjacency_matrix(g1, mode = "undirected", weighted = T)
	w<-get.edge.attribute(g2)
	w<-as.numeric(w$weight)
	g3<-delete.edges(g2, which(is.na(w)))
	g3 <- simplify(g3)
	return (g3)
}


#' feature_clustering_connected_comp
#'
#' clustering features by connected component after thresholding
#'
#' @param g3 result of feature_graph
#' @param total_table result of generate_table_for_analysis
#' @param feature_size number of features
#' @param z_index conduct z normalization on distances for the thresholding\
#' \itemize{
#' \item TRUE: use z-score
#' \item FALSE: use raw distance score
#' }
#' @param z_thr threshold of z-score
#' @return list(new_a, module_info, g6, g4)
#' \itemize{
#' \item new_a: new table of mean feature score in each clusters  
#' \item module_info: list object contains cluster member information for each feature
#' \item g6: igraph object after thresholding
#' \item g4: g4<-components(g6)
#' }
#'
#' @examples
#' new_result <- feature_clustering_connected_comp(g3, total_table_pre, feature_size_pre)
#' @export
#g3: raw graph
feature_clustering_connected_comp <- function(g3, total_table, feature_size, z_index =T, z_thr = 2.576){
	a<-total_table[,1:feature_size,drop=F]
	b<-total_table[,c((feature_size +1):dim(total_table)[2])]
	w<-E(g3)$weight
	z1<-(w-mean(w))/sd(w)

	if (z_index){
		g5 <- delete.edges(g3, which(z1 > (-1*z_thr)))
	} else{
		g5 <- delete.edges(g3, which(w > z_thr))
	}
	g6<-delete.vertices(simplify(g5), degree(g5)==0)
	similarity_w <- exp(-E(g6)$weight)
	E(g6)$weight <- similarity_w

	g4 <- components(g6)


	member <- g4$membership
	member_uniq <- unique(member)
	new_a <- c()
	module_info <- list()
	for (j in c(1:length(member_uniq))){
		member_index <- member_uniq[j]
		cluster_member <- names(g4$membership[which(g4$membership == member_index)])
		a1<-a[,cluster_member, drop = F]
		a2 <- avg_funct(a1)
		new_a <- cbind(new_a, a2)
		colnames(new_a)[j] <- paste("module", member_index, sep = "_")
		module_info[[j]] <- cluster_member
		names(module_info)[j] <- paste("module", member_index, sep = "_")

	}
	new_a <- cbind(data.frame(new_a), data.frame(b))
	return (list(new_a, module_info,g6,g4))
}


#' clustering_plot
#'
#' plot clustered feature graph
#'
#' @param g6 clustered feature igraph object 
#' @param clustering_res module_info from clustering function
#' @param name tilte of the data
#' @param layout layout for plot
#' \itemize{
#' \item "layout_with_fr"
#' \item "layout_with_graphopt"
#' }
#' @param format png or pdf
#' @param width width of the plot
#' @param height height of the plot
#' @param resolution resolution of png plot
#' @return image
#'
#' @examples
#' [1: plotting connected componenets]
#' new_result <- feature_clustering_connected_comp(g3, total_table_pre, feature_size_pre) 
#' total_table_new<-new_result[[1]]
#' module_info_new <- new_result[[2]]
#' feature_size_new <- length(module_info_new)
#' g6<-new_result[[3]]
#' clustering_res <- new_result[[4]]
#' clustering_plot(g6,clustering_res,"myplot")
#'
#' [2: plotting louvain clustering]
#' new_result1 <- feature_clustering_louvain(g6, total_table_pre, feature_size_pre, resolution = 1.1)
#' total_table<-new_result1[[1]]
#' module_info <- new_result1[[2]]
#' feature_size <- length(module_info)
#' g3<-new_result1[[3]]
#' g4<-new_result1[[4]]
#' clustering_plot(g3,g4,"myplot")
#' @export
clustering_plot <- function(g6, clustering_res, name, layout = "layout_with_fr", format = "png", width = 1600, height = 900, resolution = 60){
	v<-names(V(g6))	
	zz<-order(clustering_res$membership)
	zz1 <- clustering_res$membership[zz]
	zz2<-v[zz]


	if (layout == "layout_with_fr"){
		coords <- layout_with_fr(g6, weights =E(g6)$weight)
	} else if (layout == "layout_with_graphopt"){
		coords <- layout_with_graphopt(g6)
	}
#force directed
	
	if (format == "png"){
		png(name, width = width, height = height, res = resolution)
	} else if (format == "pdf"){
		pdf(name, width = width, height = height)
	}
	par(mar = c(2, 10, 2, 0))

	plt <- plot(g6, layout = coords,  vertex.color = rainbow(length(unique(clustering_res$membership)), alpha=0.5)[clustering_res$membership], vertex.label = match(v, zz2))
	legend("topleft",legend=paste(c(1:length(v)),zz2, sep = "__"), col = rainbow(length(unique(clustering_res$membership)), alpha=0.5)[zz1], pch = 19, inset=c(-0.05,0))
	print(plt)
	dev.off()



}



#' feautre_clustering_louvain
#'
#' Louvain clustering for feature graph
#'
#' @param g3 result of feature_graph
#' @param total_table result of generate_table_for_analysis
#' @param feature_size number of the features
#' @param resolution resolution for louvain clustering
#'
#' @return list(new_a, module_info, g3, g4)
#' \itemize{
#' \item new_a: new table of mean feature score in each clusters
#' \item module_info: list object contains cluster member information for each feature
#' \item g3: initial input igraph object
#' \item g4: g4<-components(g6)
#' }
#'
#'
#' @examples
#'
#' @export
#g3: pruned network
feature_clustering_louvain <- function(g3, total_table, feature_size, resolution=1.1){
	a<-total_table[,1:feature_size,drop=F]
	b<-total_table[,c((feature_size +1):dim(total_table)[2])]
	w<-get.edge.attribute(g3)
	w<-as.numeric(w$weight)

	g4<-cluster_louvain(g3, weights = w, resolution = resolution)
#g4<-cluster_louvain(g3, resolution = resolution)

	member <- g4$membership
	member_uniq <- unique(member)
	new_a <- c()
	module_info <- list()
	for (j in c(1:length(member_uniq))){
		member_index <- member_uniq[j]
		cluster_member <- g4$names[which(g4$membership == member_index)]
		a1<-a[,cluster_member, drop = F]
		a2 <- avg_funct(a1)
		new_a <- cbind(new_a, a2)
		colnames(new_a)[j] <- paste("module", member_index, sep = "_")
		module_info[[j]] <- cluster_member
		names(module_info)[j] <- paste("module", member_index, sep = "_")

	}
	new_a <- cbind(data.frame(new_a), data.frame(b))
	return (list(new_a, module_info,g3, g4))
}


