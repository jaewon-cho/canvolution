############################## clonal minimum spanning tree ############################
# Input:
#  clone_gety: the inferred clonal genotype based on the Louvain-Jaccard clustering output by subclone_GTM function;
#  robust_clone: the clustering result output by LJClustering function;
#  type: the data type ('SNV' or 'CNV') of input.
#        Here, 'CNV' data element is the number of copies of each chromosome segment in each cell, such as: 0, 1, 2, 3, 4, 5 , ……, 
#            where copy number 2 is normal,
#       'SNV' data element is binary or ternary, that is 0, 1 or 0, 1, 2, 
#            where 0 represents normal, 1 in binary data represents mutation under the hypothesis of infinite site, and 1,2 in ternary data represent mutation under the hypothesis of finite site;
#  pdf_name: the name of pdf with the clonal MST graph.
# Output: 
#  the pdf with the clonal MST graph;
#  el: the connected edges in the MST.


plot_MST <- function(clone_gety, robust_clone, type, pdf_name){
  if(type == 'SNV'){
    clone_gety_root <- matrix(0,(length(robust_clone)+1),ncol(clone_gety))
  }
  if(type == 'CNV'){
    clone_gety_root <- matrix(2,(length(robust_clone)+1),ncol(clone_gety))
  }
  
  clone_gety_root[2:(length(robust_clone)+1),] <- clone_gety
  
  rownames(clone_gety_root) <- c('Root',paste('subclone',c(1:length(robust_clone)),sep=''))
  rownames(clone_gety) <- paste('subclone',c(1:length(robust_clone)),sep='')
  
  dis1<- dist(clone_gety_root,p=2)
  dis1 <- as.matrix(dis1)
  dis<- dist(clone_gety,p=2)
  
  spanningtree <- spantree(dis)
  rootid <- which.min(dis1[2:nrow(dis1),1])
  
  idspantree<-spanningtree$kid
  numedge<- nrow(clone_gety)-1
  spantree<-matrix(1,numedge,2)
  spantree[,1]<-as.matrix(2:nrow(clone_gety),nrow(clone_gety),1)
  spantree[,2]<-t(idspantree)
  
  adj<- matrix(0,nrow(clone_gety),nrow(clone_gety))
  adj[spantree]<-1
  adj<-(adj+t(adj))
  
  allpath<-findpath(rootid,rootid,info = adj)
  k1 <- 0
  for(i in 1:length(allpath)){
    k1 <- k1+(length(allpath[[i]])-1)
  }
  el <- matrix(0,k1,2)
  k2 <- 1
  for(i in 1:length(allpath)){
    subpath <- allpath[[i]]
    for(j in 1:(length(subpath)-1)){
      el[k2,1] <- subpath[j]
      el[k2,2] <- subpath[j+1]
      k2 <- k2+1
    }
  }
  el <- unique(el)
  colnames(el) <- c('from','to')
  minspantree <- igraph::graph.edgelist(el, directed = TRUE)
  labelname <- rownames(clone_gety)
  col <- intpalette(c('magenta','red','limegreen','skyblue','pink','brown','gold','blue','cyan'),length(robust_clone))
  size <- c(1:length(robust_clone))
  for(i in 1:length(size)){
    size[i] <- length(robust_clone[[i]])/sum(lengths(robust_clone))*250
  }
  
  pdf(paste('RobustClone_MST_', pdf_name, '.pdf', sep=''))
  plot(minspantree, layout=layout_as_tree, vertex.color=col, edge.color='black', vertex.label.color='black', alpha=0.5,
       edge.width=4, vertex.size=size, vertex.shape= "sphere", vertex.label=labelname)
  dev.off()
  
  return(el)
}


findpath<-function(node,prev_node,info){
  adj<-which(info[node,]!=0)
  kid<-setdiff(adj,prev_node)
  if(length(kid)==0)
    return(list(node))
  t<-list()
  l<-1
  for (i in kid) {
    path <- findpath(i, node, info)
    for (j in 1:length(path)) {
      t[[l]]<-c(node,path[[j]])
      l<-l+1
    }
  }
  return(t)
}


#robust_clone: list
#[[1]]: clone1: 1,4,5,6 (cell index)
#[[2]]: clone2: 2,3,7,8

#clone_gety: clone representative mutation profile

