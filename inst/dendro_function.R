dendro_process1 <- function(name){
	dlist <- fread(paste(name, '.vcf.gz',sep = ""),  header =T, skip = '#CHROM')
	print("reading file done")
	chr_list=sapply(seq(1,22),function(x)paste0('chr',x),simplify = T)
	dlist <- dlist %>% filter(`#CHROM`%in%chr_list)
	k<-sapply(dlist[,8], function(x) (!grepl("COSMIC",x) & grepl("RS=",x)|grepl("RNAEDIT",x)))
	k1<-!k
	k2<-which(k1)
	dlist<-dlist[k2,]
	q<-sapply(dlist[,8], function(x) (grepl("COSMIC",x)))
	q1<-sapply(dlist[,8], function(x) (strsplit(x, "|", fixed = T)))

	tmp1 <- c()
	tmp2 <- c()
	for(j in c(1: length(q1))){
		tmp1 <- c(tmp1,q1[[j]][2])
		tmp2 <- c(tmp2, q1[[j]][4])
	}
	q2 <- cbind(q, tmp1, tmp2)
	colnames(q2) <- c("cosmic", "mutation_type", "gene")

	Format <- dlist%>%select(FORMAT)
	Info <- dlist %>% select(`#CHROM`,POS,REF,ALT)
	dlist <- dlist %>% select(-`#CHROM`,-POS,-ID,-REF,-ALT,-QUAL,-FILTER,-INFO,-FORMAT)

	GT_pos=1
	AD_pos=2
	DP_pos=3

	exinfo_x <- function(info){
		  return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,as.numeric(strsplit(x[AD_pos],',')[[1]][2]))},simplify=T))
	}
	exinfo_n <- function(info){
		return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,as.numeric(x[DP_pos]))},simplify=T))
	}
	exinfo_z <- function(info){
		return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,ifelse(x[GT_pos]!='0/0',1,0))},simplify=T))
	}

	X=sapply(dlist,exinfo_x)
	N=sapply(dlist,exinfo_n)
	Z=sapply(dlist,exinfo_z)

	thres=floor(0.01*ncol(X))
	sel=rowSums(!is.na(Z))>thres
	X=X[sel,]
	N=N[sel,]
	Z=Z[sel,]

	X[is.na(X)]<-0
	N[is.na(N)]<-0
	q2<-q2[sel,]
	Info <- cbind(Info, q2)

	demo_qc = FilterCellMutation(X,N,Z,Info, cut.off.VAF = 0.05, cut.off.sd = 5)
	print("variant qc done")

	demo_qc$dist = DENDRO.dist(demo_qc$X,demo_qc$N,demo_qc$Z,show.progress=FALSE)
	demo_qc$cluster = DENDRO.cluster(demo_qc$dist,label=demo_qc$label,type='fan')

	png(paste(name, "_icd.png", sep =""))
	demo_qc$icd = DENDRO.icd(demo_qc$dist,demo_qc$cluster, kmax = dim(demo_qc$X)[2])
	print(demo_qc$icd)
	dev.off()

	save(demo_qc, file = paste(name, "_demo_qc", sep = ""))
}

###########
# assigning K

dendro_process2 <- function(name, opt){
	load(paste(name, "_demo_qc", sep = ""))

	demo_qc$optK = opt
	demo_qc$DENDRO_label = cutree(demo_qc$cluster,demo_qc$optK)
	demo_qc$cluster = DENDRO.cluster(demo_qc$dist,label=demo_qc$DENDRO_label,type='fan')

	demo_cluster = jw.recalculate(demo_qc$X,demo_qc$N, demo_qc$Info, demo_qc$DENDRO_label, cluster.name= paste(name, c(1:demo_qc$optK), sep = "_"))
	write.csv(demo_cluster$Info, paste(name, "_mutation_info.csv", sep = ""), quote = F)

	save(demo_cluster, file = paste(name, "_final_dendro", sep = ""))
	return(demo_qc$DENDRO_label)
}

jw.recalculate = function(X,N,Info,DENDRO_label,cluster.name=NULL,top=NULL,epi = 0.001,m=2){
	if(is.null(cluster.name)){
		cluster.name=as.character(1:length(unique(DENDRO_label)))
	}

	N_cluster = mapply(function(k){
	sel=which(DENDRO_label==k)
	if(length(sel)==1){
		return(N[,sel])
	}
	return(rowSums(N[,sel],na.rm=T))
	},1:length(unique(DENDRO_label)),SIMPLIFY = T)
	colnames(N_cluster)=cluster.name

	X_cluster = mapply(function(k){
	sel=which(DENDRO_label==k)
	if(length(sel)==1){
		return(X[,sel])
	}
	return(rowSums(X[,sel],na.rm=T))
	},1:length(unique(DENDRO_label)),SIMPLIFY = T)
	colnames(X_cluster)=colnames(N_cluster)

	lg = function(m,k,l,g,epi){
		-k*log(m)+l*log((m-g)*epi+g*(1-epi))+(k-l)*log((m-g)*(1-epi)+g*epi)
	}

	lg_array = array(data=NA,dim=c(dim(X_cluster),3))
	lg_array[,,1] = lg(m,N_cluster,N_cluster-X_cluster,0,epi)
	lg_array[,,2] = lg(m,N_cluster,N_cluster-X_cluster,1,epi)
	lg_array[,,3] = lg(m,N_cluster,N_cluster-X_cluster,2,epi)

	Z_cluster=apply(lg_array,c(1,2),function(x){
	if(all(is.na(x))){return(c(NA,NA))}
		else{return(c(which.max(x),max(x,na.rm=T)))}
	})

	Z_cluster_lg=Z_cluster[2,,]
	Z_cluster=Z_cluster[1,,]

	for(j in c(1:dim(Z_cluster_lg)[2])){
		w<-which(Z_cluster_lg[,j] == 0)
		Z_cluster[w,j] <- NA
	}
	Z_cluster<-Z_cluster - 1
	colnames(Z_cluster)=colnames(X_cluster)
	cat('Before QC, there are total ',nrow(Z_cluster),' mutations across ',ncol(Z_cluster),' subclones \n')

# Filter
	not_sel=apply(Z_cluster,1,function(x)all(x[1] == x))
	not_sel[is.na(not_sel)]=FALSE
	Z_cluster = Z_cluster[!not_sel,]
# Z_cluster_lg = Z_cluster_lg[!not_sel,]
	N_cluster= N_cluster[!not_sel,]
	X_cluster= X_cluster[!not_sel,]
	Info_cluster=Info[!not_sel,]

	if(!is.null(top)){
		h_score = order(rowSums(abs(Z_cluster_lg)),decreasing = TRUE)[1:top]
		Z_cluster=Z_cluster[h_score,]
		hdlist_x_cluster_qc=dlist_x_cluster_qc[h_score,]
		N_cluster=N_cluster[h_score,]
		X_cluster=X_cluster[h_score,]
		Info_cluster=Info_cluster[h_score,]
	}
	cat('After QC, there are total ',nrow(Z_cluster),' mutations across ',
	ncol(Z_cluster),' subclones \n')

	return(list(X=X_cluster,N=N_cluster, Z=Z_cluster,Info=Info_cluster))
}


