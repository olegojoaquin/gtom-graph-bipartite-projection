devtools::install_github("nicolewhite/RNeo4j")
install.packages("igraph")
install.packages("MCL")
source("http://bioconductor.org/biocLite.R") 

BiocManager::install(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
install.packages("WGCNA") 
install.packages("dynamicTreeCut")
install.packages("RColorBrewer")
install.packages("gplots")
source("http://www.bioconductor.org/biocLite.R")

source("https://bioconductor.org/biocLite.R")
BiocManager::install("DOSE")
source("http://bioconductor.org/biocLite.R")
BiocManager::install("topGO")
install.packages("clusteval")

library(igraph)
library(Matrix)
#library(MCL)
library(WGCNA)
library(dynamicTreeCut)
library(parallel)
library("RColorBrewer")
library("gplots")
#library("clusteval")


#Levanto DF

library(readr)
fav_all <- read_csv("Documents/Moove-IT/Projects/YouScience/Clusterizacion-Favoritos/FAVORITES_ALL.csv")


g <- graph.empty(directed = F)
node.out <- unique(fav_all$USER_ID) #stringsAsFactor = F in data frame
node.in <- unique(fav_all$TITLE) #stringsAsFactor = F in data frame
g <- add.vertices(g,nv=length(node.out),attr=list(name=node.out),type=rep(TRUE,length(node.out)))
g <- add.vertices(g,nv=length(node.in),attr=list(name=node.in),type=rep(FALSE,length(node.in)))
edge.list.vec <- as.vector(t(as.matrix(data.frame(fav_all[,c(1,3)]))))
g <- add.edges(g,edge.list.vec)

#Type TRUE==USER_ID TYPE FALSE==TITLE

usernodes<-which(V(g)$type==TRUE)
titlenodes<-which(V(g)$type==FALSE)


ng<-length(usernodes)
nd<-length(titlenodes)
           



sqrtk<-1/sqrt(degree(g))[1:ng]

msqrtk1<-matrix(rep(sqrtk,each=nd),nrow=nd,ncol=ng)

colnames(msqrtk1)<-V(g)$name[1:ng]

A <- as_adj(g,sparse=TRUE)[(ng+1):(ng+nd),1:ng]  #matriz de adyacencia gen-enfermedad
A <- A * msqrtk1  #divido por sqrt del grado de cada gen 

aa<- Matrix::tcrossprod(A)

ddg <- graph_from_adjacency_matrix(aa,mode="undirected",weighted=TRUE,diag=FALSE)
V(ddg)$kgenes <- degree(g)[V(ddg)$name]

kikj<-apply(ends(ddg,E(ddg)),1,function(x){ 
  return(sum(1/get.vertex.attribute(ddg,"kgenes",x))/2)
})
E(ddg)$w <- E(ddg)$weight*kikj



adj<-get.adjacency(ddg,attr="w")
adj<-adj/max(adj)  #escaleo para que los pesos esten [0,1]

tom1<-GTOMdist(as.matrix(adj))

plot(ecdf(tom1[upper.tri(tom1)]))

hca <- hclust(as.dist(tom1),method="average")
#hcw <- hclust(as.dist(tom1),method="ward")

#Elijo linkage: average, porque ward me parece que fuerza clusters de tamanio similar
#hc <- hcw
hc <- hca

hc$height <- round(hc$height,6)

ctdwt0<-cutreeDynamic(hc,minClusterSize=1,method="tree",dist=tom1,
                      deepSplit=0)
ctdw1<-cutreeDynamic(hc,minClusterSize=1,method="hybrid",dist=tom1,
                     deepSplit=1)
ctdw2<-cutreeDynamic(hc,minClusterSize=1,method="hybrid",dist=tom1,
                     deepSplit=2)
ctdw3<-cutreeDynamic(hc,minClusterSize=1,method="hybrid",dist=tom1,
                     deepSplit=3) 
ctdw4<-cutreeDynamic(hc,minClusterSize=1,method="hybrid",dist=tom1,
                     deepSplit=4)
plotDendroAndColors(hc,colors=cbind(ctdw1,ctdw2,ctdw3,ctdw4),
                    groupLabels=c("H_DS1","H_DS2","H_DS3","H_DS4"),
                    dendroLabels=FALSE)



V(ddg)$ctdw1 = ctdw1



nodos <- as.vector(names(V(ddg)))
com_full <- as.vector(V(ddg)$ctdw1)

nodos_fav <- data.frame(nodos,com_full)
colnames(nodos_fav) <- c('TITLE','GTOM-COM')

write.csv(nodos_fav,'/home/joaquin/Documents/Moove-IT/Projects/YouScience/Clusterizacion-Favoritos/cluster.csv',row.names = FALSE)