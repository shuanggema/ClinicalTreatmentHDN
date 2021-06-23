rm(list=ls())

year<-"all"

source("C:/Users/hm336/Box Sync/Work/Network/Aim2/Code/Ancillary functions.R")
setwd("C:/Users/hm336/Box Sync/Work/Network/Aim2/Result")
load(paste("result_",year,".RData",sep=""))

p<-dim(X)[2]
name<-colnames(X)

beta_x<-sapply(1:p,function(j){return(insert(0,para.list[[j]][[1]],j))})
beta_y<-sapply(1:p,function(j){return(insert(0,para.list[[j]][[2]],j))})
beta_z<-sapply(1:p,function(j){return(insert(0,para.list[[j]][[3]],j))})

cutoff=0

adjacency<-diag(1,p,p)
for (i in 2:p) {
  for (j in 1:(i-1)) {
    ifelse(abs(beta_x[i,j])<=cutoff & abs(beta_x[j,i])<=cutoff &
             abs(beta_y[i,j])<=cutoff & abs(beta_y[j,i])<=cutoff &
             abs(beta_z[i,j])<=cutoff & abs(beta_z[j,i])<=cutoff,
           adjacency[i,j]<-adjacency[j,i]<-NA,
           adjacency[i,j]<-adjacency[j,i]<-1)
  }
} 

n.edge<-(sum(adjacency,na.rm = TRUE)-p)/2
n.edge
p.nzero=n.edge/(p*(p-1)/2)
p.nzero


library("WGCNA")

TOM = TOMsimilarity(adjacency)
tom_connectivity = apply(TOM,1,sum)
dissTOM = 1-TOM
colnames(dissTOM) = name

NEITree = hclust(as.dist(dissTOM), method = "average")

minModuleSize = 5
dynamicMods = cutreeDynamic(dendro = NEITree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = TRUE, minClusterSize = minModuleSize)
table(dynamicMods)

plotDendroAndColors(NEITree,colors=dynamicMods)
plotTOM = dissTOM
diag(plotTOM) = NA
TOMplot(plotTOM, NEITree, Colors =dynamicMods,
        main = "Network Heatmap Plot")
connectivity_matrix<-intramodularConnectivity(adjacency, dynamicMods)
ID = seq(1,nrow(TOM))
label = name
weight = tom_connectivity
module=dynamicMods
connectivity<-connectivity_matrix$kTotal
intra_connectivity<-connectivity_matrix$kWithin
net<-cbind(ID,label,module,weight,connectivity,intra_connectivity)
net=as.data.frame(net)
write.csv(net,file = paste("vinfo_",year,".csv", sep=""),fileEncoding = "GBK")

edges = as.data.frame(t(combn(nrow(adjacency),2)))
edges$weight=NA

for (i in 1:nrow(edges)) {
  a=edges[i,1]
  b=edges[i,2]
  edges[i,3]=adjacency[a,b]
}

edges=edges[-c(which(is.na(edges$weight))),]
#edges$weight=edges$weight*10
colnames(edges)[1:2]<-c("source","target")
write.csv(edges,file = paste("einfo_",year,".csv", sep=""),fileEncoding = "GBK")





