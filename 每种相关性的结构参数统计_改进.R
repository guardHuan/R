library(igraph)
#设置当前的工作目录
setwd('F:/system_modeling/cors_d/')
#library(infotheo)
#library(dplyr)
#获取当前目录下面的所有文件
corslist<-dir()[-6]
#normal_entropy<-c()
#tumor_entropy<-c()
seqs<-seq(0.01,1,0.01)
for(j in corslist){
  #j的值为文件的名称
  dcor_order<-read.csv(j)
  normal<-dcor_order[,c(2,3,4)]
  tumor<-dcor_order[,c(2,3,5)]
  #normal<-discretize(cors[,c(2,3,4)])
  #tumor<-discretize(cors[,c(2,3,5)])
  #normal_entropy<-c(normal_entropy,entropy(normal))
  #tumor_entropy<-c(tumor_entropy,entropy(tumor))
  #filter_normal_dcor_order<-filter_normal_dcor_order[,c(2,3,4)]
  #filter_tumor_dcor_order<-filter_tumor_dcor_order[,c(2,3,5)]
  dcor_tumor_vertex_num<-c()
  dcor_normal_vertex_num<-c()
  dcor_tumor_edges_num<-c()
  dcor_normal_edges_num<-c()
  dcor_tumor_avedegree_num<-c()
  dcor_normal_avedegree_num<-c()
  dcor_tumor_islands_num<-c()
  dcor_normal_islands_num<-c()
  dcor_tumor_transitivity_num<-c()
  dcor_normal_transitivity_num<-c()
  dcor_tumor_connectivity_num<-c()
  dcor_normal_connectivity_num<-c()
  dcor_tumor_diameter_num<-c()
  dcor_normal_diameter_num<-c()
  for(i in 1:100){
    #进行阈值筛选
    filter_tumor_dcor_order=tumor[tumor$tumor>=(i*0.01),]
    filter_normal_dcor_order=normal[normal$normal>=(i*0.01),]

    graph_tumor_dcor_order<-graph.data.frame(filter_tumor_dcor_order,directed = FALSE)
    graph_normal_dcor_order<-graph.data.frame(filter_normal_dcor_order,directed =FALSE)
    
    #number of node
    dcor_tumor_vertex_num[i]<-length(V(graph_tumor_dcor_order))
    dcor_normal_vertex_num[i]<-length(V(graph_normal_dcor_order))
    
    #number of edges
    dcor_tumor_edges_num[i]<-length(E(graph_tumor_dcor_order))
    dcor_normal_edges_num[i]<-length(E(graph_normal_dcor_order))
    
    #number of density
    dcor_tumor_avedegree_num[i]<-graph.density(graph_tumor_dcor_order)
    dcor_normal_avedegree_num[i]<-graph.density(graph_normal_dcor_order)
    
    #number of islands
    dcor_tumor_islands_num[i]<-clusters(graph_tumor_dcor_order)$no
    dcor_normal_islands_num[i]<-clusters(graph_normal_dcor_order)$no
    
    #transitivity
    dcor_tumor_transitivity_num[i]<-transitivity(graph_tumor_dcor_order, type="global")
    dcor_normal_transitivity_num[i]<-transitivity(graph_normal_dcor_order, type="global")
    
    #connectivity 
    dcor_tumor_connectivity_num[i]<-edge.connectivity(graph_tumor_dcor_order)
    dcor_normal_connectivity_num[i]<-edge.connectivity(graph_normal_dcor_order)
    
    # Diameter diam_normal<-diameter(s_normal_graph)
    dcor_tumor_diameter_num[i]<-diameter(graph_tumor_dcor_order)
    dcor_normal_diameter_num[i]<-diameter(graph_normal_dcor_order)

  }
  #dcor_threshold_vertex_tumor<-data.frame(threshold=seq(0.01,1,0.01),vertex_num=dcor_tumor_vertex_num)
  #dcor_threshold_vertex_normal<-data.frame(threshold=seq(0.01,1,0.01),vertex_num=dcor_normal_vertex_num)
  #dcor_threshold_edges_tumor<-data.frame(threshold=seq(0.01,1,0.01),edges_num=dcor_tumor_edges_num)
  #dcor_threshold_edges_normal<-data.frame(threshold=seq(0.01,1,0.01),edges_num=dcor_normal_edges_num)
  #dcor_threshold_avedegree_tumor<-data.frame(threshold=seq(0.01,1,0.01),avedegree_num=dcor_tumor_avedegree_num)
  #dcor_threshold_avedegree_normal<-data.frame(threshold=seq(0.01,1,0.01),avedegree_num=dcor_normal_avedegree_num)
  #dcor_threshold_islands_tumor<-data.frame(threshold=seq(0.01,1,0.01),islands_num=dcor_tumor_islands_num)
  #dcor_threshold_islands_normal<-data.frame(threshold=seq(0.01,1,0.01),islands_num=dcor_normal_islands_num)
  #dcor_threshold_transitivity_tumor<-data.frame(threshold=seq(0.01,1,0.01),transitivity_num=dcor_tumor_transitivity_num)
  #dcor_threshold_transitivity_normal<-data.frame(threshold=seq(0.01,1,0.01),transitivity_num=dcor_normal_transitivity_num)
  #dcor_threshold_connectivity_tumor<-data.frame(threshold=seq(0.01,1,0.01),connectivity_num=dcor_tumor_connectivity_num)
  #dcor_threshold_connectivity_normal<-data.frame(threshold=seq(0.01,1,0.01),connectivity_num=dcor_normal_connectivity_num)
  #dcor_threshold_diameter_tumor<-data.frame(threshold=seq(0.01,1,0.01),diameter_num=dcor_tumor_diameter_num)
  #dcor_threshold_diameter_normal<-data.frame(threshold=seq(0.01,1,0.01),diameter_num=dcor_normal_diameter_num)
  value<-cbind(seqs,dcor_tumor_vertex_num,
               dcor_normal_vertex_num,
               dcor_tumor_edges_num,
               dcor_normal_edges_num,
               dcor_tumor_avedegree_num,
               dcor_normal_avedegree_num,
               dcor_tumor_islands_num,
               dcor_normal_islands_num,
               dcor_tumor_transitivity_num,
               dcor_normal_transitivity_num,
               dcor_tumor_connectivity_num,
               dcor_normal_connectivity_num,
               dcor_tumor_diameter_num,
               dcor_normal_diameter_num)
  write.csv(value,paste('network-parameters/',j))
#  png(file=paste(j,'.png'), bg="transparent")
#  #画图
#  plot(dcor_threshold_diameter_tumor,type="p",col=1)
#  points(dcor_threshold_diameter_normal,type="p",col=2)
#  lines(1,lty=1)
#  legend("topright",c("1.tumor","2.normal"),col = c(1,2),lty=1,lwd = 5)
#  title(paste(j,'threshold'))
#  dev.off()
}
#在F:/system_modeling/cors_d/network-parameters/文件夹下面有七个相关性算法得到的七个表
setwd('F:/system_modeling/cors_d/network-parameters/')

corslist1<-dir()
cnames<-c('vertex','edges','averagedegree','islands','transitivity','connectivity','diameter')
cnames2<-c('corgc','dcor','kendall','MIC','multinfo','pearson','spearman')
seqs<-seq(0.01,1,0.01)
j=1
for(k in corslist1){
  parameters<-read.csv(k)[,c(-1,-2)]
  parameters[is.na(parameters)]<-0
  for(i in seq(1,14,2)){
    tumor<-parameters[,i]
    normal<-parameters[,i+1]
    png(file=paste(cnames2[j],cnames[(i+1)/2],'.png',sep=''), bg="white",width=800,height = 600)
    #画图
    plot(x=seqs,y=tumor,type="l",col=1,xlab='阈值',ylab=cnames[(i+1)/2],lwd=2)
    points(x=seqs,y=normal,type="l",col=2,lwd=2)
    legend("topright",c("1.tumor","2.normal"),col = c(1,2),lty=2)
    title(paste(cnames2[j],cnames[(i+1)/2]))
    dev.off()
  }
  j=j+1
}





