library(igraph,warn.conflicts = F)
library(infotheo)
library(dplyr)
#设置当前的工作目录
setwd('F:/allcors/allcors')
#获取当前目录下面的所有文件
corslist<-dir()
for(j in corslist){
  #j的值为文件的名称
  dcor_order<-read.csv(j)
  normal<-dcor_order[,c(2,3,4)]
  tumor<-dcor_order[,c(2,3,5)]
  
  seqs<-seq(0.01,1,0.01)
  #平均路径长度L
  dcor_tumor_l_num<-c()
  dcor_normal_l_num<-c()
  #平均度
  dcor_tumor_avedegree_num<-c()
  dcor_normal_avedegree_num<-c()
  #平均边介数betweenness
  dcor_tumor_betweenness_num<-c()
  dcor_normal_betweenness_num<-c()
  #聚类系数transitivity,表示图中节点的聚集程度
  dcor_tumor_transitivity_num<-c()
  dcor_normal_transitivity_num<-c()
  #模块度Q
  dcor_tumor_modularity_num<-c()
  dcor_normal_modularity_num<-c()
  #度分布degree distribution
  dcor_tumor_degree_num<-c()
  dcor_normal_degree_num<-c()
  #平均核数average coreness
  dcor_tumor_corness_num<-c()
  dcor_normal_corness_num<-c()
  for(i in 1:100){
    #进行阈值筛选
    filter_tumor_dcor_order=tumor[tumor$tumor>=(i*0.01),]
    filter_normal_dcor_order=normal[normal$normal>=(i*0.01),]

    graph_tumor_dcor_order<-graph.data.frame(filter_tumor_dcor_order,directed = FALSE)
    graph_normal_dcor_order<-graph.data.frame(filter_normal_dcor_order,directed =FALSE)
    
    #平均路径长度L
    dcor_tumor_l_num[i]<-average.path.length(graph_tumor_dcor_order)
    dcor_normal_l_num[i]<-average.path.length(graph_normal_dcor_order)
    
    #平均边介数
    #先求出每条边的介数
    Betweenness_tumor<-edge_betweenness(graph_tumor_dcor_order,e=E(graph_tumor_dcor_order),directed = F)
    Betweenness_normal<-edge_betweenness(graph_normal_dcor_order,e=E(graph_normal_dcor_order),directed = F)
    #所有边的介数之和再除以边数
    dcor_tumor_betweenness_num[i]<-sum(Betweenness_tumor)/length(E(graph_tumor_dcor_order))
    dcor_normal_betweenness_num[i]<-sum(Betweenness_normal)/length(E(graph_normal_dcor_order))
    
    # average degree
    #先计算每个顶点的度
    d_tumor<-degree(graph_tumor_dcor_order)
    d_normal<-degree(graph_normal_dcor_order)
    #所有顶点的度数之和再除以顶点数
    dcor_tumor_avedegree_num[i]<-sum(d_tumor)/length(V(graph_tumor_dcor_order))
    dcor_normal_avedegree_num[i]<-sum(d_normal)/length(V(graph_normal_dcor_order))
    

    #transitivity
    dcor_tumor_transitivity_num[i]<-transitivity(graph_tumor_dcor_order, type="global")
    dcor_normal_transitivity_num[i]<-transitivity(graph_normal_dcor_order, type="global")
    
    #模块度 modularity 
    #先计算出社群模型，使用随机游走方法
    wtc_tumor<-cluster_walktrap(graph_tumor_dcor_order)
    wtc_normal<-cluster_walktrap(graph_normal_dcor_order)
    #模块化指标Q
    dcor_tumor_modularity_num[i]<-modularity(graph_tumor_dcor_order,membership (wtc_tumor))
    dcor_normal_modularity_num[i]<-modularity(graph_normal_dcor_order,membership(wtc_normal))
    
    #先求所有点的核数
    corness_tumor<-coreness(graph_tumor_dcor_order)
    corness_normal<-coreness(graph_normal_dcor_order)
    #平均核数 所有点的核数之和再除以顶点数
    dcor_tumor_corness_num[i]<-sum(corness_tumor)/length(V(graph_tumor_dcor_order))
    dcor_normal_corness_num[i]<-sum(corness_normal)/length(V(graph_normal_dcor_order))
   
    #     #number of islands
    #     dcor_tumor_islands_num[i]<-clusters(graph_tumor_dcor_order)$no
    #     dcor_normal_islands_num[i]<-clusters(graph_normal_dcor_order)$no
    
  }
  
  number<-data.frame(threshold=seqs,
                     average_path_length_tumor=dcor_tumor_l_num,
                     average_path_length_normal=dcor_normal_l_num,
                     betweenness_tumor=dcor_tumor_betweenness_num,
                     betweenness_normal=dcor_normal_betweenness_num,
                     avedegree_tumor=dcor_tumor_avedegree_num,
                     avedegree_normal=dcor_normal_avedegree_num,
                     transitivity_tumor=dcor_tumor_transitivity_num,
                     transitivity_normal=dcor_normal_transitivity_num,
                     modularity_tumor=dcor_tumor_modularity_num,
                     modularity_normal=dcor_normal_modularity_num,
                     coreness_tumor=dcor_tumor_corness_num,
                     coreness_normal=dcor_normal_corness_num
                     )
  write.csv(number,paste("F:/network_parameter/",j))
  
}

seqs<-seq(0.01,1,0.01)
cor_name<-c("corgc","dcor","kendall","mic","multinfo","pearson","spearman")
parameter_name<-c("average_path_length","betweenness","avedegree","transitivity","modularity","coreness")
setwd('F:/network_parameter')
corlist1<-dir()
#i表示相关性算法名称
k=1
for(i in corlist1){
  number<-read.csv(i)
  number[is.na(number)]<-0
  number<-number[,c(-1,-2)]
  #j表示结构参数
  for(j in seq(1,11,2)){
    tumor<-number[,j]
    normal<-number[,j+1]
    png(file=paste(cor_name[k],parameter_name[(j+1)/2],'.png',sep =""), bg="transparent",width = 800,height = 600)
    plot(x=seqs,y=tumor,type="l",xlab="阈值",ylab=paste(parameter_name[(j+1)/2],"number"),col=1)
    points(x=seqs,y=normal,type="l",col=2)
    lines(1,lty=1)
    legend("topright",c("1.tumor","2.normal"),col = c(1,2),lty=1,lwd = 5)
    title(paste(cor_name[k],parameter_name[(j+1)/2]))
    dev.off()
  }
  k=k+1
}

  



