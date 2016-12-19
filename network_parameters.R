library(igraph)
s_normal_gene<-read.csv("F:/system modeling/normal_index.csv")[,-1]
s_tumor_gene<-read.csv('F:/system modeling/tumor_index.csv')[,-1]
s_normal_graph<-graph.data.frame(s_normal_gene,directed = FALSE)
s_tumor_graph<-graph.data.frame(s_tumor_gene,directed = FALSE)
g<-s_tumor_graph
plot(s_tumor_graph)

#find all shortest path
short_tumor<-shortest.paths(s_tumor_graph)
short_normal<-shortest.paths(s_normal_graph)

write.csv(short_tumor,"F:/system modeling/TCGA_COADREAD_exp_HiSeqV2-2015-02-24/short_tumor.csv")
write.csv(short_normal,"F:/system modeling/TCGA_COADREAD_exp_HiSeqV2-2015-02-24/short_normal.csv")

#compute some measures

# No of nodes
len_v_normal<-length(V(s_normal_graph))
len_v_tumor<-length(V(s_tumor_graph))

# No of edges
len_e_normal<-length(E(s_normal_graph))
len_e_tumor<-length(E(s_tumor_graph))

# Density (No of edges / possible edges)
den_normal<-graph.density(s_normal_graph)
den_tumor<-graph.density(s_tumor_graph)

# Number of islands
clus_normal<-clusters(s_normal_graph)$no
clus_tumor<-clusters(s_tumor_graph)$no

# Global cluster coefficient:
#(close triplets/all triplets)
trans_normal<-transitivity(s_normal_graph, type="global")
trans_tumor<-transitivity(s_tumor_graph,type = 'global')

# Edge connectivity, 0 since graph is disconnected
connc_normal<-edge.connectivity(s_normal_graph)
connc_tumor<-edge.connectivity(s_tumor_graph)

# Diameter of the graph
diam_normal<-diameter(s_normal_graph)
diam_tumor<-diameter(s_tumor_graph)
# distribution

deg_normal<-degree.distribution(s_normal_graph)
deg_tumor<-degree.distribution(s_tumor_graph)

plot(deg_tumor, xlab="node degree")
plot(deg_normal, xlab="node degree")

lines(deg_tumor)
lines(deg_normal)








