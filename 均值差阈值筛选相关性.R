#数据处理
clinical_data<-read.csv('F:/system modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/clinical_data.csv')
clinical_data<-data.frame(sample_id=clinical_data$sampleID,type_id=clinical_data$sample_type_id)
tumor_sample<-clinical_data[which(clinical_data$type_id==1),]
normal_sample<-clinical_data[which(clinical_data$type_id==11),]
tumor_sample_id<-as.character(tumor_sample$sample_id)
tumor_sample_id<-sub('-','.',tumor_sample_id)
tumor_sample_id<-sub('-','.',tumor_sample_id)
tumor_sample_id<-sub('-','.',tumor_sample_id)
normal_sample_id<-as.character(normal_sample$sample_id)
normal_sample_id<-sub('-','.',normal_sample_id)
normal_sample_id<-sub('-','.',normal_sample_id)
normal_sample_id<-sub('-','.',normal_sample_id)


expression_matrixs<-read.csv('F:/system modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/genomicMatrix.csv')
tumor_expression<-expression_matrixs[,colnames(expression_matrixs)%in%tumor_sample_id]
rownames(tumor_expression)<-expression_matrixs[,1]
normal_expression<-expression_matrixs[,colnames(expression_matrixs)%in%normal_sample_id]
rownames(normal_expression)<-expression_matrixs[,1]

write.csv(tumor_expression,'F:/system_modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/tumor_expression_gene.csv')
write.csv(normal_expression,'F:/system modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/normal_expression_gene.csv')

##读入数据
tumor_expression<-read.csv('F:/system_modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/tumor_expression_gene.csv')
normal_expression<-read.csv('F:/system_modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/normal_expression_gene.csv')
rownames(tumor_expression)<-tumor_expression[,1]
rownames(normal_expression)<-normal_expression[,1]
tumor_expression<-tumor_expression[,-1]
normal_expression<-normal_expression[,-1]

#####################
t_normal_expression<-as.data.frame(t(normal_expression))
t_tumor_expression<-as.data.frame(t(tumor_expression))

ttest<-mapply(function(vec1,vec2){
  t.test(x=vec1,y=vec2)$p.value
},vec1=t_tumor_expression,vec2=t_normal_expression)

ttest<-as.numeric(ttest)

#1,求平均值，获取平均值相差4以上的基因
t_normal<-t_normal_expression[,which(ttest<0.001)]
t_tumor<-t_tumor_expression[,which(ttest<0.001)]



a<-lapply(t_tumor,mean)
a<-as.numeric(a)

b<-lapply(t_normal, mean)
b<-as.numeric(b)

difs<-abs(a-b)
t_tumor_gene<-t_tumor[,which(difs>4)]
t_normal_gene<-t_normal[,which(difs>4)]


##求相关系数，用互信息
#library('infotheo')
library('synRNASeqNet')
#d_t_tumor_gene<-discretize(t_tumor_gene)
#d_t_normal_gene<-discretize(t_normal_gene)
counts<-as.matrix(t(t_tumor_gene))
tumor_gene_cor<-parMIEstimate(counts,method='KD',nchips=2)
counts<-as.matrix(t(t_normal_gene))
normal_gene_cor<-parMIEstimate(counts,method='KD',nchips=2)

#pearson get correlations
p_normal_gene_cor<-cor(t_normal_gene,method = 'pearson')
p_tumor_gene_cor<-cor(t_tumor_gene,method = 'spearman')

p_normal_gene_cor[lower.tri(p_normal_gene_cor)]<-NA
diag(p_normal_gene_cor)<-NA
p_tumor_gene_cor[lower.tri(p_tumor_gene_cor)]<-NA
diag(p_tumor_gene_cor)<-NA

library(reshape)

list_tumor_cor<-melt.array(p_tumor_gene_cor)
list_normal_cor<-melt.array(p_normal_gene_cor)
list_tumor_cor<-list_tumor_cor[!is.na(list_tumor_cor$value),]
list_normal_cor<-list_normal_cor[!is.na(list_normal_cor$value),]

list_normal_cor$genepair<-paste(list_normal_cor$X1,list_normal_cor$X2,sep = '-')
list_tumor_cor$genepair<-paste(list_tumor_cor$X1,list_tumor_cor$X2,sep='-')
dotchart(list_normal_cor$value,labels=list_normal_cor$genepair)
library(ggplot2)
ggplot(list_normal_cor,aes(x=factor(list_normal_cor$genepair),y=list_normal_cor$value))+
  geom_dotplot()

#提取相关性阈值大于0.5
library(dplyr)
normal_index=filter(list_normal_cor,list_normal_cor$value>0.5&(list_normal_cor$X1!=list_normal_cor$X2))
tumor_index=filter(list_tumor_cor,list_tumor_cor$value>0.5&(list_tumor_cor$X1!=list_tumor_cor$X2))

write.csv(normal_index,'F:/system_modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/h_normal_index.csv')
write.csv(tumor_index,'F:/system_modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/h_tumor_index.csv')

##########################################################换一种筛选基因方法再做一次。

#p.value阈值曲线
seqs<-seq(1,60,1)
genenum=c()
for(i in seqs) {
  genenum[i]<-ncol(t_tumor_expression[,which(ttest<10^(-i))])
  
}
plot(x=seqs,y=genenum)


#2.获取p.value小于1e-45的基因
t_normal_gene1<-t_normal_expression[,which(ttest<1e-45)]
t_tumor_gene1<-t_tumor_expression[,which(ttest<1e-45)]

##求相关系数，用互信息
library('synRNASeqNet')
counts<-as.matrix(t(t_tumor_gene1))
tumor_gene_cor1<-parMIEstimate(counts,method='KD',nchips=2)
counts<-as.matrix(t(t_normal_gene1))
normal_gene_cor1<-parMIEstimate(counts,method='KD',nchips=2)

library(reshape)
list_tumor_cor1<-melt.array(tumor_gene_cor1)
list_normal_cor1<-melt.array(normal_gene_cor1)

#提取相关性阈值大于0.5
library(dplyr)
normal_index1=filter(list_normal_cor1,list_normal_cor1$value>0.5&(list_normal_cor1$X1!=list_normal_cor1$X2))
tumor_index1=filter(list_tumor_cor1,list_tumor_cor1$value>0.5&(list_tumor_cor1$X1!=list_tumor_cor1$X2))

write.csv(normal_index1,'F:/system_modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/h_normal_index1.csv')
write.csv(tumor_index1,'F:/system_modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/h_tumor_index1.csv')

#获取网络参数
normal_index<-read.csv('F:/allcors/allcors/pearson_all_cors.csv')
library(igraph)
normal_index1<-normal_index[,c(2,3,4)]
normal_graph<-graph.data.frame(normal_index1,directed = FALSE)
tumor_graph<-graph.data.frame(tumor_index,directed = FALSE)

short_tumor<-shortest.paths(tumor_graph)
short_normal<-shortest.paths(normal_graph)

write.csv(short_tumor,"F:/system_modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/short_tumor.csv")
write.csv(short_normal,"F:/system_modeling/TCGA_COAD_exp_HiSeqV2-2015-02-24/short_normal.csv")

# No of nodes
length(V(normal_graph))
length(V(tumor_graph))
# No of edges
length(E(normal_graph))
length(E(tumor_graph))
# Density (No of edges / possible edges)
graph.density(normal_graph)
graph.density(tumor_graph)
# Number of islands
clusters(normal_graph)$no
clusters(tumor_graph)$no
# Global cluster coefficient:
#(close triplets/all triplets)
transitivity(normal_graph, type="global")
transitivity(tumor_graph,type = 'global')
# Edge connectivity, 0 since graph is disconnected
edge.connectivity(normal_graph)
edge.connectivity(tumor_graph)
# Diameter of the graph
diameter(normal_graph)
diameter(tumor_graph)

deg_normal<-degree.distribution(normal_graph)
deg_tumor<-degree.distribution(tumor_graph)

plot(deg_tumor, xlab="node degree")
plot(deg_normal, xlab="node degree")

lines(deg_tumor)
lines(deg_normal)







