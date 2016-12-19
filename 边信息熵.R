setwd('F:/allcors/allcors/')
corslist<-dir()

#pdnor为读取的数据文件，返回这个数据文件中每个顶点的信息熵之和
vertex_entropy <- function(pdnor){
  allvertex <- union(pdnor$tumor_name,pdnor$normal_name)#并集,得到所有顶点的集合
  vertexH <- 0#点信息熵
  for (eachvertex in allvertex){
    library(dplyr)
    vertex_p<-0.5*(sum((filter(pdnor, pdnor$tumor_name==eachvertex|pdnor$normal_name==eachvertex))[,3]))#每个顶点的概率
    if(vertex_p!=0)
      vertexH <- vertexH+vertex_p*log(vertex_p)
  }
  return(-vertexH)
}



h1<-c()
h2<-c()
v1<-c()
v2<-c()
for(j in corslist){
  data<-read.csv(j)
  normal<-abs(data[,4])
  tumor<-abs(data[,5])
  v_normal<-cbind(data[,2:3],normal/sum(normal))
  v_tumor<-cbind(data[,2:3],tumor/sum(tumor))
  
  v1<-c(v1,vertex_entropy(v_normal))
  v2<-c(v2,vertex_entropy(v_tumor))
  
  
  #bian xin xi shang
  
  p_normal<-normal/sum(normal)
  p_tumor<-tumor/sum(tumor)
  p_normal<-p_normal[p_normal!=0]
  p_tumor<-p_tumor[p_tumor!=0]
  h1<-c(h1,sum(as.numeric(lapply(p_normal, function(p){-1*p*log2(p)}))))
  h2<-c(h2,sum(as.numeric(lapply(p_tumor, function(p){-1*p*log2(p)}))))
  
  #dian xin xi shang
}


png(file='bian.png', bg="transparent")
e_e<-cbind(h1,h2,v1,v2)
rownames(e_e)<-corslist
plot(e_e[,1],col=1,type='l')
points(e_e[,2],col=2,type='l')
title('edges_entropy')
lines(2,lty=1)
cor_method<-c("1.normal","2.tumor")
legend("bottomleft",cor_method,col = c(1,2),lty=1,lwd = 5,cex=0.7)
dev.off()
write.csv(e_e,'entropy.csv')






























































































