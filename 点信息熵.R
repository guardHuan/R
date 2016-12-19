#计算点信息熵
vertex_entropy <- function(pdnor){
  allvertex <- union(pdnor$Var1,pdnor$Var2)#并集,得到所有顶点的集合
  vertexH <- 0#点信息熵
  for (eachvertex in allvertex){
    library(dplyr)
    vertex_p<-0.5*(sum((filter(pdnor, pdnor$tumor_name==eachvertex|pdnor$normal_name==eachvertex))$correlation))#每个顶点的概率
    if(vertex_p!=0)
      vertexH <- vertexH+vertex_p*log(vertex_p)
  }
  return(-vertexH)
}


d_vertexH_1 <-vertex_entropy(d_pdnor1)
d_vertexH_2 <-vertex_entropy(d_pdnor2)
d_vertexH_3 <-vertex_entropy(d_pdnor3)
d_vertexH_5 <-vertex_entropy(d_pdnor5)
d_vertexH_6 <-vertex_entropy(d_pdnor6)
d_vertexH_7 <-vertex_entropy(d_pdnor7)
d_vertexH_8 <-vertex_entropy(d_pdnor8)
d_allvertexH <- c(d_vertexH_1,d_vertexH_2,d_vertexH_3,d_vertexH_5,d_vertexH_6,d_vertexH_7,d_vertexH_8)
#绘制疾病组7种点信息熵的条形图
#barplot(d_allvertexH , main="疾病组7种点信息熵的条形图",xlab="method", ylab="点信息熵", 
#names.arg=c("pearson","spearman","Kendall","互信息","CorGc","dCor","MIC"),
#col=c(1,2,3,5,6,7,8))
h_vertexH_1 <-vertex_entropy(h_pdnor1)
h_vertexH_2 <-vertex_entropy(h_pdnor2)
h_vertexH_3 <-vertex_entropy(h_pdnor3)
h_vertexH_5 <-vertex_entropy(h_pdnor5)
h_vertexH_6 <-vertex_entropy(h_pdnor6)
h_vertexH_7 <-vertex_entropy(h_pdnor7)
h_vertexH_8 <-vertex_entropy(h_pdnor8)
h_allvertexH <- c(h_vertexH_1,h_vertexH_2,h_vertexH_3,h_vertexH_5,h_vertexH_6,h_vertexH_7,h_vertexH_8)

#绘制正常组7种点信息熵的条形图
barplot(h_allvertexH , main="正常组7种点信息熵的条形图",xlab="method", ylab="点信息熵", 
        names.arg=c("pearson","spearman","Kendall","互信息","CorGc","dCor","MIC"),
        col=c(1,2,3,5,6,7,8))


#折线图
allH <- union(d_allvertexH,h_allvertexH)
plot(c(1,2,3,5,6,7,8),d_allvertexH,type='b',
     xlab="method", ylab="点信息熵",col=2,xlim=c(1,8),ylim=c(min(allH),max(allH)))
points(c(1,2,3,5,6,7,8),h_allvertexH,type='b',
       xlab="method", ylab="点信息熵",col=9,xlim=c(1,8),ylim=c(min(allH),max(allH)))
lines(1,lty=1)
cor_method <- c("疾病组","正常组")
legend("topright",cor_method,col=c(2,9),lty=1,lwd=5,x.intersp=0.1,y.intersp=0.3,)
