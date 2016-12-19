t_normal_expression<-read.csv('C:\\Users\\Administrator\\Desktop\\基因网络\\Data\\t_normal_expression_gene.csv')
t_tumor_expression<-read.csv('C:\\Users\\Administrator\\Desktop\\基因网络\\Data\\t_tumor_expression_gene.csv')

index <- c()
library('car')
for (i in 2:ncol(t_normal_expression)){
  #先判断方差齐不齐
  m <- as.vector(t_tumor_expression[,i]) 
  n <- as.vector(t_normal_expression[,i])
  y<-c(m,n)
  #y <- as.vector(cbind(m,n))
  a <- c(rep(1,length(m)));
  b <- c(rep(2,length(n)))
  group<-as.factor(c(a,b))
  #group <- as.factor(cbind(a,b))
  
  if(leveneTest(y = y, group=group)$Pr[1] >0.05){#方差齐(有11664种基因)
    p <- t.test(m,n,var.equal = TRUE)$p.value
  }else{
    p <- t.test(m,n,var.equal = FALSE)$p.value#方差不齐
  }
  
  if (p>=0.01)
    index <- c(index,i)
}
filter_tumor <- t_tumor_expression[,-index]
filter_normal <- t_normal_expression[,-index]