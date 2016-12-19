#一读入sp数据
sp_order<-read.csv("C:\\Users\\Administrator\\Desktop\\基因网络\\Data\\cors\\sp_all_cors.csv")
#正常组归一化处理 
max<-1
min<--1
for(i in 1:length(sp_order[,4])){
  sp_order[i,4]<-(sp_order[i,4]-min)/(max-min)
}
#癌症组归一化处理 
for(i in 1:length(sp_order[,5])){
  sp_order[i,5]<-(sp_order[i,5]-min)/(max-min)
}
#将正常组和癌症组的数据分开了
sp_normal<-sp_order[,c(2,3,4)]
sp_tumor<-sp_order[,c(2,3,5)]
#按照相关性值排序
sp_normal<-sp_normal[order(sp_normal[,3],decreasing = T),]
sp_tumor<-sp_tumor[order(sp_tumor[,3],decreasing = T),]

sp_normal<-data.frame(pair=seq(1,44850,1),value=sp_normal$normal)
sp_tumor<-data.frame(pair=seq(1,44850,1),value=sp_tumor$tumor)

#二读入pearson数据
pearson_order<-read.csv("C:\\Users\\Administrator\\Desktop\\基因网络\\Data\\cors\\pearson_all_cors.csv")
#正常组归一化处理 
max<-1
min<--1
for(i in 1:length(pearson_order[,4])){
  pearson_order[i,4]<-(pearson_order[i,4]-min)/(max-min)
}
#癌症组归一化处理 
for(i in 1:length(pearson_order[,5])){
  pearson_order[i,5]<-(pearson_order[i,5]-min)/(max-min)
}
#将正常组和癌症组的数据分开了
pearson_normal<-pearson_order[,c(2,3,4)]
pearson_tumor<-pearson_order[,c(2,3,5)]
#按照相关性值排序
pearson_normal<-pearson_normal[order(pearson_normal[,3],decreasing = T),]
pearson_tumor<-pearson_tumor[order(pearson_tumor[,3],decreasing = T),]

pearson_normal<-data.frame(pair=seq(1,44850,1),value=pearson_normal$normal)
pearson_tumor<-data.frame(pair=seq(1,44850,1),value=pearson_tumor$tumor)

#三读dcor入数据
dcor_order<-read.csv("C:\\Users\\Administrator\\Desktop\\基因网络\\Data\\cors\\dcor_all_cors.csv")
#正常组归一化处理 
max<-1
min<-0
for(i in 1:length(dcor_order[,4])){
  dcor_order[i,4]<-(dcor_order[i,4]-min)/(max-min)
}
#癌症组归一化处理 
for(i in 1:length(dcor_order[,5])){
  dcor_order[i,5]<-(dcor_order[i,5]-min)/(max-min)
}
#将正常组和癌症组的数据分开了
dcor_normal<-dcor_order[,c(2,3,4)]
dcor_tumor<-dcor_order[,c(2,3,5)]
#按照相关性值排序
dcor_normal<-dcor_normal[order(dcor_normal[,3],decreasing = T),]
dcor_tumor<-dcor_tumor[order(dcor_tumor[,3],decreasing = T),]

dcor_normal<-data.frame(pair=seq(1,44850,1),value=dcor_normal$normal)
dcor_tumor<-data.frame(pair=seq(1,44850,1),value=dcor_tumor$tumor)


#四读入kendall数据
ke_order<-read.csv("C:\\Users\\Administrator\\Desktop\\基因网络\\Data\\cors\\ke_all_cors.csv")
#正常组归一化处理 
max<-1
min<--1
for(i in 1:length(ke_order[,4])){
  ke_order[i,4]<-(ke_order[i,4]-min)/(max-min)
}
#癌症组归一化处理 
for(i in 1:length(ke_order[,5])){
  ke_order[i,5]<-(ke_order[i,5]-min)/(max-min)
}
#将正常组和癌症组的数据分开了
ke_normal<-ke_order[,c(2,3,4)]
ke_tumor<-ke_order[,c(2,3,5)]
#按照相关性值排序
ke_normal<-ke_normal[order(ke_normal[,3],decreasing = T),]
ke_tumor<-ke_tumor[order(ke_tumor[,3],decreasing = T),]

ke_normal<-data.frame(pair=seq(1,44850,1),value=ke_normal$normal)
ke_tumor<-data.frame(pair=seq(1,44850,1),value=ke_tumor$tumor)

#五读入mic数据
mic_order<-read.csv("C:\\Users\\Administrator\\Desktop\\基因网络\\Data\\cors\\mic_all_cors.csv")
#正常组归一化处理 
max<-1
min<-0
for(i in 1:length(mic_order[,4])){
  mic_order[i,4]<-(mic_order[i,4]-min)/(max-min)
}
#癌症组归一化处理 
for(i in 1:length(mic_order[,5])){
  mic_order[i,5]<-(mic_order[i,5]-min)/(max-min)
}
#将正常组和癌症组的数据分开了
mic_normal<-mic_order[,c(2,3,4)]
mic_tumor<-mic_order[,c(2,3,5)]
#按照相关性值排序
mic_normal<-mic_normal[order(mic_normal[,3],decreasing = T),]
mic_tumor<-mic_tumor[order(mic_tumor[,3],decreasing = T),]

mic_normal<-data.frame(pair=seq(1,44850,1),value=mic_normal$normal)
mic_tumor<-data.frame(pair=seq(1,44850,1),value=mic_tumor$tumor)


#六读入multinfo数据
multinfo_order<-read.csv("C:\\Users\\Administrator\\Desktop\\基因网络\\Data\\cors\\multinfo_all_cors.csv")
#正常组归一化处理 
max<-5
min<--5
for(i in 1:length(multinfo_order[,4])){
  multinfo_order[i,4]<-(multinfo_order[i,4]-min)/(max-min)
}
#癌症组归一化处理 
for(i in 1:length(multinfo_order[,5])){
  multinfo_order[i,5]<-(multinfo_order[i,5]-min)/(max-min)
}
#将正常组和癌症组的数据分开了
multinfo_normal<-multinfo_order[,c(2,3,4)]
multinfo_tumor<-multinfo_order[,c(2,3,5)]
#按照相关性值排序
multinfo_normal<-multinfo_normal[order(multinfo_normal[,3],decreasing = T),]
multinfo_tumor<-multinfo_tumor[order(multinfo_tumor[,3],decreasing = T),]

multinfo_normal<-data.frame(pair=seq(1,44850,1),value=multinfo_normal$normal)
multinfo_tumor<-data.frame(pair=seq(1,44850,1),value=multinfo_tumor$tumor)

#七读入corgc数据
corgc_order<-read.csv("C:\\Users\\Administrator\\Desktop\\基因网络\\Data\\cors\\corgc_all_cors.csv")
#正常组归一化处理 
max<-1
min<-0
for(i in 1:length(corgc_order[,4])){
  corgc_order[i,4]<-(corgc_order[i,4]-min)/(max-min)
}
#癌症组归一化处理 
for(i in 1:length(corgc_order[,5])){
  corgc_order[i,5]<-(corgc_order[i,5]-min)/(max-min)
}
#将正常组和癌症组的数据分开了
corgc_normal<-corgc_order[,c(2,3,4)]
corgc_tumor<-corgc_order[,c(2,3,5)]
#按照相关性值排序
corgc_normal<-corgc_normal[order(corgc_normal[,3],decreasing = T),]
corgc_tumor<-corgc_tumor[order(corgc_tumor[,3],decreasing = T),]

corgc_normal<-data.frame(pair=seq(1,40709,1),value=corgc_normal$normal)
corgc_tumor<-data.frame(pair=seq(1,40709,1),value=corgc_tumor$tumor)



#tumor组画图
plot(sp_tumor,type="p",main="tumor组",xlab="基因对排序",ylab="相关性",col=1)
#points用来向plot上添加曲线，实质用法和plot一样
points(pearson_tumor,type="p",col=2)
points(ke_tumor,type="p",col=3)
points(dcor_tumor,type="p",col=4)
points(mic_tumor,type="p",col=5)
points(multinfo_tumor,type="p",col=6)
points(corgc_tumor,type="p",col=8)
#设置图的样例
lines(1,lty=1)
cor_method<-c("1.sp","2.pearson","3.ke","4.doc","5.mic","6.multinfo","7.corgc")
legend("topright",cor_method,col = c(1,2,3,4,5,6,8),lty=1,lwd = 5,cex=0.7)

#normal组画图
plot(sp_normal,type="p",main="normal组",xlab="基因对排序",ylab="相关性",col=1)
#points用来向plot上添加曲线，实质用法和plot一样
points(pearson_normal,type="p",col=2)
points(ke_normal,type="p",col=3)
points(dcor_normal,type="p",col=4)
points(mic_normal,type="p",col=5)
points(multinfo_normal,type="p",col=6)
points(corgc_normal,type="p",col=8)
#设置图的样例
lines(1,lty=1)
cor_method<-c("1.sp","2.pearson","3.ke","4.doc","5.mic","6.multinfo","7.corgc")
legend("topright",cor_method,col = c(1,2,3,4,5,6,8),lty=1,lwd = 5,cex=0.7)





