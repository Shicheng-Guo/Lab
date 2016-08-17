setwd("C:\\Users\\shicheng\\Documents\\GitHub\\collaboration\\blue")

barplotReAxis<-function(mp){
  mp<-c(0,mp)
  np<-c()
  end<-length(mp)-1
  for(i in 1:end){
    np.tmp<-(mp[i]+mp[i+1])/2
    np<-c(np,np.tmp)
  }
  np<-round(np)
  names(np)=names(mp)[2:length(mp)]
  return(np)
}


nbt.data<-data.matrix(read.table("BrainRNASeqBlue.txt",head=T,row.names=1,sep="\t",as.is=T))

nbt.data[1:5,1:5]

i=6
rmpk<-nbt.data[i,]
genename<-rownames(nbt.data)[i]
genename

cell.labels <- substr(names(rmpk),0,7)
cell.labels<-unlist(lapply(cell.labels,function(x) unlist(strsplit(x,"_"))[1]))
level=c(paste("In",1:8,sep=""),paste("Ex",1:8,sep=""))
group1<-factor(unlist(lapply(cell.labels,function(x) substr(x,1,3))))
group2<-factor(unlist(lapply(cell.labels,function(x) substr(x,4,nchar(x)))))
group3<-1:length(group2)
df<-data.frame(rmpk,cell.labels,group1,group2,group3)
df<-df[order(df$cell.labels),]
head(df)

filename=paste(genename,".png",sep="")
# png(filename,width = 8, height = 2, units = 'in', res = 300)
par(mar=c(2, 3, 1, 3), xpd=TRUE)
space<-rep(0,nrow(df))
space[cumsum(table(df$cell.labels))+1]<-30
space[cumsum(table(df$group1))+1]<-200
space<-space[1:(length(space)-1)]
bar.positions <- barplot(df$rmpk,border=df$group2,xaxt='n',space=space,main=genename,cex.main=0.6,cex.lab=0.4,cex.axis=0.4)
title(ylab = "Counts per gene", cex.lab = 0.4,line = 2)
label<-levels(df$group1)
table(df$group1)
pos<-bar.positions[cumsum(barplotReAxis(table(df$group1)))]
axis(side=1,at=pos,labels=label,tick=FALSE,las=2,cex.axis=0.4)
legend(x=max(bar.positions)+100,y=max(rmpk)/2,legend=levels(df$group2),bty="n",lty=1,bg="transparent",col=c(as.numeric(unique(df$group2))),cex=0.4,inset=c(-0.2,0))
# dev.off()

