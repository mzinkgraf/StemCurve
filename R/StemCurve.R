#calculate bending for timelaps vidoes

directory<-"Data"

source("R/functions.R")
require(plyr)
require(ggplot2)

#get list of directories
dir<-dir(directory,pattern="717|35s37|miRNA|T2")

#build table of metadata
meta<-data.frame(matrix(vector(), length(dir), 5, dimnames=list(c(), c("dir","genotype","tree","data","treatment"))), stringsAsFactors=F)
meta[,1]<-dir
meta[2:5]<-data.frame(matrix(unlist(strsplit(dir,"_")),nrow=length(dir),byrow=T))

#load distance to hook
hook<-read.table("Data/distance.txt",sep="\t",header=T,stringsAsFactors = F)
meta<-merge(meta,hook,by.x="dir",by.y="dir")


#loop through directories and get data


results<-list()
grep_key<-".txt"
for(e in 1:nrow(meta))
{
  #check that each directory has been fully processes (.jpg == .txt)
  all_files<-list.files(paste(directory,meta[e,1],sep="/"))
  nt<-length(grep(".txt",all_files,perl=TRUE))
  nj<-length(grep(".JPG|jpg",all_files,perl=TRUE))
  if(nt!=nj) {print(paste("Not all files in ",meta[e,1]," have been processed",sep="")) } else
  {
    results[[paste(meta[e,1])]]<-StemCurve(paste(directory,"/",meta[e,1],sep=""),grep_key,meta[e,6],meta[e,7])
  } 
  
}

 pdf(file="Data/Bending_analysis.pdf")
 for(q in 1:length(results))
 {
 p<-ggplot(results[[q]]$data,aes(x=x,y=y,colour=file)) + geom_point(size=4) + theme_bw() + theme(legend.position="none") + ggtitle(names(results[q])) + theme(plot.title = element_text(size=28, face="bold")) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18))  
    print(p)
  }
  dev.off()

#717_127_06.06.13_H2O
#last image missing add missing data for:
results[['717_127_06.06.13_H2O']]$stat[15,1:2]<-c(NA,NA)

bending_stats<-data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("dir","curv","height","time"))), stringsAsFactors=F)

for(b in 1:length(results))
{
  temp<-cbind(names(results)[b],results[[b]]$stat,c(0:14))
  names(temp)<-c("dir","curv","height","time")
  bending_stats<-rbind(bending_stats,temp)
}

mbending_stats<-merge(bending_stats,meta,by.x="dir",by.y="dir",all.x=T)

mbending_stats$treatment<-factor(mbending_stats$treatment, levels=c("H2O"))


#calculate rate of bending

for(u in 1:nrow(mbending_stats))
{
  if(mbending_stats$time[u]==0)
  {
    mbending_stats[u,10]=0
    mbending_stats[u,11]=0
    zr<-mbending_stats[u,3]
  } else {
    mbending_stats[u,10]<-mbending_stats[u,3]-mbending_stats[u-1,3]
    #mbending_stats[u,11]<-mbending_stats[u,10]+mbending_stats[u-1,11]
    mbending_stats[u,11]<-mbending_stats[u,3]-zr
  }
}
names(mbending_stats)[10:11]<-c("rate","cum_height")


boxplot(cum_height~ time, data=mbending_stats)

model1<-glm(cum_height~ time , data=mbending_stats)
summary(model1)


su<-ddply(mbending_stats, .(genotype,time), summarize, 
          height_mean=mean(height,na.rm=T),
          height_se=sd(height,na.rm=T)/sqrt(length(which(!is.na(height)))),
          hn=length(which(!is.na(height))),
          curv_mean=mean(curv,na.rm=T),
          curv_se=sd(curv,na.rm=T)/sqrt(length(which(!is.na(curv)))),
          cn=length(which(!is.na(curv))),
          rate_mean=mean(rate,na.rm=T),
          rate_se=sd(rate,na.rm=T)/sqrt(length(which(!is.na(rate)))),
          cum_mean=mean(cum_height,na.rm=T),
          cum_se=sd(cum_height,na.rm=T)/sqrt(length(which(!is.na(cum_height)))))


#pdf(file="Stem_curve.pdf",width=8,height=6, useDingbats = F)
ggplot(su,aes(x=time,curv_mean,colour=genotype))+geom_point(size=3) + theme_bw() + geom_errorbar(aes(ymin=curv_mean-curv_se,ymax=curv_mean+curv_se,width=0.2)) +xlab("Time (days)") + ylab("Stem Curvature (degree) \n") + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.text=element_text(size=16), legend.title=element_text(size=16))
#dev.off()

#pdf(file="Stem_rate.pdf",width=8,height=6, , useDingbats = F)
ggplot(su[which(su$time>0),],aes(x=time,rate_mean,colour=genotype))+geom_point(size=3) + geom_errorbar(aes(ymin=rate_mean-rate_se,ymax=rate_mean+rate_se,width=0.2)) + theme_bw()  +xlab("Time (days)") + ylab("Rate of Bending\n(lift / day / stem length)\n") + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.text=element_text(size=16), legend.title=element_text(size=16))
#dev.off()

#pdf(file="Stem_height.pdf",width=8,height=6)
ggplot(su,aes(x=time,cum_mean,colour=genotype))+geom_point(size=2) + geom_errorbar(aes(ymin=cum_mean-cum_se,ymax=cum_mean+cum_se,width=0.2)) + theme_bw()  +xlab("Time (days)") + ylab("Normalized Stem Lift\n\n(lift / stem length)\n") + theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.text=element_text(size=16), legend.title=element_text(size=16))
#dev.off()
