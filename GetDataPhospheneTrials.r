GetDataPhospheneTrials<-function(fn,headerlines,responsetype){
#Same as GetData script, but analyses only trials on which a phosphene was detected
   if(responsetype==1){
		phosF=16
		visF=32
   } else{
		phosF=32
		visF=16
   }

SubjDat<- read.table(fn,sep="\t",skip=headerlines) 

SubjDat[which(SubjDat[,8]==phosF),9]=1
SubjDat[which(SubjDat[,8]==visF),9]=0

phosproportions <- array(NA,dim=c(2,11))

for(i in 0:9) {
	phosproportions[1,i+1] = sum(SubjDat[which(SubjDat$V3==i),9],na.rm=TRUE)/length(which(SubjDat$V3==i&SubjDat$V8!='X')) 
   phosproportions[2,i+1] = length(which(SubjDat$V3==i&!is.na(SubjDat$V9)))
}
phosproportions[1,11] = length(which(SubjDat$V3==10&SubjDat$V8=='X'))/length(which(SubjDat$V3==10))
phosproportions[2,11] = length(which(SubjDat$V3==10&SubjDat$V8=='X'))

if(exists("allsubj")){
   rownames(phosproportions)<-c(paste(((nrow(allsubj)/2)+1),"phosprop",sep=""),paste(((nrow(allsubj)/2)+1),"TotalPhos",sep=""))
   allsubj = rbind(allsubj,phosproportions)
} else {
   rownames(phosproportions)<-c("phosprop","TotalPhos")
   allsubj = phosproportions
}

return(allsubj)
}