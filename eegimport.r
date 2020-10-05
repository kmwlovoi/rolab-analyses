eegimport<-function(behavfn,headerlines,responsetype,eegfn,eventfn){
#///////////////////////////////////////////
#// Version: 1.0.0                    	   //
#// Date: 04/01						  	      //
#// Desc: Imports all behavioral/EEG data //
#//		  for PhosTOJ study 	            //
#// behavfn= behav file; headerlines= # of//
#// header lines in behav file; response  //
#// type= response instructions, 1: (L= PF//
#// R= VF), 0: (L= VF, R= PF); eegfn= eeg //
#// file; eventfn= event file			      //
#///////////////////////////////////////////
#import behavioral data
SubjDat<- read.table(behavfn,sep="\t",skip=headerlines,header=FALSE) 

colnames(SubjDat) <- c("subj","dur","soa","fixtime","vistime","tmstime","rt","resp")
	


#import EEG data
srate=1000
nchannels=32
bytespersample= 2

EEGfileinfo= file.info(eegfn)
EEGdata=readBin(eegfn,'int',n=(EEGfileinfo$size/bytespersample),size=bytespersample,endian="little")

datapoints= EEGfileinfo$size/(nchannels*bytespersample)

dim(EEGdata)<-c(nchannels,datapoints) #reshape the data so that rows=channels and columns=timepoints
EEGdata=EEGdata*.1 #scale the data


# #import event info
events<- read.table(eventfn,sep=",",skip=12) #imports the trigger information from the vmrk file and skips the first 11 header lines and the first event (which is a boundary event)
events$V2<- sapply(events$V2,gsub,pattern="S  ",replacement="")
events$V2<- sapply(events$V2,gsub,pattern="S ",replacement="")

if(eventfn=='PhosTOJ/EEG_0002.vmrk')
	events <- events[c(1:180),]

return(list(SubjDat,EEGdata,events))

}

