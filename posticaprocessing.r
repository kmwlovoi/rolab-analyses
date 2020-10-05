posticaprocessing <- function (eeg,vis,adjust,behav,gc) {
#///////////////////////////////////////////
#// Version: 1.0.1                    	   //
#// Date: 01/07/19					  	      //
#// Desc: Does all eeg pre-processing and //
#//		  single-subject avging 	      //
#// eeg= eeg data; vis= vis events;       //
#// gc= good channels vector;			      //
#// Updated to include downsampling			//
#///////////////////////////////////////////

	library(signal)   # For filtering
	gc()

	eyechannels= c(31:32)
	goodchannels = which(c(1:30)%in%gc)
	gc()

	#vistrg <- cbind(as.numeric(vis[,3]),as.numeric(vis[,2]))
	adjustvis <- 501 - adjust #find the new "0", aka the start time of the visual stimulus
	vistrg <- cbind(adjustvis,as.numeric(vis[,2]))
	
#re-epoch the data with respect to the visual stimulus
eegepoched=array(0,dim=c(nrow(eeg),32,1201))
for(i in 1:nrow(vistrg)){
	eegepoched[i,,] <- eeg[i,,(vistrg[i,1]-200):(vistrg[i,1]+1000)]
}
	 
tempindex=c(1:nrow(vistrg))
for(i in 1:nrow(vistrg)) {
	timepoints=c(1:1201)
	test<- apply(eegepoched[i,gc,timepoints],1,function(x) any(abs(x)>100) ) 
	if(any(test==TRUE)) {
		tempindex= tempindex[! tempindex %in% i] 
		gc()
		}
	eegclean <- eegepoched[tempindex,,]
	gc()
 }
gc()

#Events for clean data, removes the bad trials from trial_labels
final_trg= vistrg[tempindex,]
final_trg_resp= behav[tempindex,]
gc()
	 
#CAR re-referencing
chansmean=apply(eegclean[,goodchannels,],c(1,3),mean)
eegCAR= sweep(eegclean,c(1,3),chansmean,FUN="-")
gc()


#Baseline Correction
BCmean=apply(eegCAR[,,1:200],c(1,2),mean) #takes the mean of the first 500 ms for all trials and all electrodes #note: 07/20/16 changed this to 200 ms
eegBC= sweep(eegCAR,c(1,2),BCmean,FUN="-") #performs the baseline correction procedure by subtracting the values in BCmean from each epoch per electrode
gc() 


#Find trials for C1 analysis of visual stimulus; any SOA >= 100 ms
C1trgs= which(final_trg[,2] %in% c(11,13,15,17,19))
gc()

C1trials=eegBC[C1trgs,,]


#Trial type averaging for C1 events
eegavgevents= array(0,dim=c(6,ncol(C1trials),length(C1trials[1,1,])))

for(i in c(11,13,15,17,19)){
	if(i == 11)
		eegavgevents[1,,]=apply(C1trials[which(final_trg[C1trgs,2]==11),,],c(2,3),mean)
	else if(i == 13)
		eegavgevents[2,,]=apply(C1trials[which(final_trg[C1trgs,2]==13),,],c(2,3),mean)
	else if(i==15)
		eegavgevents[3,,]=apply(C1trials[which(final_trg[C1trgs,2]==15),,],c(2,3),mean)
	else if(i==17)
		eegavgevents[4,,]=apply(C1trials[which(final_trg[C1trgs,2]==17),,],c(2,3),mean)
	else if(i==19)
		eegavgevents[5,,]=apply(C1trials[which(final_trg[C1trgs,2]==19),,],c(2,3),mean)
}
gc()
eegavgevents[6,,] <- apply(C1trials,c(2,3),mean) # take the average across all event types
rownames(eegavgevents) <- c("100 ms","120 ms","140 ms","160 ms","180 ms","avgd")

return(eegavgevents) 
}



