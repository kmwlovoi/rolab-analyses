loadtrg = function(filename) {
    trgfile = paste(filename,".trg",sep="")
    trgfileSize = file.info(trgfile)$size
    rawdat = readBin(trgfile, what = "raw", n=trgfileSize) # read file into raw
    nsamples = trgfileSize/(4+2) # 4 bytes for timestamp and 2 bytes for trigger
    dim(rawdat) =  c(4+2,nsamples)
    timestamps = readBin(con=rawdat[1:4,], integer(), size=4, n=nsamples)
    triggers = readBin(con=rawdat[5:6,], integer(), size=2, n=nsamples)
    trgdata_preshaped = cbind(floor((timestamps-timestamps[1]+(timestamps[1]%%160000))/16),triggers)
    trgdata = trgdata_preshaped[2:nsamples,]
	return(trgdata)
}    

loadeeg = function(filename) {
    eegfile = paste(filename,".eeg",sep="")
    eegfileSize = file.info(eegfile)$size
    nchannels = 16 # number of channels in eeg recording
    bytesperchan = 2 # number of bytes in eeg measurement per channel
    nsamples = eegfileSize/bytesperchan # calculate number of samples in the data
    channeldata_preshaped = readBin(eegfile, integer(), size=2, n=nsamples)
    channeldata = (matrix(channeldata_preshaped, nchannels,(length(channeldata_preshaped))/nchannels)) # rows are channels, columns are timepoints
    return(channeldata)
}

filter = function(eegdata,freqpass,hiorlo) {
    bf <- butter(4, freqpass/500, type=hiorlo)
    filtdata <- filtfilt(bf, eegdata[,])
    return(filtdata)
}

epoch = function(trgdata,eegdata,channels,epochstart,epochend) { # epochstart and epochend times are relative to 
    all_trials = array(0,dim=c(nrow(trgdata),channels,abs(epochstart)+abs(epochend)+1))
    for (s in seq_along(1:nrow(trgdata))) { # for all trials
       trgtime = trgdata[s]
       all_trials[s,,] = eegdata[1:channels,(trgtime+epochstart):(trgtime+epochend)]
    }
	return(all_trials)
}

erps = function(trgdata,all_trials,channel,conds) {
    erp = array(0,dim=c(length(conds),length(all_trials[1,1,]))) 
    par(bg="black",fg="white",col.main="white",col.lab="white",col.axis="white",font.axis=2)
    plot((colMeans(all_trials[which(trgdata[,2] == 1),channel,])-2048)/8.192,col=cm[1],type="l",xaxt="n",ylim=range(-200,200),yaxp=c(-200,200,8),font.main=2,cex.lab=1.2,xlab='Time (in ms)',ylab="Amplitude (in microVolts)")
    axis(1, at=seq(from=0,to=400,by=100),labels=seq(from=-100,to=300,by=100),col="white")
    cm = colors()
    for (s in seq_along(conds)) {
       lines((colMeans(all_trials[which(trgdata[,2] == s),channel,])-2048)/8.192,col=cm[s*9],lty=2)
       erp[s,] = colMeans(all_trials[which(trgdata[,2] == s),channel,])
    }
    return(erp)
}

scaleeeg = function(eegdata) {
	eegscaled= (eegdata-2048)/8.192
	return(eegscaled)
}


