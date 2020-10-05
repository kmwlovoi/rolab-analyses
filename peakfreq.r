#Ocular artifacts removed via ICA in eeglab
#Data is already re-referenced and segmented so that the data starts at the beginning of the first event

	library(R.matlab) # For reading in the .mat files
	library(wavelets) # For the FFT
	library(signal)   
	library(abind)    #For binding 3d matrices
	library(pracma)
	
	poly <- stats::poly
	
	
	importSet <- function(EEGsetname){
		set <- readMat(EEGsetname)
		
		fdt_file <- paste(strsplit(EEGsetname,c("."),fixed=TRUE)[[1]][1],".fdt",sep="")

		var_names <- dimnames(set$EEG)[[1]]

		n_chans <- set$EEG[[which(var_names == "nbchan")]]
		n_trials <- set$EEG[[which(var_names == "trials")]]
		times <- set$EEG[[which(var_names == "times")]]
		
		EEG <- readBin(fdt_file,
							 "double",
							 n = n_chans * n_trials * length(times),
							 size = 4,
							 endian = "little")
		
		gc()
		
		dim(EEG) <- c(n_chans, length(times),max(n_trials, 1))

		return(EEG)
	}
	
	detrendfft <- function(alltf_mean,freqsforfft){
			y= alltf_mean[2:151]
			x= freqsforfft[2:151]
			x2 <- x^2
			quadratic.model <-lm(y ~ x + x2)
			predicted <- predict(quadratic.model,list(x=freqsforfft[2:151], x2=(freqsforfft[2:151])^2))
			alltf_det = alltf_mean[2:151]-predicted
			
			return(alltf_det)	
	}
	
	returnpeaks <- function(alltf_det,freqsforfft){
			alph_range = which(freqsforfft>=7&freqsforfft<=13)
			y= alltf_det[alph_range[1]:alph_range[length(alph_range)]]
			x= freqsforfft[alph_range[1]:alph_range[length(alph_range)]]

			tmp = findpeaks(y[1:length(y)], zero="0",minpeakheight = -Inf, minpeakdistance=10,threshold = 0, npeaks = 0, sortstr = FALSE)
		
			return(tmp[which(x[tmp[,2]]>=1)[1],2]+(alph_range[1]-1)) #Return the first value (the largest peak) in the second column (the peak location) that is greater than or equal to 1 Hz; 70 instead of 71 because starts at an index of 1
	}
	
	
	#Set up the importing info
	subj = c(10:12) #EEG data wasn't recorded for S4
	label = 'eps_ICA'
	
for(i in 1:length(subj)){
   s = subj[i] #Which subject will be analyzed? 
	nchans = 6 #Number of channels for analysis; 4 if averaging across channels, 12 if using individual electrodes
	all_channels = c("Fp1","Fp2","F3","F4","C3","C4","P3","P4","O1","O2","F7","F8","T7","T8","P7","P8","Fz","Cz","Pz","Iz","FC1","FC2","CP1","CP2","FC5","FC6","CP5","CP6","TP9","TP10","31","32")		
	chosen_chans = c(7:10,19:20)
   notask_peakfreq = numeric(nchans) 
	#Pzfftovertime = array(0,dim=c(200,151))
	#Pzfreqovertime = numeric(200)
	alltf_mean = array(0,dim=c(nchans,151)) #For holding the mean log-transformed FFT
	alltf_det = array(0,dim=c(nchans,150)) #For holding the detrended FFT, excluding the 0 frequency component

	#Set up analysis parameters
	freq_min = 7 #minimum frequency of alpha band
	freq_max = 13 #maximum frequency of alpha band

	tracesfn = (paste("alpha/ICA_rem/s",s,label,".set",sep=""))
	traces <- importSet(tracesfn)
	traces <- aperm(traces,c(3,1,2))
	
	#alltf = array(0,dim=c(nrow(traces), ncol(traces[1,,]))) #For holding the FFT
	alltf_logtrans = array(0,dim=c(nrow(traces), 151)) #For holding the log-transformed FFT
	runchans = 1
	par(mfrow=c(4,3))		
	for(chan in chosen_chans){
		for(seg in 1:nrow(traces)){ #For each segment...
			freqsforfft = 500 * seq(0,1,(1/((ncol(traces[1,,])+1)/2)))  #the +1 is for the zero-frequency component
			alltf_logtrans[seg,] = 10*log10(((abs(fft(traces[seg,chan,])))^2)[1:151]) #Log-transform the FFT for the frequency analysis, just save data up to 50 Hz
		}			

		
		alltf_mean[which(chosen_chans==chan),] = colMeans(alltf_logtrans)
		alltf_det[which(chosen_chans==chan),] <- detrendfft(alltf_mean[which(chosen_chans==chan),],freqsforfft)
		#peak <- returnpeaks(alltf_det[chan-3,,f],freqsforfft) #OLD method-- we only want to use the detrended data when necessary though, so:
		peak <- returnpeaks(alltf_mean[which(chosen_chans==chan),],freqsforfft) #Remove the 0 freq component
		notask_peakfreq[which(chosen_chans==chan)] = freqsforfft[1:151][peak] #Extract the frequency, ignoring the zero-frequency component

	}



	assign(paste("s",s, "notask","peakfreq", sep='_'), notask_peakfreq)
	assign(paste("s",s, "notask","meanfft",sep='_'),alltf_mean)
	assign(paste("s",s, "notask","meanfftdet",sep='_'),alltf_det)


	rm(traces,chan,runchans,s,notask_peakfreq,nchans,alltf_logtrans,seg,peak,freq_max,freq_min,chosen_chans,all_channels,tracesfn)
}



