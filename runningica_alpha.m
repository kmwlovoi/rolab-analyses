EEGs = { 's1_epsrs','s2_epsrs','s3_epsrs','s5_epsrs','s6_epsrs','s7_epsrs','s8_epsrs','s9_epsrs'};
plusica = '_ica';
extn = '.mat';
path = 'C:\\Users\\kwebster\\odrive\\OneDrive For Business\\PhosTOJStudy\\alpha\\preproc\\';
pathnew = 'C:\\Users\\kwebster\\odrive\\OneDrive For Business\\PhosTOJStudy\\alpha\\ICA_notrem\\';

for i=1:length(EEGs)
	fn = strcat(path,EEGs(i),extn)
	EEG = pop_importdata('dataformat','matlab','nbchan',32,'data',fn,'setname',EEGs(i),'srate',1000,'pnts',3001,'xmin',-3,'chanlocs','{ ALLEEG(1).chanlocs ALLEEG(1).chaninfo ALLEEG(1).urchanlocs }');
	EEG = eeg_checkset( EEG );
	EEG = eeg_checkset( EEG );
	EEG = pop_runica(EEG, 'extended',1,'interupt','on');
	EEG = eeg_checkset( EEG );
	fnnew = strcat(EEGs(i),plusica,extn)
	EEG = pop_saveset( EEG, 'filename',char(fnnew),'filepath',pathnew);
	ALLEEG = pop_delset( ALLEEG, [2] );
end
