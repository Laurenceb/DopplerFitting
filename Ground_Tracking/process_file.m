function status=process_file(offset)
	####### config
	output_rate=25000;%output will be 25ksps I and Q, enough for +-28ppm at 868mhz
	name_append="_processed.wav";
	processing_chunk=1e6;%process 1M samples at a time
	if(nargin()==0)
		offset=0;%the offset input is the band center offset in Hz, usually the SDR should be tuned so thats its close to zero, but other values can be passed
	endif
	####### end config
	status=0;
	files=ls("*.wav");%find the first wav file in pwd
	if(index(files(1,:),'.'))
		file=files(1,:);%first element is blank for some reason?
	else
		file=files(2,:);
	endif
	file=strtrim(file);
	info=stat(file);
	size_=info.size-44;%size of data payload in bytes	
	size_/=2;	%size in samples*channels (at this point we don't know if its a 2 or 4 channel file, 4 channel has 2 channels of timestamping)
	[a,spr]=wavread(file,[1,1]);%spr is the incoming sample rate
	channels=length(a);%number of data channels
	if(channels!=2 && channels!=4)
		printf("Error, %d channels, should be 2 or 4\n",channels);
		status=1;
		return;
	endif
	printf("Processing file:'%s', with %d channels at %.1fsps\n",file,channels,spr);fflush(stdout);
	if(channels==4)%4channels - timestamps are present, so create a csv file
		timefile=[file(1:end-4),name_append(1:end-4),".csv"];%csv file to store timestamps and raw sample offsets
		starttime=timeconvert(a(3:4));
		timesstamps=[0,starttime];%the first column is relative sampling clock time, the second is system time (starting with seconds since offset)
		a=wavread(file,[2,2]);
		timesstamps=[timestamps;0,timeconvert(a(3:4))];
	endif
	samples=floor(size_/channels);%number of samples
	output=[file(1:end-4),name_append];%name of the output file
	phase_samp=2*pi*offset/spr;%phase per input sample
	endphase=1;		%phase at the end of each chunk (for carry over)
	phasors=exp(-i*[0:processing_chunk-1].*phase_samp);
	[a,b]=butter(4,output_rate*0.75/spr);%4th order Butterworth low pass covering 75% of the output specral range
	filterstate=zeros(4,1);%holds filter state between iterations
	sampleoffset=1;		%used to track decimator offset
	stridefactor=spr/output_rate;%the sample stride
	filtered=zeros(ceil(samples/stridefactor),1);%holds the filtered data
	m=1;
	for n=1:processing_chunk:samples
		this_samples=min([samples-n,processing_chunk]);%number of samples to be processed in this iteration
		sp=wavread(file,[n,n+this_samples-1]);
		sp_=sp(:,1).+i.*sp(:,2);%complex samples
		sp_.*=phasors(1:this_samples)'.*endphase;%mix up/down
		[sp_,filterstate]=filter(a,b,sp_,filterstate);%low pass filter the data
		endphase.*=phasors(this_samples);%the new beginning phase for the next block of incoming data
		strides=[sampleoffset:stridefactor:this_samples];
		sampleoffset=strides(end)+stridefactor-processing_chunk;%wrap around to next chunk
		filtered(m:m+length(strides)-1)=sp_(round(strides));%sample the filtered baseband
		m+=length(strides);
		if(channels==4)
			times=timeconvert(sp(:,3:4));%in units of nanoseconds
			indices_=find(diff(times)>0).+1;
			indices=indices_.+n;%indices into the unresampled sample set where new blocks of samples arrive
			timestamps=[timestamps;indices./spr,times(indices_)];%first column is sampling clock time of arrival (nanoseconds from 0), second system time
		endif
		printf("Processing... %.1f%%		\r",(n+(processing_chunk/2))/samples*100.0);fflush(stdout);%progress dialog
	endfor
	filtered=[real(filtered),imag(filtered)];%move back to real notation
	wavwrite(filtered,output_rate,output);%write the filtered file
	if(channels==4)
		csvwrite(timefile,timestamps);
	endif
	printf("\nDone\n");
	return;
endfunction

function time_=timeconvert(samples)
	samples.*=2^15;
	samples(find(samples<0))+=2^15;
	time_=samples(:,1)*2^16+samples(:,2);
endfunction
