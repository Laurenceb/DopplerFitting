function [intphase,discrm,phasors,lock_indicator,sampletime,strengths]=pll_tracker_s(filename,trim)
#######################
[samples,spr]=wavread(filename);
if(size(samples)(2)==4)		#there is timestamping data in the file
	timestamps=1;
else
	timestamps=0;
endif
#get the tuning frequency from the filename
name_=strsplit(name,'_');
name_=name_{end};
name_=name_(1:end-10);
freq=1000*sscanf(name_,"%d",1);#Frequency in Hz
#get the timestamps
timestamps=convert_wav_int32(samples(:,3:4));
starttime=timestamps(1);	#this is unix time for the starting time
timestamps(1)=timestamps(2);	#other times are integer nanoseconds 
timestamps./=1000;		#time is now in microseconds (a more reasonable unit)
negstaps=find(diff(timestamps)<0).+1;#fix the 1 second wraparound, look for pointswhere the time change is negative. Offset by 1 due to 
steps=zeros(size(timestamps));
steps(negsteps)=1;
steps=cumsum(steps);
timestamps.+=steps*1e6;		#now in second units
arrival_indices=find(diff(timestamps));#sample indices where new samples arrived
arrival_times=timestamps(arrival_indices);
[X,realtime]=detrend(arrival_indices,arrival_times,1);#time as a function of sample index
used1=find(X<quantile(X,0.2));	#indices that are below the 20% quantile
[X,realtime]=detrend(arrival_indices(used1),arrival_times(used1),1);
quantiles=[0.7:-0.1:0.1];	#loop through potential quantile values for the second iteration and pick the optimal one 
quantiles=quantile(X,quantiles);
for n=1:length(quantiles)
	use_=find(X<quantile(n));
	errors(n)=std(X(use_))/sqrt(length(use_));#the estimated error for each quantile setting
endfor
[minima,index]=min(errors);
used2=find(X<quantile(X,quantiles(index)));#second iteration, pick the optimal quantile of remaining data in order to minimise the error on the real time fix
used2=used1(used2);
[X,realtime]=detrend(arrival_indices(used2),arrival_times(used2),2);#final timebase uses 2nd order polynomial
numsamples=size(samples,1);
realtime=[realtime,realtime(end)+diff(realtime)(end)*(numsamples-arrival_indices(used2(end))/diff(arrival_indices(used2))(end))];#extropolate a new sample at end
realtime=[realtime(1)-diff(realtime)(1)*(arrival_indices(used2(1))/diff(arrival_indices(used2))(1)),realtime];#new sample at start
arrival_indices=[1,arrival_indices(used2),numsamples];#for interpolation purposes
sampletime=interp1(arrival_indices,realtime,[1:numsamples]);#the actual time for each sample in microseconds
sampletime=samplestime*1e-6+starttime;%this is now in unix time units
#convert to complex signals
signal=signal(:,1).+i.*signal(:,2);
%test code
std(samples)
samples.+=randn(size(samples))/1;#0.2;#/0.12;
std(samples)
samples=samples(1:3e6);
samples(1:1e5)=randn(1e5,1);
samples(end-1e6:end)=randn(1e6+1,1)/1;
%end of test code
samples=fliplr(samples')';	#the processing is time reversed to allow accurate launch time to be determined
duration=length(samples)/spr;	#in seconds (note that this doesnt use true time units)
#[b, a] = butter (4, [0.226,0.231]);
#samples=filter(b,a,samples);
increment=0;				
while 1
	aquisition_chunk=samples((1:round(0.2*spr)).+increment);#use a 0.2 second period at the end of the recording to try to find the carrier
	increment+=spr;#loop using 1 second indexing
	if(increment/spr+0.2>=duration)
		printf("No signal found\r\n");fflush(stdout);
		return;
	endif
	#[m,indx]=max(abs(fft(aquisition_chunk)(1:floor(0.05*spr))));#find the highest frequency component - wont work due to wideband signal?
	numsmp=size(aquisition_chunk,1);
	window=blackmanharris(numsmp)';		#typical airspy tcxo + tx xtal + doppler range at 868Mhz carrier
	search_indices=[1:floor(0.2*10000/spr)];#section of the FFT output to use
	carriers=(search_indices.-1).*5;	#in Hz
	numshift=length(search_indices);
	search_indices=[search_indices,(numsmp.-search_indices.+1)];#also use the negative frequencies
	carriers=[carriers,-(search_indices.*5)];#negative frequencies in Hz
	strengths=abs(fft(aquisition_chunk.*window)(search_indices));#only look in the +-10kHz range
	strengths=circshift(strengths,-numshift);		#to avoid possible wrap around errors around zero
	carriers=circshift(carriers,-numshift);
	#carriers=[500:5:0.35*spr].-5;		
	#indx=1;
	#loop through the carrier offsets looking ofr the best signal
	#for m=2000:10:min([(spr/2),8000])#run through the range of the search frequencies
	#	[phasors,discrm,strength,intphase,lock_indicator]=track(m,aquisition_chunk,spr);
	#	strengths(indx)=strength;
	#	carriers(indx)=m;
	#	indx+=1;			#store the strength
	#	printf("\rCarrier:%f Hz PLL energy:%f   ",m,strength);fflush(stdout);
	#endfor
	#plot(strengths);#pause();
	mx=diff(-sign(diff(strengths)));
	strengths(find(mx<=0).+1)=0;		#blank out everything accept the local maxima
	[mx,maxindx]=max(strengths);
	strengths(isnan(strengths))=0;
	target=quantile(strengths',1-(4/size(strengths)(1)));#the fourth highest value
	if(isnan(target))
		target=0;
	endif
	if(mx>target*2.5)			#need the peak to be more than 2.5 times the value of the 4th highest point
		printf("\nFound carrier at:%.1f Hz, PLL median energy is:%.0f\n",carriers(maxindx)+freq,mx);
		fflush(stdout);
		plot(strengths)
		break;
	else
		printf("\rDid not find a carrier, max: %.1f Hz, energy:%.0f, index %.0f seconds	",carriers(maxindx),mx,increment/spr);
		fflush(stdout);
	endif
endwhile
#just to the correct point in the data to begin tracking with the PLL
samples_=samples(increment-spr:end);
#now we know where the carrier is, we can run over the file using the correct carrier initialisation. Note that this is time reversed
[phasors,discrm,strength,intphase,lock_indicator]=track(carriers(maxindx),samples_,spr);
intphase=-intphase.+intphase(1);	#intphase initialised at zero
locked=zeros(size(lock_indicator));
locked(find(lock_indicator<0))=-1;	#points where pll is not locked
locked(find(lock_indicator>0.35))=1;	#points where pll is well locked 
[a,b]=butter(3,2/spr);			#1hz lowpass
locked=filtfilt(a,b,locked);
locked(find(locked>0.4))=1;
locked(find(locked<=0.4))=0;
lock_indicator=locked;			#1==lock, 0==unlocked
#the last 1 value indicates the point where lock was lost
locklost=max(find(lock_indicator));
numsamples_=length(samples_);		#we are using the trimmed samples
timelost=sampletime(numsamples-locklost);#time at which lock was lost. Subtract due to the reversed order
#estimate the initial carrier frequency by measuring change in integrated phase over 0.2s 
deltaphase=intphase(locklost)-intphase(round(locklost-0.2*spr));#integrated phase is in cycle units
deltatime=sampletime(locklost)-sampletime(round(locklost-0.2*spr));#calculate the time difference in real time
launchcarrier=deltaphase/deltatime;
#trim the data if requested
if(trim)
	indices=[1:locklost+round(0.2*spr)];	#these are the indices that will be kept
	phasors=phasors(indices);	#trim the data
	discrm=discrm(indices);
	intphase=intphase(indices);
	lock_indicator=lock_indicator(indices);
	sampletime=sampletime(indices);
endif
phasors=fliplr(phasors);
discrm=fliplr(discrm);
intphase=-fliplr(intphase);
lock_indicator=fliplr(lock_indicator);
#print out some debug info 
printf("Carrier at: %.1f Hz\n",launchcarrier+freq);
printf("Launch at: %s\n",ctime(timelost));fflush(stdout);
#now resample the data to 20hz, aligned to 50ms incriments UTC
time_strt=round(sampletime(1)*20)./20;
times=[time_strt:1/20:sampletime(end)];#50ms intervals
phasors=interp1(sampletime,phasors,times);
discrm=interp1(sampletime,discrm,times);
intphase=interp1(sampletime,intphase,times);
lock_indicator=interp1(sampletime,lock_indicator,times);#resample all of these
endfunction

function [phasors,discrm,strength,integrated_phase,lock_indicator]=track(carrier,samples,spr)
#filter spr
fspr=200;				#desired loop filter rate
fratio=round(spr/fspr);			#the pll loop filter runs at ~200hz
fspr=spr/fratio;			#the actual speed that the pll loop filter will be running at 
#the discriminator filter
[discr_b,discr_a]=butter(3,100/spr);	#a 50hz low pass for the discriminator
LBW_=30;
LBW=LBW_*2*pi;				#in Hz (convert to radians/s?)
#LBW=LBW_*2*pi/fspr;			#in units of radians/iteration
zeta=0.707;				#critical damping
k=2*pi;					#0.25 from Kai Borre example, here dco is 2*pi and discrim is 1->2/pi, giving k_dco*k_disc=2*pi->4
T=1/fspr;
#the baseband filter
LBW_=9;					#this is lower than the Kai Borre fudge factor version
[baseband_b,baseband_a]=butter(1,1.7*LBW_/(0.5*spr));#this could be narrowed up to exclude noise, theoretically can go as low as 1 or therabouts
% Solve natural frequency (lowered slightly due to the damping)
Wn = LBW*8*zeta / (4*zeta.^2 + 1);
% solve for t1 & t2
#tau1 = k / (Wn * Wn);
#tau2 = 2.0 * zeta / Wn;
#PLL P term
#PLLP=1+(tau2/tau1)/fspr;
#PLL I term
#PLLI=1/tau1;
#from http://read.pudn.com/downloads154/sourcecode/app/685304/GNSS_SDR/A%20Software-Defined%20GPS%20and%20Galileo%20Receiver.pdf p104
PLLP=1/k*8*zeta*Wn*T/(4+4*Wn*zeta*T+(Wn*T)^2)
PLLI=PLLP*Wn*T/(2*zeta)
PLLI=44;
PLLP=10;
#window=[[1/15:1/15:1],[1:-1/15:1/15]];
#carrier that is passed is the central frequency
integrated_phase=0;
PLL_Integrator=carrier;
phase=1;				#phase of the LO
phase_chunk=phase.*exp(-2*i*pi*PLL_Integrator*[0:fratio-1]/spr);
#preallocate
phasors=zeros(1,size(samples)(1));
discrm=phasors;
rawdis=phasors;
integrated_phase=phasors;
lock_indicator=phasors;
lastphase=0;
for n=1:fratio:size(samples)(1)-fratio
	phasor_=samples(n:n+fratio-1)'.*phase_chunk;
	if n>1
		[phasor_,filterstate]=filter(baseband_b,baseband_a,phasor_,filterstate);#filter the two baseband signals
	else
		[phasor_,filterstate]=filter(baseband_b,baseband_a,phasor_);#filter the two baseband signals
	endif
	#discriminator_=sum(imag(phasor_))*sign(sum(real(phasor_)))/sum(abs(phasor_));
	discriminator_=arg((phasor_(end)));
	discriminator_(find(isnan(discriminator_)))=0;
	#if n>1
	#	[disc_,disc_filterstate_]=filter(discr_b,discr_a,discriminator_,disc_filterstate_);#filter the discriminator
	#else
	#	[disc_,disc_filterstate_]=filter(discr_b,discr_a,discriminator_);#we need to run the iterative filter as appropriate
	#endif
	disc_=discriminator_;
	PLL_Integrator+=PLLI*T*disc_;#/(fspr);#run the PI filter, I is shared
	carrier=PLLP*disc_+PLL_Integrator;#the carrier is adjusted as appropriate
	phase=phase_chunk(end);	#the phasors for the upper and lower carriers are generated to match with previous phase
	phase/=abs(phase);	#normalize the phasors
	#use the phasors to generate chunks
	phase_chunk=phase.*exp((-2*i*pi*carrier).*[1:fratio]./spr);
	#the outputs are inband phasor, relative signal strength, and integrated phase
	phasors(n:n+fratio-1)=phasor_;#the complex baseband
	discrm(:,n:n+fratio-1)=disc_;#the (filtered) discriminator output
	rawdis(:,n:n+fratio-1)=discriminator_;
	integrated_phase(n:n+fratio-1)=((carrier).*[1:fratio]./spr).+lastphase;#this is at the full sample rate, in units of cycles
	lastphase=integrated_phase(n+fratio-1);#integrate this
	#(mean(real(phasor_)).^2-mean(imag(phasor_)).^2)/mean(abs(phasor_)).^2
	lock_indicator(n:n+fratio-1)=repmat(1-2*abs(mean(discriminator_))/pi,1,fratio);#greater than zero indicates a lock
endfor
#the signal "strength"; a scalar indicating the relative amount of inphase energy in the output channel
strength=(median(abs(real(phasors)))/median(abs(imag(phasors)))-1)...
*median(abs(phasors))/std(samples);
[b,a]=butter(3,30/spr);			#a 15hz low pass filter for this
#lock_indicator=filter(b,a,sum(abs(rawdis-discrm)));
lock_indicator=filter(b,a,lock_indicator);
endfunction

%this function recovers a 32bit integer stored across two channels of a wav file
function integer=convert_wav_int32(ints)
	integer=ints(:,1)*2^15;
	integer(find(integer<0)).+=2^16;
	integer+=ints(:,2)*2^(15+16);
	integer(find(ints(:,2)<0)).+=2^32;
endfunction

