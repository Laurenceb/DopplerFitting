function [intphase,discrm,phasors_l,phasors_h,strengths,lock_indicator]=pll_tracker(filename)
#Config information here
#tone_space=830.6503;	#old setting
tone_space=557.8995;	#used on the latest recordings
#######################
[samples,spr]=wavread(filename);
if(size(samples)(2)==2)
	#samples(:,1).-=samples(:,2)*1.7;
	#samples(:,2).*=3.4;
	#samples_=randn(10001,1).+i*(+0.9-0.09*i)*randn(10001,1);#fundge
	#samples=[real(samples_),0.1*imag(samples_)];
	#simple correlation based calibration
	samples_=samples;	#store this for reference later
	samples=samples(:,1).+(i)*samples(:,2);#move to complex domain
	samples_fr=fft(samples);#move to frequency domain
	total_energy=samples_fr'*samples_fr;
	energy=samples_fr.*conj(samples_fr);
	[sorted,order]=sort(energy);#the bin energies are now in ascending order
	sumtotal=cumsum(sorted);
	tokeep=find(sumtotal>0.9*total_energy);#try to keep only the top 10% of the entire spectral energy
	printf("Keeping %f bins for correlation search\n",size(tokeep)(1));fflush(stdout);
	tokeep=order(tokeep);	#the FFT index bins that will be kept
	keepme=zeros(size(samples_)(1),1);
	keepme(tokeep)=1;
	halfsize=floor(size(samples_)(1)/2);#number of elements in half of the frequnecy space
	keepme(2:halfsize).+=fliplr(keepme(end+2-halfsize:end)')';#mirror keepme, so we keep if there is a spike at positive or negative
	keepme(end+2-halfsize:end)=fliplr(keepme(2:halfsize)')';
	keepme(find(keepme>0))=1;#fix range of keepme to 0,1
	samples_fr.*=keepme;	#dump the bins without spikes
	samples=ifft(samples_fr);samples=[real(samples),imag(samples)];#back to the time domain
	a=1/sqrt((samples(:,1)'*samples(:,1))/((samples(:,2)'*samples(:,2))-((samples(:,1)'*samples(:,2))/(samples(:,2)'*samples(:,2)))^2));#real part of true a
	if(imag(a))
		printf("Error: excessive mixing term - complex root found\n");fflush(stdout);
		a=1;
	else
		a+=-i*a*((samples(:,1)'*samples(:,2))/(samples(:,2)'*samples(:,2)));#calculate the imaginary part of the true a
	endif
	a=1/a;
	#try to do IQ calibration in frequency space - this need more work
	#r_=fft(real(samples_));
	#i_=fft(i*imag(samples_));
	#power=abs(r_).+abs(i_);
	#reject=find(power>quantile(power,0.5));
	#r_(reject)=0;i_(reject)=0;
	#numsamp=floor(size(r_)(1)/2);
	#r_=[r_(2:numsamp),conj(fliplr((r_(end+2-numsamp:end))')')];
	#i_=[i_(2:numsamp),conj(fliplr((i_(end+2-numsamp:end))')')];
	#arg(dot(r_(:,1),i_(:,2))-dot(r_(:,2),i_(:,1)))*180/pi
	#a_=dot(i_(:,1),i_(:,2));
	#b=dot(i_(:,1),r_(:,2))+dot(i_(:,2),r_(:,1));
	#c=dot(r_(:,1),r_(:,2));
	#a=(-b+sqrt((b^2)-(4*a_*c)))/(2*a_);
	printf("Phase:%f degrees, Amplitude:%f\n",arg(a)*180/pi,abs(a));fflush(stdout);
	#correct the samples
	samples=samples_;
	#plot(log(abs(fft(samples(:,1).+(i)*samples(:,2)))).*(20/log(10)),'r');hold on;
	samples=samples(:,1).+(i*a)*samples(:,2);
	#plot(log(abs(fft(samples))).*(20/log(10)));pause();	
endif
samples=fliplr(samples')';		#run in reversed time, this means the nastiest data will be towards the end
#samples.+=randn(size(samples))/7;	#test to add some noise to the data
samples(end-2e5:end)=randn(2e5+1,1)/7;	#mor noise test
numsamples=size(samples)(1);
times=[0:1/spr:(numsamples-1)/spr]; 	#time index for each sample
aquisition_chunk=samples(1:round(0.1*spr));#use a 0.1 second period at the end of the recording to try to find the carrier
#[m,indx]=max(abs(fft(aquisition_chunk)(1:floor(0.05*spr))));#find the highest frequency component - wont work due to wideband signal
#indx=1;
#loop through the carrier offsets looking for the best signal
#for m=(tone_space/2)+300:20:min([(spr/2)-(tone_space/2),77400])#run through the range of the search frequencies
#	[phasors_l,phasors_h,discrm,strength,intphase,lock_indicator]=track(m,tone_space,aquisition_chunk,spr);
#	strengths(indx)=strength;
#	carriers(indx)=m;
#	indx+=1;			#store the strength
#	printf("\rCarrier:%f Hz PLL energy:%f   ",m,strength);fflush(stdout);
#endfor
carriers=[0:1/0.1:0.1*spr];		#the carrier frequencies
strengths=abs(fft(aquisition_chunk));
strengths.*=circshift(strengths,round(0.1*tone_space));#multiply with the other bin - this technique works if tone spacing = bitrate*n, where n is an integer
#plot(strengths);pause();
[mx,maxindx]=max(strengths);
target=quantile(strengths,1-(4/size(strengths)(2)));#the fourth highest value
if(mx>target*3)				#need the peak to be more than three times the value of the 4th highest point
	printf("\nFound carrier at:%f Hz, PLL median energy is:%f\n",carriers(maxindx),mx);
	fflush(stdout);
else
	printf("\nDid not find a carrier, max signal was: %f Hz, PLL median energy:%f\n",carriers(maxindx),mx);
	break;
endif
#now we know where the carrier is, we can run over the file using the correct carrier initialisation
[phasors_l,phasors_h,discrm,strength,intphase,lock_indicator]=track(carriers(maxindx),tone_space,samples,spr);
phasors_l=fliplr(phasors_l);
phasors_h=fliplr(phasors_h);
discrm=fliplr(discrm);
intphase=fliplr(intphase(2:end));
intphase=-intphase.+intphase(1);	#intphase initialised at zero
lock_indicator=fliplr(lock_indicator);
lock_indicator=1.5.-lock_indicator;	#positive values no indicate a lock
endfunction

function [phasors_l,phasors_h,discrm,strength,integrated_phase,lock_indicator]=track(carrier,tone_space,samples,spr)
#filter spr
fratio=220;
fspr=spr/fratio;			#the pll loop filter runs at ~200hz
#the baseband filter
[baseband_b,baseband_a]=butter(4,tone_space*1.5/spr);#this could be narrowed up to exclude noise, theoretically can go as low as 1
#the discriminator filter
[discr_b,discr_a]=butter(3,100/spr);	#a 50hz low pass for the discriminator
#PLL P term
PLLP=35;
#PLL I term
PLLI=40;
LBW=10;
zeta=0.707;				#critical damping
k=0.25;					#from Kai Borre example
% Solve natural frequency
Wn = LBW*8*zeta / (4*zeta.^2 + 1);
% solve for t1 & t2
tau1 = k / (Wn * Wn);
tau2 = 2.0 * zeta / Wn;
#PLL P term
PLLP=tau2/tau1;
#PLL I term
PLLI=1/tau1;
#carrier that is passed is the central frequency
integrated_phase=0;
PLL_Integrator=carrier-tone_space/2;	#we define the lowest frequency as the carrier
phase_low=1;				#phase of the LO
phase_high=exp(-i*pi/2);		#second LO for the upper carrier
phase_low_chunk=phase_low.*exp(-2*i*pi*PLL_Integrator*[0:fratio-1]/spr);
phase_high_chunk=phase_high.*exp(-2*i*pi*(PLL_Integrator+tone_space)*[0:fratio-1]/spr);
for n=1:fratio:size(samples)(1)-fratio
	phasor_low_=samples(n:n+fratio-1)'.*phase_low_chunk;
	phasor_high_=samples(n:n+fratio-1)'.*phase_high_chunk;
	if n>1
		[phasor_low,low_filterstate]=filter(baseband_b,baseband_a,phasor_low_,low_filterstate);#filter the two baseband signals
		[phasor_high,high_filterstate]=filter(baseband_b,baseband_a,phasor_high_,high_filterstate);
	else
		[phasor_low,low_filterstate]=filter(baseband_b,baseband_a,phasor_low_);#filter the two baseband signals
		[phasor_high,high_filterstate]=filter(baseband_b,baseband_a,phasor_high_);
	endif
	discriminator_low=asin((imag(phasor_low).*sign(real(phasor_low)))./abs(phasor_low));
        discriminator_high=asin((imag(phasor_high).*sign(real(phasor_high)))./abs(phasor_high));#in radians
	discriminator_low(find(isnan(discriminator_low)))=0;
	discriminator_high(find(isnan(discriminator_high)))=0;
	#if n>1
	#	[disc_l,disc_filterstate_l]=filter(discr_b,discr_a,discriminator_low,disc_filterstate_l);#filter the discriminator
	#	[disc_h,disc_filterstate_h]=filter(discr_b,discr_a,discriminator_high,disc_filterstate_h);
	#else
	#	[disc_l,disc_filterstate_l]=filter(discr_b,discr_a,discriminator_low);#we need to run the iterative filter as appropriate
	#	[disc_h,disc_filterstate_h]=filter(discr_b,discr_a,discriminator_high);
	#endif
	disc_l=discriminator_low;
	disc_h=discriminator_high;
	PLL_Integrator+=PLLI*(mean(disc_l.+disc_h))/(2*fspr);#run the PI filter, I is shared
	carrier_l=PLLP*mean(disc_l)+PLL_Integrator;#the carrier is adjusted as appropriate
	phase_low=phase_low_chunk(end);	#the phasors for the upper and lower carriers are generated to match with previous phase
	phase_low/=abs(phase_low);	#normalize the phasors
	carrier_h=PLLP*mean(disc_h)+PLL_Integrator+tone_space;
	phase_high=phase_high_chunk(end);
	phase_high/=abs(phase_high);	#normalize the phasors
	#use the phasors to generate chunks
	phase_low_chunk=phase_low.*exp((-2*i*pi*carrier_l).*[1:fratio]./spr);
	phase_high_chunk=phase_high.*exp((-2*i*pi*carrier_h).*[1:fratio]./spr);
	#the outputs are inband phasor, relative signal strength, and integrated phase
	phasors_l(n:n+fratio-1)=phasor_low;#the complex baseband
	phasors_h(n:n+fratio-1)=phasor_high;
	discrm(:,n:n+fratio-1)=[disc_l;disc_h];#the (filtered) discriminator output
	rawdis(:,n:n+fratio-1)=[discriminator_low;discriminator_high];
	integrated_phase(n:n+fratio-1)=((carrier_l+tone_space/2).*[1:fratio]./spr).+integrated_phase(end);#this is at the full sample rate, in units of cycles
endfor
#the signal "strength"; a scalar indicating the relative amount of inphase energy in the output channel
strength=(median([abs(real(phasors_l)),abs(real(phasors_h))])/median([abs(imag(phasors_l)),abs(imag(phasors_h))])-1)...
*median([abs(phasors_l),abs(phasors_h)])/std(samples);
[b,a]=butter(3,30/spr);			#a 15hz low pass filter for this
lock_indicator=filter(b,a,sum(abs(rawdis-discrm)));
endfunction


