[maxalt,position,velocity,difference]=process_central()
	addpath("../Inversion");
	s=load("-text","../config.txt",'station1','station2','station3','balloon','carrier');#the config data is stored in Octave text format (can be hand edited)
	#[intphase,discrm,phasors,lock_indicator,sampletime,strengths]=pll_tracker_s(filename,trim)<-this is done by the ground stations
	t=load("-binary","../station1.dat",'intphase','sampletime','phasors');#first ground station data
	u=load("-binary","../station2.dat",'intphase','sampletime','phasors');#second ground station data
	starttime=[t.sampletime(1),u.sampletime(1)];
	if(exist('station3'))
		v=load("-binary","../station3.dat",'intphase','sampletime','phasors');#optional third station
		starttime=[starttime,v.sampletime(1)];
	endif
	#now align the data to the largest UTC incriment. Abort if the difference is too high (>150ms)
	if(range(starttime)>0.15)
		printf("Error, range of start times is %f \n",range(starttime));
		return;
	else
		start=max(starttime)
		timedata=[t.sampletime(1:10)',u.sampletime(1:10)'];
		if(exist('station3'))
			timedata=[timedata,v.sampletime(1:10)'];
		endif
		indices=[min(find(timedata(:,1)==start)),min(find(timedata(:,2)==start))]
		if(exist('station3'))
			indices=[indices,min(find(timedata(:,1)==start))];
		endif
		sizes=[length(t.sampletime),length(u.sampletime)];
		if(exist('station3'))
			sizes=[sizes,length(v.sampletime)];
		endif
		windowsamples=min(sizes-indices);
	endif
	observed_delta_range=u.intphase(indices(2):indices(2)+windowsamples)-u.intphase(indices(1):indices(1)+windowsamples)
	if(exist('station3'))
		observed_delta_range=[observed_delta_range;v.intphase(indices(3):indices(3)+windowsamples)-u.intphase(indices(1):indices(1)+windowsamples)];
	endif
	observed_delta_range.*=2.99792e8/(2*pi*s.carrier);#convert to actual distances
	times=[1:windowsamples].*diff(t.sampletime(4:5));
	if(exist('station3'))
		[maxalt,position,velocity,difference]=rangefit(observed_delta_range_,times,s.station1,s.station2,s.balloon,s.station3);
	else
		[maxalt,position,velocity,difference]=rangefit(observed_delta_range_,times,s.station1,s.station2,s.balloon);
	end
endfunction
