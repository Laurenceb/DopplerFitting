function [maxalt,position,velocity]=rangefit_sa(observed_delta_range_,times,station1,station2,balloon,station3)
	#input the observed range delta between the two stations (set to zero at times(1), and LLA of the stations and balloon)
	#these are as row vectors in WGS-84
	#the delta range is range_station_1-range_station2 with zero offset at the start of the data
	#it is possible to use either three or two stations as the input, for three stations the offsets are pairwise as rows, 
	threestations=0;
	#if(nargin>5 && size(observed_delta_range_,1)==2)#there is three station data (so two deltas)
	#	station3(1:2).*=pi/180;
	#	station3=xform( 'pos', station3, -1, 0 );
	#	observed_delta_range_=[observed_delta_range_(1,:),observed_delta_range_(2,:)];
	#	threestations=1;
	#endif
	station1(1:2).*=pi/180;	#convert to radians
	station2(1:2).*=pi/180;
	balloon(1:2).*=pi/180;
	xform( 'init', 1 );
	pos=xform( 'pos', balloon, -1, 0 );#move balloon into ECEF
	xform( 'new', pos, 0, 1 ); #new NED at the balloon
	vel=xform( 'vel', [0,0,-1400], 1, 0 );#assume 1.2km/s vertical velocity at the balloon, move it into ECEF
	station1=xform( 'pos', station1, -1, 0 );
	station2=xform( 'pos', station2, -1, 0 );
	posvel=[pos,vel];	#this is for the Rockoon initialisation
	#config ------------------------------------
	# SA controls
	# setting ub and lb to same value restricts that parameter, and the algorithm does not search
	ub = (posvel.+[1e4,1e4,1e4,300,300,300])';
	lb = (posvel.-[1e4,1e4,1e4,300,300,300])';
	nt = 20;
	ns = 5;
	rt = 0.25; # 0.5==careful - this is too low for many problems
	maxevals = 5000;
	neps = 5;
	functol = 1e-5;
	paramtol = 5e-3;
	verbosity = 2; # 2==full, 1==only final results. Inc
	minarg = 1;
	control = { lb, ub, nt, ns, rt, maxevals, neps, functol, paramtol, verbosity, 1}; 
	[posvel, obj_value, convergence] = samin("error_calc", {posvel', station1,station2,times,observed_delta_range_}, control);
		xform( 'new', pos, 0, 1 ); #new NED
		position=xform( 'pos' , posvel(1,1:3), 0, -1 );#the solved initial position (LLA)
		velocity=xform( 'vel' , posvel(1,4:6), 0, 1 );#and velocity (NED)
		surface_projection=xform( 'pos' , [position(1:2),0], -1, 0 );#the geoid surface under the Rockoon in ECEF
		[dop,pos]=generate_doppler(posvel(1,1:3),posvel(1,4:6),surface_projection,times);
		maxalt=max(pos);	#maxalt is the maximum altitude of the Rockoon
		printf("\nApogee at %f km\n",maxalt*1e-3);
		position(1:2).*=180/pi;		#this is now in degrees
endfunction


function f=error_calc(pos_vel,station1,station2,times,baseline_)
	pos_vel=pos_vel';
	[dop,baseline1]=generate_doppler(pos_vel(1,1:3),pos_vel(1,4:6),station1,times);
	[dop,baseline2]=generate_doppler(pos_vel(1,1:3),pos_vel(1,4:6),station2,times);
	baseline=baseline1-baseline2;
	baseline.-=mean(baseline);
	f=mean(abs(baseline.-baseline_));
	printf("\rMean error is: %f meters     ",f);fflush(stdout);
endfunction
