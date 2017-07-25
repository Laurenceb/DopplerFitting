function monte_attempt(stat1,stat2,balloon)
	times_=[0:260];
	stat1_=xform( 'pos' , stat1.*[pi/180,pi/180,1], -1, 0 );
	stat2_=xform( 'pos' , stat2.*[pi/180,pi/180,1], -1, 0 );
	for n=1:10
		v=generate_velocity();
		p=generate_position(balloon,v);
		xform( 'new', p, 0, 1 ); #new NED co-ordinate tranform at the starting point
		v=xform( 'vel' , v, 1, 0 );#velocity into ECEF
		posvel=[p,v];
		a=find_apogee(posvel,times_);
		printf("Run %d -------------\n",n);printf(" Actual apogee=%f km\n",a);fflush(stdout);
		[d,p_one,pos]=generate_doppler(posvel(1,1:3),posvel(1,4:6),stat1_,times_);
		[d,p_two,pos]=generate_doppler(posvel(1,1:3),posvel(1,4:6),stat2_,times_);
		observed_delta_range_=p_one-p_two;
		observed_delta_range_.-=mean(observed_delta_range_);
		[maxalt,position,velocity,difference]=rangefit(observed_delta_range_,times_,stat1,stat2,balloon);
		position
	endfor
endfunction


function v=generate_velocity()#generates a velocity in local NED space
	#this will run a Monte-Carlo with random trajectories
	v_r=randn(1,3).*[0.15,0.15,0.1];
	v_v=1300.*(1-v_r(3));#the speed of the vehicle
	v=[v_v*v_r(1:2),-v_v*sqrt(1-sum(v_r(1:2).^2))];#there is a random tip off error as well
endfunction

function p=generate_position(p_,v)#the input argument is the LLA and NED velocity, output is ECEF with a random perturbation and some vertical offset added
	p_(1:2).*=pi/180;
	north_factor=6630e3;
	east_factor=north_factor*cos(p_(1));
	error_=randn(1,3).*[100,100,250];#the position error distribution in meters in NED space (these are on top of the velocity drift)
	error_(1)/=north_factor;
	error_(2)/=east_factor;	#convert to position change in radians
	error_(3).+=2000;	#add 2 Km vertical offset
	error_(1:2).+=((v(1:2)./abs(v(3))).*2000)./[north_factor,east_factor];#now account for the velocity induced drift over the 2km ascent
	p=xform( 'pos' , p_.+error_, -1, 0 );	
endfunction

function a=find_apogee(posvel,times_)
	position=xform( 'pos' , posvel(1,1:3), 0, -1 );#the solved initial position (LLA)
	xform( 'new', posvel(1,1:3), 0, 1 ); #new NED co-ordinate tranform at the starting point
	velocity=xform( 'vel' , posvel(1,4:6), 0, 1 );#transform the velocity (NED)
	surface_projection=xform( 'pos' , [position(1:2),0], -1, 0 );#the geoid surface under the Rockoon in ECEF
	[d,p,pos]=generate_doppler(posvel(1,1:3),posvel(1,4:6),surface_projection,times_);# need altitude above surface at each point, pos is in ECEF
	pos=xform('pos',pos,0,-1);
	maxalt=max(pos(:,3));	#maxalt is the maximum altitude of the Rockoon
	a=maxalt*1e-3;
endfunction
