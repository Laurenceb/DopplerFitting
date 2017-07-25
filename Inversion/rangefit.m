function [maxalt,position,velocity,difference]=rangefit(observed_delta_range_,times,station1,station2,balloon,station3)
	#input the observed range delta between the two stations (set to zero at times(1), and LLA of the stations and balloon)
	#these are as row vectors in WGS-84
	#the delta range is range_station_1-range_station2 with zero offset at the start of the data
	#it is possible to use either three or two stations as the input, for three stations the offsets are pairwise as rows, 
	threestations=0;
	if(nargin>5 && size(observed_delta_range_,1)==2)#there is three station data (so two deltas)
		station3(1:2).*=pi/180;
		station3=xform( 'pos', station3, -1, 0 );
		observed_delta_range_=[observed_delta_range_(1,:),observed_delta_range_(2,:)];
		threestations=1;
	endif
	station1(1:2).*=pi/180;	#convert to radians
	station2(1:2).*=pi/180;
	balloon(1:2).*=pi/180;
	xform( 'init', 1 );
	pos=xform( 'pos', balloon, -1, 0 );#move balloon into ECEF
	balloon=pos;		#store for reference
	xform( 'new', pos, 0, 1 ); #new NED at the balloon
	vel=xform( 'vel', [0,0,-1300], 1, 0 );#assume 1.2km/s vertical velocity at the balloon, move it into ECEF
	station1=xform( 'pos', station1, -1, 0 );
	station2=xform( 'pos', station2, -1, 0 );
	posvel=[pos,vel];	#this is for the Rockoon initialisation
	#approx entry point finder - ring velocity
	for n=1:36
		vect=[350*sin(n*pi/(18)),350*cos(n*pi/(18)),0];
		vect=xform( 'vel', vect, 1, 0 );
		vel=posvel(1,4:6).+vect;
		[dop,range1]=generate_doppler(posvel(1,1:3),vel,station1,times);
		[dop,range2]=generate_doppler(posvel(1,1:3),vel,station2,times);
		r=(range1-range2);		
		r.-=mean(r);
		delta=observed_delta_range_.-r;
		difference(n)=max(abs(delta));
		printf("%f%%\r",n/(36)*100);fflush(stdout);
	end
#	[a,b]=min(difference);
	b=sign(-diff(difference).*circshift(diff(difference)',-1)')(1:end-1);#changes in direction are positive
	g=diff(diff(difference));
	g(find(g<0))=0;
	b(find(b<0))=0;
	b.*=g;
	[a,c]=max(b);
	b(c)=0;	#wipe out the peak
	c+=1;	#the largest peak
	[a,c_]=max(b);#second highest peak
	c(2)=c_+1;
	vect_=[350*sin(c(1)*pi/(18)),350*cos(c(1)*pi/(18)),0];
	vect_=xform( 'vel', vect_, 1, 0 );
	vect=[350*sin(c(2)*pi/(18)),350*cos(c(2)*pi/(18)),0];
	vect=xform( 'vel', vect, 1, 0 );
	%From Monte-Karlo runs, c(2) is the best choice in around 65% of cases
	startposvel=[posvel;posvel];
	startposvel(1,4:6).+=vect_;
	startposvel(2,4:6).+=vect;
	################
	pos_=pos;
	for l=2:-1:1
		pos=pos_;
		posvel=startposvel(l,:);
		gain=0.15;
		err=1e5;
		if(threestations)
			gain=0.25;
		endif
		gain_initial=gain;
		convergence_factor=0.0;
		lasterror=0;
		convergepoint=0;
		for n=1:200
			[error_jacobian,baseline]=jacobian_calc(posvel,station1,station2,times);
			#cond(error_jacobian)
			if(threestations==1)
				[error_jacobian_,baseline_]=jacobian_calc(posvel,station1,station3,times);
				error_jacobian=[error_jacobian;error_jacobian_];
				baseline=[baseline,baseline_];
				posvel+=(pinv(error_jacobian)*(observed_delta_range_-baseline)')'.*gain;
			else
				if(err>1600)
					posvel(4:6)+=(pinv(error_jacobian(:,4:6))*(observed_delta_range_-baseline)')'.*gain;
				else
					posvel+=(pinv(error_jacobian)*(observed_delta_range_-baseline)')'.*gain;
				endif
			endif
			err=max(abs(observed_delta_range_-baseline));
			if(lasterror!=0)
				if(n!=2)
					if((err/lasterror)<4)
						convergence_factor+=0.2*((1-err/lasterror)-convergence_factor);#fact approaches 1 for perfect approach
					endif
				else
					convergence_factor=(1-err/lasterror);
				endif
			endif
			lasterror=err;
			printf("\rMax error is: %f meters      ,convergence factor:  %f		",err,convergence_factor);
			#posvel
			fflush(stdout);
			#if(err<30000) #-just some test code here
			#	gain=0.025;
			#end
			if(err<500)
				gain=0.75;
			endif
			if(err<6 && n>80)
				break;
			endif
			if(err<5 && n>50)
				break;
			endif
			if(err<2.5 && n>40)
				break;
			endif
			if(err<1 && n<25)
				break;
			endif
			if(err>1600 && threestations!=1)
				posvel(1:3)=pos;		#when error is high, adjust only velocity
			endif
			speed=sqrt(sum(posvel(4:6).^2));
			posoffset=posvel(1:3).-pos;
			if(speed>2000 || speed<500 )#||abs(posoffset(1))>10000||abs(posoffset(2))>10000 || abs(posoffset(3))<10000)#impossible situation, restart+ randomness
				posvel=[pos,vel]+randn(size(posvel)).*[500,500,500,50,50,50];
			endif
			#angle=sqrt(sum(posvel(4:5).^2))/abs(posvel(6));
			#if(angle>0.35)
			#	angle=0.35/angle;
			#	posvel(4:5).*=angle;
			#endif
			if(convergence_factor<(gain_initial/2) && n>8 && convergence_factor>0)
				convergence_factor=1;
				if(convergepoint && abs(convergepoint-err)/err<0.01)	#we converged twice to very similar error, likely to be the best poss solution
					break;
				end
				convergepoint=err;
				posvel=((posvel-[pos,vel]).*(randn(size(posvel))*0.05.+1)).+[pos,vel];#random noise
				posvel(1:3)=pos;	#not clear if restricting the position help at all?
				gain=gain_initial;
			endif
		endfor
		if n<200
			%expt dot product
			tstvec=posvel(1:3).-balloon;
			tstvec./=norm(tstvec);
			tstvec_=posvel(4:6)./norm(posvel(4:6));
			printf("\r\nDirection coefficient is:%f\r\n",dot(tstvec,tstvec_));fflush(stdout);
			%
			position=xform( 'pos' , posvel(1,1:3), 0, -1 );#the solved initial position (LLA)
			xform( 'new', posvel(1,1:3), 0, 1 ); #new NED
			velocity=xform( 'vel' , posvel(1,4:6), 0, 1 );#and velocity (NED)
			surface_projection=xform( 'pos' , [position(1:2),0], -1, 0 );#the geoid surface under the Rockoon in ECEF
			[d,p,pos]=generate_doppler(posvel(1,1:3),posvel(1,4:6),surface_projection,times);# need altitude above surface at each point, pos is in ECEF
			pos=xform('pos',pos,0,-1);
			maxalt=max(pos(:,3));	#maxalt is the maximum altitude of the Rockoon
			printf("\nApogee at %f km\n",maxalt*1e-3);
			position(1:2).*=180/pi;		#this is now in degrees
			if(dot(tstvec,tstvec_)>0.985)
				return;
			else
				printf("\r\nPoor solution, retrying with second entry point\r\n");fflush(stdout);
			endif
		else
			printf("\nDid not converge\n");
		endif
	endfor
endfunction


#The jacobian is calculated by brute force
function [error_jacobian,baseline]=jacobian_calc(pos_vel,station1,station2,times)
	[dop,baseline1]=generate_doppler(pos_vel(1,1:3),pos_vel(1,4:6),station1,times);
	[dop,baseline2]=generate_doppler(pos_vel(1,1:3),pos_vel(1,4:6),station2,times);
	baseline=baseline1-baseline2;
	baseline.-=mean(baseline);				#We only know the relative offset
	[dop,range1]=generate_doppler(pos_vel(1,1:3)+[5,0,0],pos_vel(1,4:6),station1,times);
	[dop,range2]=generate_doppler(pos_vel(1,1:3)+[5,0,0],pos_vel(1,4:6),station2,times);
	r=(range1-range2);
	b(:,1)=((r.-mean(r))-baseline)./5;
	[dop,range1]=generate_doppler(pos_vel(1,1:3)+[0,5,0],pos_vel(1,4:6),station1,times);
	[dop,range2]=generate_doppler(pos_vel(1,1:3)+[0,5,0],pos_vel(1,4:6),station2,times);
	r=(range1-range2);
	b(:,2)=((r.-mean(r))-baseline)./5;
	[dop,range1]=generate_doppler(pos_vel(1,1:3)+[0,0,5],pos_vel(1,4:6),station1,times);
	[dop,range2]=generate_doppler(pos_vel(1,1:3)+[0,0,5],pos_vel(1,4:6),station2,times);
	r=(range1-range2);
	b(:,3)=((r.-mean(r))-baseline)./5;
	[dop,range1]=generate_doppler(pos_vel(1,1:3),pos_vel(1,4:6)+[0.1,0,0],station1,times);
	[dop,range2]=generate_doppler(pos_vel(1,1:3),pos_vel(1,4:6)+[0.1,0,0],station2,times);
	r=(range1-range2);
	b(:,4)=((r.-mean(r))-baseline).*10;
	[dop,range1]=generate_doppler(pos_vel(1,1:3),pos_vel(1,4:6)+[0,0.1,0],station1,times);
	[dop,range2]=generate_doppler(pos_vel(1,1:3),pos_vel(1,4:6)+[0,0.1,0],station2,times);
	r=(range1-range2);
	b(:,5)=((r.-mean(r))-baseline).*10;
	[dop,range1]=generate_doppler(pos_vel(1,1:3),pos_vel(1,4:6)+[0,0,0.1],station1,times);
	[dop,range2]=generate_doppler(pos_vel(1,1:3),pos_vel(1,4:6)+[0,0,0.1],station2,times);
	r=(range1-range2);
	b(:,6)=((r.-mean(r))-baseline).*10;
	error_jacobian=b;
endfunction
