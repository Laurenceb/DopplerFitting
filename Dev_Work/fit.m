function state_average=fit()
CUTOFF=0.001;					#reduction in uncertainty we require
N=[0,0,0,0;						#noise generation matrix
   0,0,0,0;
   0,0,50,0;
   0,0,0,50];
load("-binary","dopplers","dopplers");
dopplers=dopplers(41:size(dopplers)(2));
dopplers+=0.1*randn(size(dopplers));#noise equivalent to over 50 degrees on a PLL
########CONFIGURATION#######################################################
STATION_POS=[0;6370000];				#on the surface
DELTA_T=0.05;						#100milliseconds
DOPPLER_T=1;						#15 second intervals
TIMEOUT=50;
BURN_IN_LIMIT=100;
DEBUG=0;
BURN_IN_ERROR=2;
###########################################################################
counter=0;
iterations=0;
burnt_in=0;
reduction_factor=BURN_IN_ERROR/N(4,4);
f=1;
#state=[STATION_POS+[150000;90000];2000;1000];		#massive guesstimate of position
#initial_state=[ 1.1623e+05;6.4643e+06;2720.7;1272.8];
initial_state=[1.2e+05;6.47e+06;2700;1300];
#state=[1.18e+05;6.44e+06;2500;1400];
#state(3:4)=[2.8158e+03;1.6492e+03];
state=initial_state;
state_average=initial_state;
state_error=N(4,4);
old_q=simulate_ballistic(state,dopplers,STATION_POS,DELTA_T,DOPPLER_T);#initialise the q value
states=state;
timeout=TIMEOUT;
while(!(iterations>BURN_IN_LIMIT && timeout==TIMEOUT) && state_error>0.01)#main loop
	counter+=1;
	do	
		new_state=state+(N*randn(4,1)*f);	#we are just trying the basic 2D case for the moment - i.e. state vector is 4D
		errors=new_state-initial_state;
	until(abs(errors(1))<=N(1,1) && abs(errors(2))<=N(1,1) && abs(errors(3))<=N(3,3) && abs(errors(4))<=N(4,4))#domain of solutions
	q=simulate_ballistic(new_state,dopplers,STATION_POS,DELTA_T,DOPPLER_T);#try the new state
	a=old_q/q;
	if(a>=1)					#accept the state
		state=new_state;
		old_q=q;
		f*=1.1;				#acceptance rate scales the noise
	else
		if(3*rand<a)				#accept it with a probability proportional to a
			state=new_state;
			old_q=q;
			f*=1.1;
		else
			f*=0.9;
		endif	
	endif
	iterations+=1;
	states=[states,state];
	if(q<1)
		state;
	endif
	if(f>1)
		f=1;
	endif
	if(DEBUG)	
		printf("noise:%f  mean square error:%f v_up:%f\n",f,old_q,state(4));
		fflush(stdout);
	endif
	if((N(4,4)*f<sqrt(old_q)))
		counter+=1;
		if(counter>10)
			f=10*sqrt(old_q)/N(4,4);
		endif
	else
		counter=0;
	endif
	if(burnt_in)
		state_average+=(new_state-state_average)*state_error/(q+state_error);
		state_error*=(1-(state_error/(q+state_error)));
	endif
	if(f<reduction_factor)
		burnt_in=1;
		timeout-=1;
	endif	
endwhile
plot(states(3,:),states(4,:));				#2D plot of the progress
state=state_average;
endfunction

	
#simulates a ballistic trajectory up to the apogee
function meansquareoffset=simulate_ballistic(state_initial,dopplers,station_position,delta_t_model,delta_t_doppler)
	n=size(dopplers)(2);				#this is the number of doppler (normalised to relative speed) values
	position=state_initial(1:2);			#setup values for sim 
	velocity=state_initial(3:4);
	nopoints=n;
	time=0;
	meansquareoffset=0;
	while(n)
		if(time>delta_t_doppler || time==0)
			time=0;
			n-=1;
			offset=position-station_position;
			diff=(dopplers(nopoints-n)-dot((offset/norm(offset)),velocity))^2;
			meansquareoffset+=diff;
			fflush(stdout);
		endif
		accelg=-59736*6673*10^6*(norm(position)^-3);#assume ballistic flight, no earth rotation as yet
		velocity+=accelg*position*delta_t_model;
		position+=velocity*delta_t_model;
		time+=delta_t_model;
	endwhile
	meansquareoffset/=nopoints;
endfunction
