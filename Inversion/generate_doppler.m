function [doppler,range,positions]=generate_doppler(position,velocity,station,times_)
	#position, station, and velocity are in ECEF units, as column vectors
	#note that this ignores high order gravity effects, but this should only give <10m error
	#lunar gravity is <1m error
	#the rocket config, some of this could be made iterative/callable
	mass=0.048;
	radius=0.0122;
	drag_coefficient=0.2;#0.151;
	#########
	dt=diff(times_)(1);
	n=1;
	mass_earth=5.972e24;
	radius_earth=6370e3;
	ge=6.674e-11;
	we=7.292e-5;
	light_speed=2.9979e8;
	velocity+=[-position(2),position(1),0].*we;#Add on velocity from Earth rotation
	#########
	range=zeros(1,length(times_));
	doppler=zeros(1,length(times_));
	positions=zeros(length(times_),3);
	for t=times_
		x=[position,velocity];
		velocity_air=[-position(2),position(1),0].*we;
		la=xform( 'pos', position, 0, -1 );
		density=US_standard_atmosphere(la(3));
		if(norm(x(4:6))>320 && norm(x(4:6))<500)
			density*=2;
		endif
		g_vector=-x(1:3)*ge*mass_earth/(norm(x(1:3))^3);
		g_vector-=(x(4:6)-velocity_air).*(pi*radius^2*0.5*drag_coefficient*density*norm(x(4:6)-velocity_air)/mass);
		k_1=[x(4:6),g_vector].*dt;

		x_1=x+k_1/2;
		g_vector=-x_1(1:3)*ge*mass_earth/(norm(x_1(1:3))^3);
		g_vector-=(x_1(4:6)-velocity_air).*(pi*radius^2*0.5*drag_coefficient*density*norm(x_1(4:6)-velocity_air)/mass);
		k_2=[x_1(4:6),g_vector].*dt;

		x_2=x+k_2/2;
		g_vector=-x_2(1:3)*ge*mass_earth/(norm(x_2(1:3))^3);
		g_vector-=(x_2(4:6)-velocity_air).*(pi*radius^2*0.5*drag_coefficient*density*norm(x_2(4:6)-velocity_air)/mass);
		k_3=[x_2(4:6),g_vector].*dt;

		x_3=x+k_3;
		g_vector=-x_3(1:3)*ge*mass_earth/(norm(x_3(1:3))^3);
		g_vector-=(x_3(4:6)-velocity_air).*(pi*radius^2*0.5*drag_coefficient*density*norm(x_3(4:6)-velocity_air)/mass);
		k_4=[x_3(4:6),g_vector].*dt;

		position+=(k_1(1:3)+2*k_2(1:3)+2*k_2(1:3)+k_4(1:3))/6;
		velocity+=(k_1(4:6)+2*k_2(4:6)+2*k_2(4:6)+k_4(4:6))/6;
		#transfer from ECEF into our inertial frame
		position_station=([cos(we*t),-sin(we*t),0;sin(we*t),cos(we*t),0;0,0,1]*station')';
		velocity_station=[-position_station(2),position_station(1),0].*we;
		range(n)=norm(position_station-position);
		doppler(n)=dot((velocity_station-velocity),(position_station-position))/(light_speed*norm(position_station-position));
		positions(n,:)=position;
		n+=1;
	endfor
endfunction


function density=US_standard_atmosphere(altitude)
	if(altitude>9e4)
		density=0;
	else
		density=stdatmo(altitude,10);#+10C temperature offset
	endif
endfunction
