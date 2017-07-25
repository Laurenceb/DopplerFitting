#quick test code - grid velocity
#for n=1:21
#	for m=1:21
#		vect=[(n-11)*20,(m-11)*20,0];
#		vel=posvel(1,4:6).+vect;
#		[dop,range1]=generate_doppler(posvel(1,1:3),vel,station1,times);
#		[dop,range2]=generate_doppler(posvel(1,1:3),vel,station2,times);
#		r=(range1-range2);		
#		r.-=mean(r);
#		delta=observed_delta_range_.-r;
#		difference(n,m)=max(abs(delta));
#	end
#	printf("%d\r",n);fflush(stdout);
#end 
#return;

#return;
	################ Course velocity fitter for less than 3 stations ####
#	verticals=[-300:100:200];#one hudred increments over potential range
#	horizontals=[-4,0;-3,0;-2,0;-1,0;0,0;1,0;2,0;3,0;4,0;-4,-1;-3,-1;-2,-1;-1,-1;0,-1;1,-1;2,-1;3,-1;4,-1;-3,-2;-2,-2;-1,-2;0,-2;1,-2;2,-2;3,-2;-2,-3;-1,-3;0,-3;1,-3;2,-3;-1,-4;0,-4;1,-4;-4,1;-3,1;-2,1;-1,1;0,1;1,1;2,1;3,1;4,1;-3,2;-2,2;-1,2;0,2;1,2;2,2;3,2;-2,3;-1,3;0,3;1,3;2,3;-1,4;0,4;1,4];#approximate circle around zero
#	v_cylinder=[horizontals.*100,zeros(size(horizontals,1),1)];
#	v_cylinder=repmat(v_cylinder,size(verticals,2),1);
#	v_cylinder(:,3).+=repmat(-verticals,size(horizontals,1),1)(:);#add on vertical v, repeating elements for each set of horizontal velocities
##	v_cylinder=xform( 'vel', v_cylinder, 1, 0 );
#	for n=1:size(v_cylinder,1)	#loop through calculating weighted delta error
#		vel=posvel(1,4:6).+v_cylinder(n,:);
#		[dop,range1]=generate_doppler(posvel(1,1:3),vel,station1,times);
#		[dop,range2]=generate_doppler(posvel(1,1:3),vel,station2,times);
#		r=(range1-range2);		
##		r.-=mean(r);
##		delta=observed_delta_range_.-r;
#		difference(n)=max(abs(delta));#mean((delta).^2);
#printf("%f\n",(difference(n)));fflush(stdout);
#	endfor
#	[a,b]=min(difference);
#	v_offset=xform( 'vel', v_cylinder(b,:), 0, 1 );
#	printf("Found minima of %f (m squared) at quantised offset %f,%f,%f\n",a,v_offset(1),v_offset(2),v_offset(3));fflush(stdout);
#	posvel(1,4:6)=vel;#.+=v_cylinder(b,:);
	#posvel(1,4:6).+=[150,-100,0];
