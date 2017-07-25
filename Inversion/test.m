a=fit();
for i=1:20
	a=[a,fit()];
	i
	fflush(stdout);
endfor
plot(a(1,:),a(2,:));
