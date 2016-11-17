function t=rarmijo(fct, fk, gk, dk, xk)
tau=0.5;
sigma=0.001;
t=1;
iter=0;
while(iter<100)

	x=xk+t*dk;
	f=feval(fct,x);  %fct(x)
	if(f-fk <=t*sigma*gk'*dk)		
		return;
	end
	t=tau*t;
	iter=iter+1;
end;
