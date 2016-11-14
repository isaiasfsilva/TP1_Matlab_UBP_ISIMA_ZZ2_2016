function [f,g]=fct(x)
	f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;	
	if(nargout==1)
		return;
	end
	g=zeros(2,1);
	g(1)=100*(-4*x(2)*x(1) + 4 * x(1)^3) +2 * x(1) -2;
	g(2)=200*(x(2)-x(1)^2);
	
