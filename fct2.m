function [f,g]=fct2(x)
%Deuxième fonction pour la mise au point du code
%Pour x0 = [(7/6)^1(1/2);0] la soluction optimale est x* = [0;0] ou x* = [(7/12)^1(1/2);-(7/12)^1(1/2)]
    x1=x(1);
    x2=x(2);
%
	f = (x1 + x2)^2 + (2*(x1^2 + x2^2 -1) - 1/3)^2;	
	if(nargout==1)
		return;
	end
	g=zeros(2,1);
	g(1)=16*x1^3 + 16*x1*x2^2 - (56/3)*x1 + 2*x2+2*x1;
	g(2)=16*x2^3 + 16*x2*x1^2 - (56/3)*x2 + 2*x1+2*x2;
end	
