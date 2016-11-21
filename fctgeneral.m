function [f,g]=fctgeneral(x)
%Première fonction pour tester la robustesse du code
%Pour x0 = (1;1;...;1)^T
        n=length(x);
        f=0;
        x1 = x(2:n);
        x2 = x(1:n-1);
        i = 2:n;
        f = sum(i* (x1.*2 - x2) .^2);
%
        if(nargout==1)
	        return;
        end
        
        g=zeros(n,1);
        
        i = [1:n]';
        g = 8 * i .* vertcat(0,x1) -4*i.*vertcat(0,x2) -4*vertcat(i(2:n),0) .* vertcat(x1,0) + 2* vertcat(i(2:n),0) .* x;
end
	
