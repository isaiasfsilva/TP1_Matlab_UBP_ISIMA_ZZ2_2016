function [f,g]=fctgeneral2(x)
%Deuxième fonction pour tester la robustesse du code
%Pour x0 = (1,2,...,n)^T
    n=length(x);
    f=0;
%   
    f = 0.00001*(sum((x-1).^2))  + (norm(x)^2 - 0.25)^2;  
%      
    if(nargout==1)
       return;
    end
    g=zeros(n,1);
%     
    common_term = 4*norm(x)^2 - 0.99998 ;
    g = x * common_term - 0.00002; %optization
end   
   
