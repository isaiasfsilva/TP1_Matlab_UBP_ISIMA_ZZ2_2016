function [xK,F,k]=bfgs(myfct, xK,eps)
% Algorithme BFGS pour résoudre le méthode de quasi-Newton
%SINTAXE:
%[xK,F,k]=bfgs(myfct, xK,eps)
%
%La fonction reçois comme parâmetres le système a résoudre, le x initial et la tolérance.
%Elle retourne le vecteur soluction du système, la valeur du système pour la soluction trouvé et la quantité d'itérations necessaires.
%La quantité maximale d'itérations est de 100, alors, la fonction retourne le x obtenu, même sans être le optimal.
%
%ARGUMENTS:
%fct	fonction (système) sur laquelle s'applique le méthode
%xK	valeur de X initial
%eps	tolérance
%
%SORTIE:
%xK	vecteur solution
%F	valeur de la fonction en xK
%k	quantité d'itérations
k=0;
%xK=[0;0];
n=length(xK);
[F,gK]=feval(myfct,xK);%X est une matriz
%eps=0.001;

%gK = -dK;
H=eye(n); 

while(norm(gK) > eps && k<1000)
        [R, p] = chol(H);
	if (p == 0) % H est définie positive
		dK = -1* H * gK;
	end %if 
	
	t=rarmijo(myfct,F,gK,dK,xK); 
	
	
	xK=xK+t*dK; %OK
	
	[F,gNextK]=feval(myfct,xK); 
	deltK=t*dK; %ok
	gamK=gNextK - gK; % OK

	if(sum(deltK == 0)==0 && sum(gamK==0)==0)
		H = H - (1/(deltK'*gamK))*(deltK*gamK'*H + H*gamK*deltK') + (1 + (gamK'*H*gamK)/(deltK'*gamK)) * ((deltK*deltK')/(deltK' * gamK));
	end %if
	
	k=k+1; %OK
	gK=gNextK;
end
