function [xK,F,k]=bfgs(myfct, xK,eps, iprint)
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

if(length(xK)>100)
   fprintf('\n#\tERROR! TAILLE DE X PLUS GRAND QUE 100\n#\tDesolé, mais le taille du X ne peut pas être plus grand qui 100 éléments\n#\n'); 
   xK=0; F=0; k=0;
   return;
end
if(nargin==3)
        iprint=0;
end
k=0;
%xK=[0;0];
n=length(xK);
[F,gK]=feval(myfct,xK);

H=eye(n); 
t=0;
iterations = 1 ;
while(norm(gK) > eps && k<100)
        [R, p] = chol(H); % Complexite O(n^3) - Pas bonne :/
	if (p == 0) % H est définie positive
		dK = -1* H * gK; %Dérivative
	end %if 
	
	[t, armijo_iter]=rarmijo(myfct,F,gK,dK,xK); %Recherche linear avec Armijo méthod
	
	if(armijo_iter==100)
	        fprintf('#Recherche linéaire est mal passé. Désolé. \n\n');
                return;
	end
	
	xK=xK+t*dK; %Calcul of X
	
	[F,gNextK]=feval(myfct,xK); 
	
	deltK=t*dK; %ok
	
	gamK=gNextK - gK; % OK

	if(all(gamK~=0) && all(deltK~=0)) % Logical and
		H = H - (1/(deltK'*gamK))*(deltK*gamK'*H + H*gamK*deltK') + (1 + (gamK'*H*gamK)/(deltK'*gamK)) * ((deltK*deltK')/(deltK' * gamK));
	end %if
	
	if(k==0 && iprint==1)
	     fprintf('\n\n#### Value de Gradient (gK) et X temporaire ####\n#          Premier Iteration          #\n#####################################\n');
	     gK
	     xK	     
        end
        
	if(iprint==2)
	     fprintf('\n\n#### Value de Gradient (gK) et X temporaire ####\n#            Iteration %d          #\n#####################################\n',k);
	     gK 
	     xK	     
        end
        iterations = iterations + armijo_iter+1;
        
        %Vérifie si le numéro d'iterations est plus grand que 5 fois le numéro maximum d'iterations.
        if(iterations>500)
                fprintf('#Numero d"evaluations plus grand que 5 fois le numéro maximum d"iterations ( 500)\n#            ABORTING          \n\n');
                return;
        end
	k=k+1; %OK
	gK=gNextK;
end


if(iprint==1)
	   fprintf('\n\n#### Value de Gradient (gK) et X temporaire ####\n#           Dérniére Iteration          #\n#####################################\n');
	   gK
	   xK
	     
end
        
        
