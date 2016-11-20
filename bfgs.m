function [xK,F,k]=bfgs(fct, xK,eps, itr, iprint)
% Algorithme BFGS pour résoudre le méthode de quasi-Newton
%SINTAXE:
%[xK,F,k]=bfgs(myfct, xK,eps, itr, iprint)
%
%La fonction reçois comme parâmetres le système a résoudre, le x initial et la tolérance.Retourne le vecteur soluction du système,
%la valeur du système pour la soluction trouvé et le nombre d'itérations effectués.
%La fonction retourne le x obtenu pendant l'exécution de l'algorithme même sans être le optimal, car il y a possibilité d'échecs.
%Il faut considerer que a chaque itération il y a la construction d'une solution en essayant obtenir la solution optimale, alors
%chaque itération retourne une soluction faisable.
%
%ARGUMENTS:
%fct	fonction (système) sur laquelle s'applique le méthode
%xK	valeur de X initial
%eps	tolérance
%itr	nombre maximum d'itérations permis
%iprint	(optionnel) indicateur d'affichage du gradient de la fonction et du vecteur solution
%	par défaut:	0 -> Rien à afficher
%			1 -> Affichage des première et dernière itérations
%			2 -> Affichage par toutes les itérations
%
%SORTIE:
%xK	vecteur solution xK
%F	valeur de la fonction en xK
%k	nombre d'itérations exécutées
%
%INFORMATIONS ADDITIONNELS
%Cette fonction permet le traitement pour matrices dont la taille maximum est fixée à 100. Le contrôle de performance du méthode
%est fait pour le contrôle de la norme du gradient par rapport a une tolerance passée comme paramètre et du nombre d'itérations.
%Aussi on se limite aux matrices définies positives pour modifier le valeur des paramètres lorsque l'algorithme BFGS est éxécutée.
%La vérification de la matrice comme définie positive est faite avec le paramètre optionnel de la fonction chol, si la valeur est 
%nulle, alors la matrice est définie positive. Cependant, cette approche est trés coûteuse, car dans le pire cas peut-être O(n³),
%pour une matrice carée de taille n. La performance a été optimisée pour être plus efficace et peu coûteuse, et seulement le cas
%de la vérification d'être définie positive que n'est pas suffisant pour cet objectif.
%L'exécution de l'algoritme est aussi submise au contrôle d'itérations et d'évaluation. S'il y a un nombre d'évaluations de la 
%fonction et son gradient supérieur a 5 fois le nombre maximum d'itérations l'exécution doit être arrétée.
%
%
%L'initialisation des variables que sont utilisées pendant l'exécution de l'algorithme
%Compteur d'itérations
k=0;
%La taille du vecteur passé comme paramètre
n=length(xK);
%Vérification de la taille du vecteur passé comme paramètre. Il faut qu'elle soit au plus 100, sinon, il est interdit exécuter 
if(n > 100)
   fprintf('\n#\t ÉCHEC! TAILLE DE X PLUS GRAND QUE 100\n#\t Il faut que la taille du vecteur X soit au plus 100! \n#\n'); 
  % xK=0; F=0; k=0; // acredito que apenas imprime a mensagem de erro, a gente não tem que alterar o valor do x e etc
   return;
end
%Vérification de la commande d'affichage. Si le nombre de arguments est 4, il n'a aucune besoin d'afficher les valeurs lorsque
%l'algorithme est exécutée. Alors, l'attribut iprint reçoit la valeur par défaut.
if(nargin == 3)
        iprint = 0;
end
%Évaluation de la fonction et de son gradient
[F,gK]=feval(fct,xK);
%Initialisation de la matrice pour sauvegarder l'approximation de l'inverse du Hessien utilisé pendant l'exécution du méthode BFGS
H=eye(n);
%Variable 
t=0;
iterations = 1 ;
while(norm(gK) > eps && k<itr)
	if(iterations> 5 * itr)
                fprintf('#Numero d"evaluations plus grand que 5 fois le numéro maximum d"iterations\n#            ABORTING          \n\n');
                return;
        end
        [R, p] = chol(H); % Complexite O(n^3) - Pas bonne :/
	if (p == 0) % H est définie positive
		dK = -1* H * gK; %Dérivative
	end %if 
	
	[t, armijo_iter]=rarmijo(fct,F,gK,dK,xK); %Recherche linear avec Armijo méthod
	
	if(armijo_iter==100)
	        fprintf('#Recherche linéaire est mal passé. Désolé. \n\n');
                return;
	end
	
	xK=xK+t*dK; %Calcul of X
	%Évaluation de la fonction et de son gradient pour le prochain vecteur xK
	[F,gNextK]=feval(fct,xK); 
	
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
        
	k=k+1; %OK
	gK=gNextK;
end


if(iprint==1)
	   fprintf('\n\n#### Value de Gradient (gK) et X temporaire ####\n#           Dérniére Iteration          #\n#####################################\n');
	   gK
	   xK
	     
end
        
        
