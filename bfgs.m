function [xK,F,k]=bfgs(fct, xK,eps, itr, iprint)
% Algorithme BFGS pour résoudre le méthode de quasi-Newton
%SINTAXE:
%[xK,F,k]=bfgs(myfct, xK,eps, itr, iprint)
%
%La fonction reçoit comme paramètres le système à résoudre, le x initial et la tolérance.Retourne le vecteur solution du système,
%la valeur du système pour la solution trouvée et le nombre d'itérations effectués.
%La fonction retourne le x obtenu pendant l'exécution de l'algorithme même sans être optimal, car il y a possibilité d'échecs.
%Il faut considerer qu' à chaque itération il y a la construction d'une solution en essayant obtenir la solution optimale, alors
%chaque itération retourne une solution réalisable.
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
%QUELQUES EXAMPLES
%[xout, f,k] = bfgs(@fctgeneral,ones(40,1),0.00000001,500,0)
%
%[xout, f,k] = bfgs(@fctgeneral2,[1:90]',0.00001,500,0)
%
%[xout, f,k] = bfgs(@fctgeneral2,[1:10]',0.00001,500,0)
%
%[xout, f,k] = bfgs(@fct,[0,0]',0.00001, 500,0)
%
%[xout, f,k] = bfgs(@fct2,[(7/6)^0.5,0]',0.00001, 500,0)
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
   error('ÉCHEC! TAILLE DE X PLUS GRAND QUE 100 - Il faut que la taille du vecteur X soit au plus 100!'); 
  % xK=0; F=0; k=0; // acredito que apenas imprime a mensagem de erro, a gente não tem que alterar o valor do x e etc
   return;
end
%Vérification de la commande d'affichage. Si le nombre de arguments est 4, il n'a aucune besoin d'afficher les valeurs lorsque
%l'algorithme est exécutée. Alors, l'attribut iprint reçoit la valeur par défaut.
if(nargin == 4)
        iprint = 0;
end
%Évaluation de la fonction et de son gradient
[F,gK]=feval(fct,xK);
%Initialisation de la matrice pour sauvegarder l'approximation de l'inverse du Hessien utilisé pendant l'exécution du méthode BFGS
H=eye(n);
%Variable pour enregistrer le retourn de la recherche linéaire.
t=0;
%Variable pour contrôler le nombre d'évaluations
iterations = 1 ;
%L'algorithme pour le méthode BFGS
%Tandis que le norme du gradient de la fonction est supérieur a la tolerance et le nombre d'iterations est inférieur au maximum
%l'algorithme est exécutée.
%
while(norm(gK) > eps && k<itr)
%Vérification de la contrainte d'évaluations. Il est interdit d'avoir plus évaluations que 5 fois le nombre d'itérations
	if(iterations> 5 * itr)
                error('ÉCHEC! NOMBRE D EVALUATIONS SUPERIEUR AU PERMIS - Il y a plus évaluations que 5 fois le nombre maximum d itérations         ABORTING        ');
                return;
        end %if contrôle d'évaluations
%
%Vérification de la matrice H comme définie positive. Si p == 0 alors la matrice est définie positive. 
%En fait la fonction utilisée a le rôle de factoriser une matrice, mais sa sortie optionnel donne directement la définition de la
%matrice comme définie positive dans le cas où le paramètre optionnel p est nul. Comme dit avant, c'est très coûteux.
%
        [R, p] = chol(H); 
	if (p ~= 0) 

%On actualise la valeur de dK si et seulement si on a une matrice définie positive
        error('ÉCHEC! LA MATRICE N EST PAS DEFINIE POSITIVE. IL FAUT ARRETER L EXECUTION.    ABORTING     ')		 
	end %if vérification d'être definie positive et d'actualisation de la dérivée
%dk -> variable pour enregistrer le gradient de la fonction approché pour le Hessien
       dK = -1* H * gK;
%
%Appel à fonction que décrire la règle d'Armijo pour la recherche linéaire et donne la valeur du terme pour améliorer la solution	
	[t, armijo_iter]=rarmijo(fct,F,gK,dK,xK);
%Les sorties sont le terme de la recherche linéaire (t) et le nombre d'évaluations faites pendant l'exécution de la fonction
%
%Vérification de la performance de la recherche linéaire. Si le nombre d'iterations dans la recherche linéaire est supérieur a
%100 alors la recherche est mal passée. Alors il faut finir l'algorithme a cause de n'avoir pas comme obtenir une solution faisable
	if(armijo_iter==100)
	        error('ÉCHEC! NOMBRE D ITÉRATIONS DANS LA RECHERCHE LINÉAIRE SUPERIEUR AU PERMIS - Recherche linéaire est mal passée.');
                return;
	end % if de vérification de la performance de la recherche linéaire
%
%Calcul du nouveau vecteur solution par rapport a la recherche linéaire et le gradient de la fonction sur la solution actuelle
	xK=xK+t*dK;
%Évaluation de la fonction et de son gradient pour le nouveau vecteur xK
	[F,gNextK]=feval(fct,xK); 
%
%Calcul des paramètres pour la définition du Hessien.
%Différence entre la nouvelle solution et la solution actuelle
	deltK=t*dK; 
%Différence entre les gradients de la nouvelle solution et de la solution actuelle
	gamK=gNextK - gK;
%
%Calcul du Hessien approché sur la contrainte que les solutions sont différents. C'est à dire, on peut améliorer la solution
%en recherchant la solution optimale
%Le calcul du Hessien est donné pour une formule obtenue dans toutes les documents que traitent du méthode BFGS
	if((deltK'*gamK)>0)
	        %error('ÉCHEC! LA MATRIX NE PEUT ETRE CALCULEE. IL Y A EU UNE DIVISION PAR ZERO.       ABORTING      ');
                H = H - (1/(deltK'*gamK))*(deltK*gamK'*H + H*gamK*deltK') + (1 + (gamK'*H*gamK)/(deltK'*gamK)) * ((deltK*deltK')/(deltK' * gamK));
	end %if calcul du Hessien pour nouvelles solutions
	
	if(mod(k,n)==0)
	        H=eye(n);
	end

%Vérification des specifications d'affichage pour la première itération dans le cas du paramètre iprint égal a 1
%Si l'algorithme est dans la première itération et comme il y a besoin d'affichage des valeurs du gradient et de la solution
% alors cettes valeurs sont affichées
	if(k==0 && iprint==1)
	     fprintf('\n\n#### Value de Gradient (gK) et X  ####\n#          Premier Iteration          #\n#####################################\n');
	     gK
	     xK	     
        end  %fin de la vérification du paramètre d'affichage iprint == 1
%
%Vérification des specifications d'affichage pour la itération courante dans le cas du paramètre iprint égal a 2
%Comme il y a besoin d'affichage des valeurs du gradient et de la solution alors cettes valeurs sont affichées        
	if(iprint==2)
	     fprintf('\n\n#### Value de Gradient (gK) et X  ####\n#            Iteration %d          #\n#####################################\n',k);
	     gK 
	     xK	     
        end  %fin de la vérification du paramètre d'affichage iprint == 2
%
%Calcul du nombre d'évaluations de la fonction et de son gradient
        iterations = iterations + armijo_iter+1;
%Augmentation dans le compteur d'itérations
	k=k+1;
%Attribution du gradient de la nouvelle solution a la variable correct car gNext est une variable temporaire
	gK=gNextK;
end %fin du boucle while

%Vérification des specifications d'affichage pour la dernière itération dans le cas du paramètre iprint égal a 1
%Comme l'algorithme est dans la dernière itération et comme il y a besoin d'affichage des valeurs du gradient et de la solution
% alors cettes valeurs sont affichées  
if(iprint==1)
	   fprintf('\n\n#### Value de Gradient (gK) et X  ####\n#           Dérniére Iteration          #\n#####################################\n');
	   gK
	   xK
	     
end %fin de la vérification du paramètre d'affichage iprint == 1
%
%
end % fin de la fonction
        
        
