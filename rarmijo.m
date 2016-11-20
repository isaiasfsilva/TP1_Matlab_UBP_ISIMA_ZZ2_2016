function [t, iter]=rarmijo(fct, fk, gk, dk, xk)
%Fonction pour la recherche linéaire à partir de la règle d'Armijo. Prend t que vérifie
%la contrainte de la recherche linéaire définie dans la règle. Le t est le pas pour lequel
%on peu approcher de la solution optimale de façon accpelérée.
%
%SINTAXE:
%[t, iter]=rarmijo(fct, fk, gk, dk, xk)
%
%ARGUMENTS:
%fct	fonction (système) sur laquelle s'applique le méthode
%fk	valeur de X initial
%gk	gradient de xk
%dk	gradient de xk approché pour le Hessien
%xk	valeur du vecteur x
%
%SORTIE:
%t	vecteur solution xK
%iter	nombre d'évaluations du vecteur et de son gradient
%
%INFORMATIONS ADDITIONNELS
%Pour construire le pas il y a des paramètres tau et sigma dans un intervalle prefixé.
%Pour cette fonction ils sont définis comme tau égal a 0.5 et sigma égal a 0.001.
%
%
%L'initialisation des variables utilisées pendant la recherche linéaire pour modifier t
tau=0.5;
sigma=0.001;
%Initialisation du pas comme 1
t=1;
%initialisation du compteur d'évaluation comme nul
iter=0;
%Boucle de la recherche linéaire limité au maximum de 100 itérations
while(iter<100)
%Le valeur de la solution x est calculé vers l'inverse de son gradient approché pour le
%Hessien à partir du pas t
	x=xk+t*dk;
%Évaluation de la fonction
	f=feval(fct,x);
%Vérification de la condition pour modifier le pas dans la recherche linéaire d'accord
%la règle d'Armijo
	if(f-fk <=t*sigma*gk'*dk)
%Si la valeur de la fonction sur laquelle est cherché le pas est adequée, le pas est optimale
		return;
	end %end vérification de la contrainte de la règle d'Armijo
%Augmentation du pas vers la solution optimale
	t=tau*t;
%Augmentation du compteur d'évaluations
	iter=iter+1;
end; %fin de la fonction
