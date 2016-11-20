disp('_____________________________________________________________')
disp('|                        TP1 Matlab                           |')
disp('| Université Blaise Pascal - ISIMA - 2016/2 - ZZ2 - Filière 4 |')
disp('|                                                             |')
disp(' -------------------------------------------------------------')
disp(' ')
disp(' >> Étudiants:')
disp('    - Isaías FARIA SILVA')
disp('    - Paula METZKER SOARES')
disp('  ')
  
disp('  >> Novembre - 2016 - All rights reserved.')
fprintf('\n\n Salut, Je vais executer les examples suivantes:\n\n');

disp('>> Exemple Basique 1 : ')
X0 = [0,0]'
[xout, f,k] = bfgs(@fct,X0,0.00001, 500,0)

disp('>> Exemple Basique 2 : ')
X0 = [(7/6)^0.5,0]'
[xout, f,k] = bfgs(@fct2,X0,0.00001, 500,0)

disp('>> Exemple Robustesse 1: ')
disp('X0 = [1, 1, 1, ... ,1 ,1 ,1 ] %taille 40')
[xout, f,k] = bfgs(@fctgeneral,ones(40,1),0.00000001,500,0)

disp('>> Exemple Robustesse 2: ')
disp('X0 = [1, 2, 3, ... ,n-1 ,n ] %n=90')
[xout, f,k] = bfgs(@fctgeneral2,[1:90]',0.00001,500,0)

disp('Et Voilà, J"espère qui tous si pass bien')
