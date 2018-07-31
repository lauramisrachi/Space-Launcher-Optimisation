
function [Y] = fonction_merite_opt(f_Xk, c_Xk, rho, lambda)
   
    % fonction mérite optimisé (i.e. sans appel à f ni c) spécialisé pour le problème min f(x) s.c. c(x) = 0
    % Input :
    % f_Xk : valeur de la fonction objectif en Xk 
    % c_Xk : valeur de la fonction contrainte en Xk
    % rho : Paramètre  de pénalisation
    % lambda : Multiplicateur de Lagrange du problème
    % Output :
    % Y : Valeur de la fonction mérite en X
    
    Y = f_Xk + rho * norm(c_Xk) ; 

end