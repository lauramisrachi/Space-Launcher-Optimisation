function [Y] = fonction_merite(X, f, c, rho, lambda)
   
    % Fonction mérite spécialisé pour le problème min f(x) s.c. c(x) = 0
    % Input :
    % X : Point d'application
    % f : fonction objectif 
    % c : Fonction contrainte 
    % rho : Paramètre  de pénalisation
    % lambda : Multiplicateur de Lagrange du problème
    % Output :
    % Y : Valeur de la fonction mérite en X
    
    Y = f(X) + rho * norm(c(X)) ; 

end
