function [Y] = c_theta(R0, Rt, Rc, k, me, mu, theta)

    % fonction contraintes pour le problème de trajectoire pour Ariane 6 (qui fournit les angles theta)

    % Input
    % R0 : Position initiale du lanceur (Vecteur de R^2)
    % Rt : Rayon de la Terre
    % Rc : Hauteur à atteindre (rayon de la Terre + hauteur de l'orbite
    %       cible)
    % k : Constantes constructeurs k1, k2, k3 (Vecteur de R^3)
    % me : vecteur de R^3, contenant les masses d'ergos me1, me2, me3 
    % mu : Masse du satellite
    % theta : Angles de poussée (Vecteur de R^4)

    % Output
    % Y : vecteur de R^2, valeur des contraintes

    global N_call_c_test;

    [R,V] = simulateur(R0, Rt,  k, me, mu,  theta);
    
    normR = norm(R(4, :));
    normV = norm( V(4, :) );
    
    Y(1) = ( normR - Rc ) / Rc ; 
    Y(2) = dot(  R(4, :),  V(4, :) ) / (normR*normV);
    Y = Y';
    N_call_c_test = N_call_c_test + 1;

end

