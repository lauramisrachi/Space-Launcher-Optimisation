function [V ] = f_theta( R0, Rt, k, me, mu, theta)

    % Fonction objectif pour le problème de trajectoire pour Ariane 6 (qui fournit les angles theta)
    % Input
    % R0 : Position initiale du lanceur (Vecteur de R^2)
    % Rt : Rayon de la Terre
    % k : Constantes constructeurs k1, k2, k3 (Vecteur de R^3)
    % me : vecteur de R^3, contenant les masses d'ergos me1, me2, me3 
    % mu : Masse du satellite
    % theta : Angles de poussée (Vecteur de R^4)

    % Output:
    % V : Opposé de la norme du vecteur vitesse final

        global N_call_f_test;

        [~,V] = simulateur(R0, Rt,  k, me, mu,  theta);
        V = V(4, :);
        V = -norm(V) * 1E-3 ; % * 1E-3
        N_call_f_test = N_call_f_test + 1;

end
