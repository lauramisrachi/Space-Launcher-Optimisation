function [Delta_V] = c_etagement(X, mu, k, ve, Vp)

    % Fonction contraintes pour le problème d'étagement d'Ariane 6

    % Input
    % X vecteur de R^3, contenant les masses d'ergos me1, me2, me3 
    % mu : Masse du satellite
    % k : Constantes constructeurs k1, k2, k3

    % Output
    % Delta_V, valeur réelle des contraintes

        global N_call_c_test;

        Mi1 = mu + X(1) * (1 + k(1)) + X(2) * (1 + k(2)) + X(3) * (1 + k(3)) ; 
        Mf1 = Mi1 - X(1);

        Mi2 = Mf1 - k(1)*X(1);
        Mf2 = Mi2 - X(2);

        Mi3 = Mf2 - k(2)*X(2);
        Mf3 = Mi3 - X(3);

        Delta_V = (ve(1) * log(Mi1/Mf1) +  ve(2) * log(Mi2/Mf2) +  ve(3) * log(Mi3/Mf3) - Vp);


        N_call_c_test = N_call_c_test + 1;

end

