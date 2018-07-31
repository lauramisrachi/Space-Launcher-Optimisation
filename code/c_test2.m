function [Delta_V] = c_test2(X)

    % Fonction contrainte pour le problème d'étagement d'Ariane 1

    % Input
    % X vecteur de R^3, contenant les masses d'ergos me1, me2, me3 

    % Output
    % Delta_V : valeur réelle des contraintes

        global N_call_c_test;

        % Paramètres

        k1 = 0.1101 ; 
        k2 = 0.1532 ; 
        k3 = 0.2154 ; 
        mu = 1700 ; 
        ve1 = 2647.2 ; 
        ve2 = 2922.4 ; 
        ve3 = 4344.3 ; 
        Delta_vn = 11527 ; 

        Mi1 = mu + X(1) * (1 + k1) + X(2) * (1 + k2) + X(3) * (1 + k3) ; 
        Mf1 = Mi1 - X(1);

        Mi2 = Mf1 - k1*X(1);
        Mf2 = Mi2 - X(2);

        Mi3 = Mf2 - k2*X(2);
        Mf3 = Mi3 - X(3);

        Delta_V = ve1 * log(Mi1/Mf1) +  ve2 * log(Mi2/Mf2) +  ve3 * log(Mi3/Mf3) - Delta_vn;

        N_call_c_test = N_call_c_test + 1;

end

