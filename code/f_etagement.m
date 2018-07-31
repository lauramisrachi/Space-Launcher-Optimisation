function [ M_i1 ] = f_etagement( X, mu, k )

        % Fonction objectif pour le problème d'étagement d'Ariane 6
        % Input:
        % X : vecteur de R^3, contenant les masses d'ergos me1, me2, me3
        % mu : Masse du satellite
        % k : Constantes constructeur k1, k2, k3 (Vecteur de R^3)

        % Output:
        % M_i1 : Masse totale du lanceur au décollage

        global N_call_f_test;

        M_i1 = (mu + X(1) * (1 + k(1)) + X(2) * (1 + k(2)) + X(3) * (1 + k(3)) );
        
        N_call_f_test = N_call_f_test + 1;

end

