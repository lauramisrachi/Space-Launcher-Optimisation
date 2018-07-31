function [ M_i1 ] = f_test2( X )

    % Fonction objectif pour le problème d'étagement d'Ariane 1
    % Input:
    % X vecteur de R^3, contenant les masses d'ergos me1, me2, me3

    % Output:
    % M_i1 : valeur réelle donnant la masse du lanceur au décollage

    global N_call_f_test;

    % Paramètres 
    k1 = 0.1101 ; 
    k2 = 0.1532 ; 
    k3 = 0.2154 ; 
    mu = 1700 ; 


    M_i1 = (mu + X(1) * (1 + k1) + X(2) * (1 + k2) + X(3) * (1 + k3)) ; 

    N_call_f_test = N_call_f_test + 1;

end

