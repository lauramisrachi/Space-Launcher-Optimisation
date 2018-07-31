function [ Y, in1, in2 ] = g_newton( X, k, ve, Vc )

    % Fonction que l'on veut annuler pour r�soudre th�oriquement le pb
    % d'�tagement 

    w = k ./ (1 + k) ; 
    in1 = (1 / w(1))*(1 - ve(3)/ve(1) * (1 - w(3) * X)) ; 
    in2 = (1 / w(2))*(1 - ve(3)/ve(2) * (1 - w(3) * X)) ;

    Y = ve(1) * log(in1) + ve(2) * log(in2) + ve(3) * log(X) - Vc ; 



end

