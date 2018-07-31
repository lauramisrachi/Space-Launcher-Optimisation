function [Y] = f_test(X)
    
     % Fonction objectif du problème test MHW4D
     % X : Vecteur d'entrée de taille 5
     
     global N_call_f_test;
     
     Y = (X(1) - 1)^2 + (X(1) - X(2)) ^2 + (X(2) - X(3))^3 + (X(3)-X(4))^4 + (X(4) - X(5))^4;
     
     N_call_f_test = N_call_f_test + 1;
end

