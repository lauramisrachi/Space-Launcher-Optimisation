function Y = c_test(X)
    
     % Fonction représentant les contraintes du problème test MHW4D 
     % Input :
     % X : Vecteur d'entrée de taille 5
    
     global N_call_c_test;
     
     Y(1) = X(1) + X(2)^2 + X(3)^2 - 3*sqrt(2)-2 ; 
     Y(2) = X(2) - X(3)^2 + X(4) - 2*sqrt(2) + 2;
     Y(3) = X(1)*X(5) - 2  ;
     
     N_call_c_test = N_call_c_test + 1;

    Y = Y';
end

