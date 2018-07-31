function [alpha, success] = armijo(f,c, domain, f_Xk, c_Xk, gradf_xk, Xk, d, c1, rho, lambda, kmax) 
        
        % Algorithme d'Armijo spécialisé pour l'algorithme SQP min f(x) s.c. c(x) = 0
        % Input :
        % f : Fonction objectif
        % c : Fonction contrainte
        % f_Xk : Valeur de la fonction objectif en Xk
        % c_Xk : Valeur de la fonction contrainte en Xk
        % gradf_xk : Gradient de f au point Xk 
        % Xk : Point Xk actuel
        % d : Direction de recherche
        % c1 : Paramètre c1 pour Armijo
        % rho : Paramètre de pénalisation 
        % lambda : Multiplicateur de Lagrange du problème
        % kmax : Nombre d'itération maximal pour Armijo
        % Output :
        % alpha : Pas optimal dans la direction d
        
        global N_call_armijo; 
        
        alpha = 1;
        success = 1;
        k = 0 ; 
        
        Fd = fonction_merite_opt(f_Xk, c_Xk, rho, lambda);
        Lin = c1* (gradf_xk' * d - rho * norm(c_Xk));
        
        Y = abs( Fd * 1E5 );
        
        while ( (Y > Fd  + alpha * Lin)  && k < kmax )             
            k = k+1;
            alpha = alpha/2 ;
            X = Xk + alpha * d  ; 
            %reprojection
            X = domain_proj (domain, X);
            Y = fonction_merite(X, f, c, rho, lambda);
            N_call_armijo = N_call_armijo + 1;
           
        end
        
        if (k == kmax) 
            success=0;
        end
end