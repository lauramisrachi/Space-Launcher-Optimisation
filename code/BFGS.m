function [H] = BFGS (x, gradf, xprec, gradf_prec, Hprec, init)
        
        % Calcule une itération de l'algorithme BFGS. Retourne H, une  approximation de la hessienne de f au point x
        %
        %
        % Input 
        % x : Point de calcul de l'approximation de la hessienne
        % gradf : Gradient de f au point x
        % xprec : Point précédent 
        % gradf_prec : Gradient de f au point précédent
        % Hprec : Approximation au point xprec de la hessienne de f
        % init : Booléen : vaut 1 si BFGS doit renvoyer une matrice
        %        identité
        
        if(init == 1) 
                n = length(x);
                H=eye(n,n);
                return;
        end
        
        d = x - xprec;
        y = gradf - gradf_prec;
        
      
        if (y'*d <= 0) 
            
            H = Hprec ;
            
        else
            
            H = Hprec + ( y*y' / (y' * d)  ) - (  (Hprec*d*d'*Hprec) / (d'*Hprec*d)   );
            
        end
    

end 