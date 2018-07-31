function [X, K, DL_hist] = SQP( f, c, domain, X0, lambda0, Kmax, h, eps,  algo, debug, dist_reproj)

        % Algorithme d'optimisation (SQP) du problème min f(x) sous contrainte c(x)=0
        % 
        % Input 
        % f : Fonction objectif du problème
        % c : Fonction contrainte
        % domain : matrice des contraintes min/max pour chaque coordonnées de X 
        %            de taille  [taille vecteur X ; 2 ]
        % X0 : Point de départ de l'algorithme
        % lambda0 : Point de départ de l'algorithme
        % Kmax : nombre d'itération maximal
        % h : précision pour l'approximation du gradient de f et c
        % eps : précision pour le critère d'arrêt
        % algo : Choix de l'algorithme pour l'approximation de la Hessienne
        %       Valeurs possible : algo = "BFGS" ou algo = "SR1"
        % debug : variable pour le debuggage du programme (mettre à 0)
        % dist_reproj : Booléen : si égal à 1, on projette la direction de
        %               Newton d
        % Output :
        % X : minimiseur du problème
        % K : nombre de pas utilisé
   
        k = 1 ; 
        n = length(X0) ; 
        m = length(c(X0)) ;  % APPEL INUTILE !!!!
        c1 = 0.1;
        rho= 1.0;
        DL_hist = [] ; 

        if (strcmp(algo, 'BFGS') == 0 && strcmp(algo,'SR1') == 0 )

                fprintf('Algorithme %s inconnu \n\n', algo);
                   X=X0;
                return; 
        end

        reinit = 0;  % est-ce qu'il s'agit d'une boucle de réinitialisation du hessien ?
        init=1; % est-ce la première boucle ?
        rhoMax = 1E10;

        X=X0;
        Xprec=X0;
        [gprec, Aprec] = gradient(f, c, n, m, Xprec, h);
        A=Aprec;
        Q= eye(n,n);
        lambda = lambda0;
        DLprec = Q;
        
       %  while (k < Kmax && ( norm(X - Xprec) > eps  || init==1  ) && rho < rhoMax)
        while (k < Kmax && ( norm(DLprec) > eps  || init==1  ) && rho < rhoMax)
            
           % fprintf("SQP : itération %i \n", k);
            Aprec = A;
            [g, A, b, f_xk, c_xk] = gradient(f, c, n, m, X, h);
           
            DL = g +  A * lambda;
            DLprec = gprec + Aprec * lambda;
            
                if (strcmp(algo, 'BFGS'))
                        Q = BFGS (X, DL , Xprec, DLprec, Q, init);
                else
                        Q = SR1 (X, DL , Xprec, DLprec, Q, init);
                end
          
            
            init = 0;
            
            Q = modifHessien(Q);
            
            invQ = Q\eye(n,n);
            lambda = linsolve((-1)*A' * invQ * A, b + A'*invQ*g) ;

            d = linsolve((-1)*Q, g + A*lambda) ;
            
            Xprec=X;
            gprec = g;
            X = X + d; 
           
            % reprojection domaine
            X = domain_proj (domain, X);
            
            % reprojection distance dans certains cas
            if dist_reproj == 1
                d = X-Xprec; % !!!!
            end
            

            if (fonction_merite(X, f, c, rho, lambda) > fonction_merite_opt(f_xk, c_xk, rho, lambda))

                   if ( (gprec' * d - rho * norm(c_xk) ) < 0)
                       
                       [alpha, success] = armijo(f,c, domain,f_xk, c_xk,  g, Xprec, d, c1, rho, lambda, 10);
                      
                       if(success == 1)
                           X = Xprec + alpha*d;
                           X = domain_proj (domain, X);
                           reinit = 0;
                       else  
                           if(reinit == 0)
                               X = Xprec;
                               Q =eye(size(Q));
                           else 
                               rho = rho * 10.0;
                           end
                               
                           reinit=1;
                       end
                   else    
                        if (reinit == 1)
                            rho = rho * 10.0;
                        else
                            Q =eye(size(Q));
                            reinit = 1;
                        end
                   end
                   
            else
                reinit = 0;
                   
            end

            k = k+1;
            
            % R�cup�rer les normes du gradient /x du Lagrangien
            if nargout > 2 
                DL_hist = [DL_hist; norm(DL)] ;     
            end
            
            
            % Options de d�buggage:
            
            % Probl�me angles
            if(debug == 1)
                angles = X * 180 / pi;
                fprintf('angles : %f %f %f %f, hauteur : %f, vitesse : %f \n', angles(1), angles(2), angles(3), angles(4), c(X), f(X));
            
            % Probl�me MWH4D    
            elseif debug == 2 
                fprintf('iter = %d,   x = [%f, %f, %f, %f, %f],   lambda = [%f,%f,%f],     f = %f,    c = [%f, %f, %f],    rho =  %f,     norm(grad(L)) = %f  \n\n', k-1, X, lambda, f(X), c(X), rho, norm(DL)) ;
               
                
            % Probl�me �tagement Ariane 1 
            elseif debug == 3
                fprintf('iter = %d,   x = [%f, %f, %f],   lambda = [%f],     f = %f,    c = [%f],    rho =  %f,     norm(grad(L)) = %f  \n\n', k-1, X, lambda, f(X), c(X), rho, norm(DL)) ; 
                
            end

        end
        
        
        K = k;
        
        % Probl�me �tagement Ariane 6 (juste les derni�res it�rations)
        if (debug == 4)
            fprintf('Pour les masses d ergols: kiter = %d, x = [%f, %f, %f], lambda = [%f], f = %f,    c = [%f],    rho =  %f,     norm(grad(L)) = %f  \n\n', K, X, lambda, f(X), c(X), rho, norm(DL));
        end 
        
        % Probl�me angle Ariane 6 (juste les derni�res it�rations)
        if (debug == 5)
            fprintf('Pour les angles : kiter = %d, x = [%f, %f, %f, %f], lambda = [%f, %f], f = %f,    c = [%f, %f],    rho =  %f,     norm(grad(L)) = %f  \n\n', K, X, lambda, f(X), c(X), rho, norm(DL));
        end

end

