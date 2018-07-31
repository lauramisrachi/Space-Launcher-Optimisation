function [X] =  domain_proj (D, X) 
                
                % Fonction de projection dans l'hypercube D
                % X : vecteur à projeter
                % D : matrice des contraintes min/max pour chaque coordonnées de X 
                %     de taille  [taille vecteur X ; 2 ]
                    
                    
                  for (i = 1:length(X)) 
                            if(X(i) < D(i,1))
                                    X(i) = D(i,1);
                            else if (X(i) > D(i,2))
                                     X(i) = D(i,2);
                            end
                      
                  end
                 

end