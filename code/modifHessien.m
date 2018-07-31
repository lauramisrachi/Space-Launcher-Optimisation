function [Res] = modifHessien(Q) 

        % fonction de modification du Hessien, pour le rendre symétrique et inversible
        % Input : 
        % Q : Hessienne à modifier
        % Output :
        % Res : Hessien modifié
        
        Q = 0.5*(Q+Q');
       
        [n, m] = size(Q);
        v = eig(Q);
         
        if (min(v) < 0) 
                l = -min(v) + 0.1;
                Res = Q + l*eye(n,m);
        else
                Res = Q;
        end
         

end