function [R, V, M, Rhist, Vhist, Mhist, Time] = simulateur(R0, Rt,  k, me, mu,  theta) 
        
        % Simulateur de trajectoire
        % Input
        % R0 : Position initiale du lanceur (Vecteur de R^2)
        % Rt : Rayon de la Terre
        % k : Constantes constructeurs k1, k2, k3 (Vecteur de R^3)
        % me : vecteur de R^3, contenant les masses d'ergos me1, me2, me3 
        % mu : Masse du satellite
        % theta : Angles de poussée (Vecteur de R^4)
        % Output :
        % R : Matrice 4, 2 des positions intermédaires 
        % V : Matrice 4, 2 des vitesses intermédaires
        % M : Matrice 4, 2 des masses intermédaires
        % Rhist : Matrice :,2 contenant toutes les positions du lanceur au
        %         cours du vol
        % Vhist : Matrice :,2 contenant toutes les vitesses du lanceur au
        %         cours du vol
        % Mhist : Vecteur :,1 contenant l'historique de la masse du lanceur au
        %         cours du vol
        % Time : Vecteur :,1 contenant les pas de temps de la simulation
 

        alpha = [15;10;10];
        ve = [2600; 3000; 4400];
        t = zeros(4,1);
        t(1) = 0;
        
        R(1, :) = R0;
        V(1,:) = 100 * [cos(theta(1)); sin(theta(1))];
        X = me;
        M(1, :) = mu + X(1) * (1 + k(1)) + X(2) * (1 + k(2)) + X(3) * (1 + k(3)) ; 
        
        rho0 = 1.225;
        H = 7000;
        cx= 0.1;
    
        Rhist = R;
        Vhist = V;
        Mhist = M;
        Time(1) = 0.0;
        
        for(j=1:3)
                
                T = alpha(j) * M(j, :);
                q = T / ve(j);
                dt = me(j) / q;
                t(j+1) = t(j) + dt;
                tspan = [t(j) (floor(t(j))+1):1:floor(t(j+1)) t(j+1) ];
                
                func = @(t, y)( fode (t, y, T, q, cx, rho0, Rt, H,  theta(j+1)) );
                
                Zprec(1:2) = R(j,:);
                Zprec(3:4) = V(j,:);
                Zprec(5) = M(j, :);
                
                opts = odeset ('RelTol', 1E-4, 'AbsTol', 1) ;
                [Temp, Z] =   ode45(func, tspan, Zprec, opts);
                n = size(Z,1);
                R(j+1,:) = Z(n, 1:2);
                V(j+1,:)= Z(n, 3:4);
                M(j+1, :) = Z(n, 5) - k(j)*me(j);  
                
                Rhist = [Rhist ; Z(:, 1:2)]; 
                Vhist = [Vhist ; Z(:, 3:4)]; 
                Mhist = [Mhist ; Z(:, 5) ];
                Time = [Time ; Temp] ; 
        end
        
end

function [dy] = fode (t, y, T, q, cx, rho0, Rt, H,  theta)
        
        mu_g = 3.986E14; 
        
        R= [y(1);y(2)];
        V=[y(3); y(4)];
        m = [y(5)];
      
        rho = rho0 * exp(- (norm(R) - Rt) / H);
         gamma = asin (dot(R/norm(R), V/norm(V) )) ;
                
        
        D = -cx * rho * V * norm(V);
        W = -mu_g * (m / norm(R) ^3) * R;
        
          u(1) = ((cos(gamma+theta) + sin(gamma+theta)));
          u(2) = ((-cos(gamma+theta) +  sin(gamma+theta)));

        Tvec = T * u';
        
        dy(1) = y(3);
        dy(2) = y(4);
        dy(3:4) = (1/m) * (Tvec  + W + D);
        dy(5) = -q;
        dy = dy';
end






