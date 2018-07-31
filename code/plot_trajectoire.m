function plot_trajectoire(Time, R, M, V,  Rt, Rc ) 
        
        % Fonction qui permet l'affichage de la trajectoire du lanceur
        % Input :
        % Time : Vecteur contenant les pas de temps de la simulation
        % R : Matrice [Taille de Time x 2] contenant l'historique des
        %      positions du lanceur
        % M : Vecteur de même taille que Time contenant l'historique de la
        %       masse lors du vol
        % V : Matrice [Taille de Time x 2] contenant l'historique de
        %     la vitesse du lanceur
        % Rt : Rayon de la Terre
        % Rc : Rayon de l'orbite à atteindre (Rt + hauteur de l'orbite)
        

        figure(1);
        subplot(2,2,1);
        % Terre
        x=[0:0.01:2* pi];
        plot (cos(x)*Rt, sin(x)*Rt);
        
        hold on;
        plot(R(:,1), R(:, 2), '-');
        
        hold on;
        plot (cos(x)*Rc, sin(x)*Rc, '-');
        
        hold off;
        title('Trajectoire espace des phases');
        
        %figure(2) ;
        subplot(2,2,2);
        plot(Time, M);
        title('Masse du lanceur');
        
        %figure(3) ;
        subplot(2,2,3);
        for( i = 1:size(Time))
            normV(i) = norm(V(i, :));
        end
        
        plot(Time, normV);
        title('Vitesse du lanceur');
      
        
        %figure(4);
        subplot(2,2,4);
        
        for( i = 1:size(Time))
            dist(i) = norm(R(i, :)) - Rt;
        end
        
        plot(Time, dist);
        hold on;
        plot([1:4], [Rc-Rt;Rc-Rt;Rc-Rt;Rc-Rt]);
        hold off;
        
        title('Hauteur du lanceur');
        
end