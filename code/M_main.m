% Programme principal de résolution du problème Ariane 6

clear;
clc;

global N_call_f_test;
global N_call_c_test;
global N_call_armijo;

format long;

%w = warning('off', 'all');

Rt = 6378E3;
Rc = Rt + 250E3;
mu_g = 3.986E14;
Vc = sqrt(mu_g / Rc);

theta = [5;5;0;0];
theta = rad_convert(theta);
R0 = [Rt;0];

k(1) = 0.10; 
k(2) = 0.15 ; 
k(3) = 0.20 ;
k=k';

Ve(1) = 2600;
Ve(2) = 3000;
Ve(3) = 4400;
Ve = Ve';

mu = 1500; %2500

%me = [145349; 31215; 7933];
me = [10000; 5000; 3000];


h_etagement = 1;
tol_etagement = 1E-3;

h_theta = pi / 180; %1E-2
tol_theta = 1E-3;

eps= 1;
kmax = 200;%20
kmax_loc = 300; %50

kiter=1;
dv = 0.2*Vc;
Vp = Vc;

Vp = Vp + dv;

dom_angles = [-50, 50] * pi / 180;

tic;
[R, V, M, Rhist, Vhist, Mhist, Time] = simulateur(R0, Rt,  k, me, mu,  theta);
toc;

dv_hist(1) = dv ;%historique des dv
h_hist(1) = norm(R(4,:)) - Rc; %historique des h

while (kiter < kmax && (norm(dv) > eps) ) 
           
        
        fprintf('********************* ITERATION %i ****************************** \n \n', kiter);
        
        
         N_call_f_test = 0;
         N_call_c_test = 0;
         N_call_armijo = 0;
         
         plot_trajectoire(Time, Rhist, Mhist, Vhist, Rt, Rc);
         
      
         c_et = @(X) ( c_etagement(X, mu, k, Ve, Vp) ) ;
         f_et = @(X) (f_etagement(X, mu, k));
         
         kconfig = 0;
         fprintf('Début SQP pour les masses d ergol\n');
         [me, kconfig] = SQP(f_et, c_et, [ [5000, 100000] ; [3000, 50000] ; [100, 10000] ], me, [0], kmax_loc, h_etagement*ones(3,1), tol_etagement, 'BFGS', 0, 1); % Configuration
         fprintf('fin. (%i étapes) \n', kconfig);
         
        
         simul = @(theta) ( f_theta(R0, Rt, k, me, mu, theta)  );
         contr = @(theta) ( c_theta(R0, Rt, Rc, k,  me, mu, theta) ) ;
            
         ktheta=0;
        
         fprintf('Début SQP pour les angles theta\n');
         [theta, ktheta] = SQP(simul, contr,[ [-45, 45]*pi/180;  [-45, 45]*pi/180 ;  [-90, 90]*pi/180  ; [-180, 180]*pi/180  ]  , theta, [0;0], kmax_loc, h_theta*ones(4,1), tol_theta,  'BFGS', 0, 1);   % THETA
          fprintf('fin. (%i étapes) \n', ktheta);
          
         [R, V, M, Rhist, Vhist, Mhist, Time] = simulateur(R0, Rt,  k, me, mu,  theta);
        
        dv = Vc - norm(V(4,:));
        Vp = Vp + dv*0.1;
       
         
         Cst = contr(theta);
         
         fprintf('Itération %i : \n', kiter);
         fprintf('dv : %f \n', norm(dv));
         fprintf('Vitesse orbitale : %f \n', norm(V(4,:)));
         fprintf('c(1) : %f \n', Cst(1));
         fprintf('c(2) : %f \n', Cst(2));
         fprintf('me : ');
         disp(me);
         fprintf('theta');
         disp(theta);
         fprintf('Nombre d\'' itérations : %i \n', kconfig+ktheta);
         fprintf('Nombre d\'' appels à la fonction f : %i \n', N_call_f_test );
         fprintf('Nombre d\'' appels à la fonction c : %i \n', N_call_c_test );
         fprintf('Nombre de boucles dans Armijo : %i \n\n', N_call_armijo );
        
          kiter = kiter +1;
          
          dv_hist(kiter) = dv ;%historique des dv
          h_hist(kiter) = abs(norm(R(4,:)) - Rc); %historique des h
          
          %%% pour les masses d'ergoles
          %kiter, me, kconfig (nbr it�rations internes) 
end

% Plot de la trajectoire
plot_trajectoire(Time, Rhist, Mhist, Vhist, Rt, Rc);

% Plot de l'erreur en vitesse et en hauteur
figure1 = figure;
subplot1 = subplot(1,2,1,'Parent',figure1);
hold(subplot1,'on');
plot(abs(dv_hist),'Parent',subplot1,'LineWidth',2,'LineStyle','-.');
xlabel({'Nombre d''it�rations k'});
ylabel({'||dv||    (m.s^{-1})'});
box(subplot1,'on');
set(subplot1,'FontSize',12);
subplot2 = subplot(1,2,2,'Parent',figure1);
hold(subplot2,'on');
plot(abs(h_hist),'Parent',subplot2,'LineWidth',2,'LineStyle','-.');
xlabel({'Nombre d''it�rations k'});
ylabel({'||dh||   (m)'});
box(subplot2,'on');
set(subplot2,'FontSize',12);
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
    'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

text(0.5, 1,'\bf Norme de l''erreur en vitesse et en altitude en fonction du nombre d''it�ration' ...
    ,'HorizontalAlignment' ,'center','VerticalAlignment', 'top')

