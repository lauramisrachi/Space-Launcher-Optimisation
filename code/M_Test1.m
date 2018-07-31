feature('DefaultCharacterSet', 'UTF8')
% Validation du code d'optimisation avec le problème MHW4D

clear;
clc;
global N_call_f_test;
global N_call_c_test;
global N_call_armijo;

N_call_f_test = 0;
N_call_c_test = 0;
N_call_armijo = 0;


n=5; % Nombre de variables du problème
m=3; % Nombre de contraintes du problème

X= [-1; 2; 1; -2; -2]; % X de départ

h = 1E-5 * ones(n,1); % Précision gradient
Kmax = 20; % Nombre max d'itérations
eps = 1E-15; % Précision critère d'arrêt
c1 = 0.1; % paramètre Armijo
L = [0;0;0]; % Lambda0

Xstar = [-1.2366;2.4616;1.1911;-0.2144;-1.6165]; % Valeur exacte 

 fprintf('Validation du code d''optimisation avec le problème 1 : \n \n');
 
 [X_SQP, KSQP, DL_hist] = SQP(@f_test, @c_test, [ [-5 5] ; [-5 5] ;[-5 5] ;[-5 5] ;[-5 5] ], X, L, Kmax, h, eps, 'BFGS', 0, 0);
 
 % Display
 fprintf('Résultat obtenu : \n ');
 disp(X_SQP');
 fprintf('Résultat exact : \n ');
 disp(Xstar');
 fprintf('Norme de l\'' erreur : %f \n', norm(X_SQP - Xstar) );
 fprintf('Nombre d\'' itérations : %i \n', KSQP );
 fprintf('Nombre d\'' appels à la fonction f : %i \n', N_call_f_test );
 fprintf('Nombre d\'' appels à la fonction c : %i \n', N_call_c_test );
 fprintf('Nombre de boucles dans Armijo : %i \n', N_call_armijo );
 
 
 % Plot du gradient du Lagrangien en fonction du nombre d'iterations
 figure1 = figure;
 axes1 = axes('Parent',figure1);
 hold(axes1,'on');
 plot(1:1:KSQP - 1,DL_hist,'LineWidth',2,'LineStyle','-.');
 xlabel('Num�ro d''it�tation k');
 title('Norme du gradient par rapport � x du Lagragien en fonction du nombre d''it�ration',...
    'FontSize',14,...
    'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
 ylabel('||\nabla_x(L)|| ','EdgeColor',[1 1 1]);

 box(axes1,'on');
 set(axes1,'FontSize',12);