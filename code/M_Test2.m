% Validation du code d'optimisation avec le problème d'étagement pour
% Ariane 1

clear;
clc;
global N_call_f_test;
global N_call_c_test;
global N_call_armijo;

N_call_f_test = 0;
N_call_c_test = 0;
N_call_armijo = 0;


n=3; % Nombre de variables du problème
m=1; % Nombre de contraintes du problème

X= [140; 310; 7933]; % X de départ

h = 1E-2 * ones(n,1); % Précision gradient
Kmax = 200; % Nombre max d'itérations
eps = 1E-7; % Précision critère d'arrêt
c1 = 0.1; % paramètre Armijo
L = [0]; % Lambda0

Xstar = [145349; 31215; 7933]; % Valeur exacte 
X = [100000;50000;10000];

 fprintf('Validation du code d''optimisation avec le problème d''étagement : \n \n');
 
 [X_SQP, KSQP, DL_hist] = SQP(@f_test2, @c_test2,  [ [10E4, 10E5] ; [10E3, 10E5] ; [10E2, 10E4] ], X, L, Kmax, h, eps, 'BFGS', 3, 1);
 
 fprintf('Résultat obtenu : \n ');
 disp(X_SQP');
 fprintf('Résultat exact : \n ');
 disp(Xstar');
 fprintf('Norme de l\'' erreur : %f \n', norm(X_SQP - Xstar) );
 fprintf('Nombre d\'' itérations : %i \n', KSQP );
 fprintf('Nombre d\'' appels à la fonction f : %i \n', N_call_f_test );
 fprintf('Nombre d\'' appels à la fonction c : %i \n', N_call_c_test );
 fprintf('Nombre de boucles dans Armijo : %i \n', N_call_armijo );

 
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