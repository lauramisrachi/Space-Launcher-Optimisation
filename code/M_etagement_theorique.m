%% Script pour la résolution analytique du problème d'étagement %%

% Cas tests : 
% - Ariane 1
% - Ariane 6

prompt = 'Tapez 1 pour le cas test Ariane 6 ou 2 pour le cas test Ariane 1 \n' ; 
cas_test = input(prompt) ; 

if (cas_test ~= 1) || (cas_test ~= 2)
    disp('Entrez 1 ou 2 pour le cas test')
end

%% Données Ariane 6

if cas_test == 1
    
    Rt = 6378E3;
    Rc = Rt + 250E3;
    mu_g = 3.986E14;
    Vc = sqrt(mu_g / Rc);

    k(1) = 0.10; 
    k(2) = 0.15 ; 
    k(3) = 0.20 ;
    k=k';

    Ve(1) = 2600;
    Ve(2) = 3000;
    Ve(3) = 4400;
    Ve = Ve';
    
    mu = 1500 ; 
    
end

%% Données Ariane 1 

if cas_test == 2
    
    Vc = 11527 ; 

    k(1) = 0.1101; 
    k(2) = 0.1532 ; 
    k(3) = 0.2154 ;
    k=k';

    Ve(1) = 2647.2;
    Ve(2) = 2922.4;
    Ve(3) = 4344.3;
    Ve = Ve';
    
    mu = 1700 ; 
    
end



%% Algorithme de Newton approché (méthode de la sécante)

% On définit une fonction anonyme en X
g_newton_x = @(X) g_newton(X, k, Ve, Vc) ; 

Kmax = 20 ;

% Ce point de départ a été choisi car il assure que le log est bien défini 
X = 2.5 ;
epsilon = 0.00001 ; 
k_iter = 0 ; 
gX1 = 1 ; 
a = 1 ; 

% Boucle 

while (abs(gX1) > epsilon) && (k_iter < Kmax)
         
    [gX, in1, in2] = g_newton_x(X) ;
    
    if (in1 <= 0) || (in2 <= 0)
        X = 2.5 ;   
        [gX, ~, ~] = g_newton_x(X) ;
    end
    
    X1 = X - gX/a ; 
    [gX1, in1, in2] =  g_newton_x(X1) ; 
    
    if (in1 <= 0) || (in2 <= 0)
        X1 = 2.6 ;  
        [gX1, ~, ~] =  g_newton_x(X1) ;
    end
    
    a = (gX - gX1) / (X - X1) ; 
    k_iter = k_iter + 1 ;
    X = X1 ; 
    
    
    
end


fprintf('Le point d annulation de la fonction g est obtenu pour X = %f \n\n',  X) ; 
fprintf('Le nombre d itération de la méthode de la sécante vaut %d \n\n',  k_iter) ; 
fprintf('La valeur de la fonction a annuler au point obtenu vaut %2.8f \n\n',  gX1) ; 

%% Résolution des problèmes Ariane 

[~, x(1), x(2)] = g_newton_x(X) ;
x(3) = X ;

M = ones(4,1); 
M(4) = mu ; 

for j = [3,2,1]
    M(j) = M(j + 1) / ((k(j) + 1)/x(j) - k(j)) ;
end

M = M(1:3) ;
me = ones(3,1) ; 
for j = 1:3
    me(j) = M(j) * (1 - (k(j) + 1)/x(j) + k(j)) / (1 + k(j)) ;
end

if cas_test == 1
    fprintf('Le résultat théorique pour les masses d ergol d Ariane 6 vaut: [%f, %f, %f] \n\n',  me) ; 
end

if cas_test == 2
    fprintf('Le résultat théorique pour les masses d ergol d Ariane 1 vaut: [%f, %f, %f] \n\n',  me) ; 
end


