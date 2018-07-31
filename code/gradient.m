function [gradf, Jac_c, b, fX, cX] = gradient(f, c, n, m, X, h)

    % Fonction qui calcule le gradient de f au point X

    % Input  
    % X : vecteur de R^n
    % n : nombre de variables pour f et c 
    % m : nombres de contraintes du probl�me
    % f : fonction de R^n dans R 
    % c : fonction de R^n dans R^m repr�sentant les contraintes du probl�me
    % h : précision pour les différences finies (vecteur de taille n)

    % Output
    % gradf : vecteur de R^n du gradient de f en X
    % Jac_c : matrice de dimension m*n qui correspond � la jacobienne de c en X
    % b : -c(X);
    % fX : Valeur de f en x
    % cX : Valeur de c en x


    % Initialisation
    gradf = zeros(n,1) ;
    Jac_c = zeros(m,n) ; 
    fX = f(X);
    cX = c(X) ; 
    b = -cX;
    % Boucle

    for i = 1:n

        Xi = X + h(i) * ones_i(n,i) ; 

        gradf(i) = (f(Xi) - fX)/(h(i)) ; 

        Yi = c(Xi) ; 

        for j = 1:m

            Jac_c(j,i) = (Yi(j) - cX(j)) / (h(i)) ;    


        end

    end
    Jac_c = Jac_c';

    end

    function vec = ones_i(n,i)

    % Input 
    % n : la dimension du vecteur
    % i :  la position o� l'on veut placer un 1

    % Output
    % le vecteur e_i avec un 1 en position i et des 0 ailleurs. 

    if n < i

        disp('Attention il y a un probl�me d appel � la fonction ones_i');

    end
    vec = zeros(n,1) ; 
    vec(i) = 1 ; 

end