% Validation du code du simulateur

clear;
clc;

global N_call_f_test;
global N_call_c_test;
global N_call_armijo;

N_call_f_test = 0;
N_call_c_test = 0;
N_call_armijo = 0;

Rt = 6371E3;
Rc = Rt + 250E3;;

R0 = [Rt;0];

theta = [5;5;0;0];
theta = rad_convert(theta);

k(1) = 0.10; 
k(2) = 0.15 ; 
k(3) = 0.20 ;

mu = 1500 ;

me = [145349; 50215; 7933];


[R, V, M, Rhist, Vhist, Mhist, Time] = simulateur(R0, Rt,  k, me, mu,  theta);
plot_trajectoire(Time, Rhist, Mhist, Vhist, Rt, Rc);
