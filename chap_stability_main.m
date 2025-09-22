clear; clc;
P = chap_params_sporttourer();
D = chap_dependent_params(P);

[E,Afun,B] = chap_build_EAB(P,D);

Vx_vec = linspace(1,60,120); % m/s
R = chap_eigen_scan(P,D,Vx_vec,E,Afun);

chap_plots(R);
