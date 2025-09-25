function D = chap_dependent_params(P)
% Compute dependent geometry and static loads

% *** TODO: 
% find x_f & x_b as a function of other bike parameters
% previous chapter ?
% ***

% Front/rear load distribution *** TODO: check 7.12 extra terms! ***
D.Nf = P.m * P.g * (P.w - P.b) / P.w;
D.Nr = P.m * P.g * P.b / P.w;

% Geometric helpers (Eq.7.21) *** TODO: check e_f (x_f?), e_b (x_b?) ***
D.bf = P.w + (P.e_f + P.a_n - P.h_f * sin(P.eps)) / cos(P.eps);
D.bb = P.w + (P.e_b + P.a_n - P.h_b * sin(P.eps)) / cos(P.eps);
D.zb = P.l_beta + ((P.a_n + P.e_b) * sin(P.eps) - P.h_b) / cos(P.eps);
end
