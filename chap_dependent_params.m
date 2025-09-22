function D = chap_dependent_params(P)
% Compute dependent geometry (Eq.7.21) and static loads

% Front/rear load distribution
D.Nf = P.m * P.g * (P.w - P.b) / P.w;
D.Nr = P.m * P.g * P.b / P.w;

% Geometric helpers (Eq.7.21)
D.bf = (P.w + P.e_f + P.a_n - P.h_f * sin(P.eps)) / cos(P.eps);
D.bb = (P.w + P.e_b + P.a_n - P.h_b * sin(P.eps)) / cos(P.eps);
D.zb = (P.l_beta + (P.a_n + P.e_b) * sin(P.eps) - P.h_b) / cos(P.eps);
end
