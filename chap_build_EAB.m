function [E,Afun,B] = chap_build_EAB(P,D)
% Build E, A(Vx,ax), B from chapter equations (pp.255–257)

% --- E matrix (symmetric descriptor, Eq.7.20 page 255) ---
E = eye(10);
% Fill mechanical inertia block (rows/cols 1–5) per page 255:
% TODO: copy each term (m, mb, mh, etc.) exactly; here schematic:
E(1,1)=P.m;  E(1,2)=P.m*P.b;  E(1,3)=P.m*P.h;
E(2,2)=P.m*P.b^2 + P.Izz;  E(2,3)=P.m*P.b*P.h - P.Ixz;
E(3,3)=P.m*P.h^2 + P.Ixx;
% etc... (Codex: complete symmetric fill with page 255 entries)

% Tyre states (diag stiffness factors in E as per page 255)
E(6,6)=P.kar;  E(7,7)=P.kaf;

% φ, δ, β rate-def rows
E(8,8)=1; E(9,9)=1; E(10,10)=1;

% ***
E(1,1) = P.m;                         % m
E(1,2) = P.m*P.b;                     % m b
E(1,3) = P.m*P.h;                     % m h
E(1,4) = P.m*P.e_f;                   % m ξ_f
E(1,5) = -P.m*D.zb;                   % -m z_b

E(2,1) = E(1,2);
E(2,2) = P.m*P.b^2 + P.Izz;           % m b^2 + I_zz
E(2,3) = P.m*P.b*P.h - P.Ixz;         % m b h - I_xz
E(2,4) = P.m*P.e_f*P.b;               % m ξ_f b     % TODO: confirm exact term on page
E(2,5) = -P.m*P.b*D.zb;               % -m b z_b    % TODO: confirm exact power/extra terms

E(3,1) = E(1,3);
E(3,2) = E(2,3);
E(3,3) = P.m*P.h^2 + P.Ixx;           % m h^2 + I_xx
E(3,4) = P.m*P.e_f*P.h;               % m ξ_f h     % TODO: confirm
E(3,5) = -P.m*P.h*D.zb;               % -m h z_b    % TODO: confirm

E(4,1) = E(1,4);
E(4,2) = E(2,4);
E(4,3) = E(3,4);
E(4,4) = P.Ifzz + P.m*P.e_f^2;        % I_fzz + m ξ_f^2
E(4,5) = -P.m*P.e_f*D.zb;             % -m ξ_f z_b  % TODO: confirm

E(5,1) = E(1,5);
E(5,2) = E(2,5);
E(5,3) = E(3,5);
E(5,4) = E(4,5);
E(5,5) = P.m*D.zb^2 + P.Izz;          % m z_b^2 + I_zz   % TODO: confirm exact inertia here

% ---- Tyre “rate” states on the diagonal (page shows k_{αr}, k_{αf}) ----
E(6,6) = P.kar;                        % k_{αr}
E(7,7) = P.kaf;                        % k_{αf}

% ---- Rate-definition rows for φ, δ, β  (unity on diagonal, see page’s last three “1”) ----
E(8,8) = 1;
E(9,9) = 1;
E(10,10) = 1;

% ***

% --- A matrix as function of Vx, ax ---
Afun = @(Vx,ax) buildA(P,D,Vx,ax);

% --- Input matrix (page 257) ---
B = [0;0;0;1;0;0;0;0;0;0];
end

function A = buildA(P,D,Vx)
A = zeros(10);

% Kinematic wheel spin speeds
omega_f = Vx / P.Rf;
omega_r = Vx / P.Rr;

% Convenience (seen in page formulas)
Nf = D.Nf;  Nr = D.Nr;
c = cos(P.eps); s = sin(P.eps);

% The page uses X_f as a lever term coupling front lateral force to (δ,β) rows.
% Until we can read the exact X_f definition, keep units-consistent placeholder:
Xf = P.kaf * Nf * P.a_n;   % TODO: replace with exact chapter X_f expression

% ---------------- Row 1: V_y dynamics ----------------
A(1,2)  = -P.m * Vx;                                     % a_{1,2}
A(1,6)  = P.kar * Nr;                                    % a_{1,6}
A(1,7)  = P.kaf * Nf;                                    % a_{1,7}
A(1,8)  = P.kgf * Nf + P.kgr * Nr;                       % a_{1,8}
A(1,9)  =  Xf * c;                                       % a_{1,9}
A(1,10) = -Xf * s;                                       % a_{1,10}

% ---------------- Row 2: ω_ψ (yaw rate) ----------------
A(2,2)  = -P.m * P.b * Vx;                               % a_{2,2}
A(2,3)  = P.Iwr*omega_r + P.Iwf*omega_f;                 % a_{2,3}
A(2,4)  = P.Iwf*omega_f * s;                             % a_{2,4}
A(2,5)  = P.Iwf*omega_f * c;                             % a_{2,5}
A(2,6)  = P.kar * Nr;                                    % a_{2,6}   % TODO: confirm presence/value on page
% a_{2,7}, a_{2,8}, a_{2,9}, a_{2,10}:  % TODO add if present on your page

% ---------------- Row 3: ω_φ (roll rate) ----------------
% Page shows gravity/tire-camber couplings into roll; fill once legible:
% A(3,*) = ... % TODO: transcribe from page 256–257

% ---------------- Row 4: ω_δ (steer rate) ----------------
A(4,4)  = -P.c_delta;                                    % a_{4,4}
% Other steer couplings (trail/head-angle) usually: a_{4,7}, a_{4,8}, a_{4,9}, a_{4,10}
% A(4,7) = ... ; A(4,8) = ... ; A(4,9) = ... ; A(4,10) = ...  % TODO

% ---------------- Row 5: ω_β (fork bending rate) ----------------
% Gyroscopic & geometric couplings seen on page (e.g., a_{5,2}, a_{5,3}, a_{5,4}, a_{5,7}, a_{5,8})
% A(5,2) =  P.m * P.b * D.zb * Vx  -  P.Iwf * omega_f * c;     % example seen
% A(5,3) = ... ;  A(5,4) = ... ;  A(5,7) = ... ;  A(5,8) = ...  % TODO

% ---------------- Row 6: α_r (rear slip angle) ----------------
% Entries involving Vx, geometry, and possibly ω_ψ, V_y:
% A(6,1) = ... ; A(6,2) = ... ; etc.   % TODO (page gives explicit forms)
% Often: A(6,1) ~ - (kα_r/N_r) * (...) etc.  % leave TODO to avoid guessing

% ---------------- Row 7: α_f (front slip angle) ----------------
% Page shows dependences on V_y, ω_ψ, δ, β, ε:
% A(7,1) = ... ; A(7,2) = ... ; A(7,4) = ... ; A(7,5) = ... ; A(7,9) = ... ; A(7,10) = ...  % TODO

% ---------------- Rows 8–10: rate definitions ----------------
A(8,3) = 1;                                              % φdot = ω_φ
A(9,4) = 1;                                              % δdot = ω_δ
A(10,5)= 1;                                              % βdot = ω_β

% ***

% Fill only non-null entries given in pp.256–257
A(1,2)  = -P.m * Vx;
A(1,6)  = P.kar * D.Nr;
A(1,7)  = P.kaf * D.Nf;
A(1,8)  = P.kgf*D.Nf + P.kgr*D.Nr;
% TODO: Define Xf correctly from chapter, placeholder:
Xf = P.a_n * P.kaf * D.Nf;
A(1,9)  = Xf * cos(P.eps);
A(1,10) = -Xf * sin(P.eps);

A(2,2)  = -P.m * P.b * Vx;
A(2,3)  = P.Iwr*omega_r + P.Iwf*omega_f;
A(2,4)  = P.Iwf*omega_f*sin(P.eps);
A(2,5)  = P.Iwf*omega_f*cos(P.eps);
A(2,6)  = P.kar*D.Nr;
% … etc — Codex: continue filling with all a2,j, a3,j, a4,j, a5,j,
% and a6–a7 tyre coupling entries as listed in your photos.

% Steering damping
A(4,4)  = -P.c_delta;

% Rate-def rows (Eq.7.20, page 257)
A(8,3)=1;  A(9,4)=1;  A(10,5)=1;

end
