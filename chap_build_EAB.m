function [E,Afun,B] = chap_build_EAB(P,D)
% Build E, A(Vx,ax), B from chapter equations (pp.255–257)

% --- E matrix (symmetric descriptor, Eq.7.20 page 255) ---
E = eye(10);

% Fill mechanical inertia block (rows/cols 1–5) per page 255:
% TODO: copy each term (m, mb, mh, etc.) exactly; here schematic:
E(1,1) = P.m;                               % m
E(1,2) = P.m*P.b;                           % m.b
E(1,3) = P.m*P.h;                           % m.h
E(1,4) = P.m*P.e_f;                         % m.ξ_f
E(1,5) = -P.m*D.zb;                         % -m.z_b

E(2,1) = E(1,2);
E(2,2) = P.m*P.b^2 + P.Izz;                 % m\.b^2 + I_zz
E(2,3) = P.m*P.b*P.h - P.Ixz;               % m.b h - I_xz
E(2,4) = P.m*P.e_f*P.b + P.Ifzz*cos(P.eps); % m.ξ_f.b + Ifzz.cos(eps)
E(2,5) = -P.m*P.b*D.zb - P.Ixx*sin(P.eps);  % -m.b.z_b - Ixx.sin(eps)

E(3,1) = E(1,3);
E(3,2) = E(2,3);
E(3,3) = P.m*P.h^2 + P.Ixx;                 % m.h^2 + I_xx
E(3,4) = P.m*P.e_f*P.h + P.Ifzz*sin(P.eps); % m.ξ_f.h + Ifzz.sin(eps)
E(3,5) = -P.m*P.h*D.zb + P.Ixx*cos(P.eps);  % -m.h.z_b + Ixx.cos(eps)

E(4,1) = E(1,4);
E(4,2) = E(2,4);
E(4,3) = E(3,4);
E(4,4) = P.m*P.e_f^2 + P.Ifzz;              % m.ξ_f^2 + I_fzz
E(4,5) = -P.m*P.e_f*D.zb;                   % -m.ξ_f.z_b

E(5,1) = E(1,5);
E(5,2) = E(2,5);
E(5,3) = E(3,5);
E(5,4) = E(4,5);
E(5,5) = P.m*D.zb^2 + P.Izz;                % m z_b^2 + I_zz

% ---- Tyre “rate” states on the diagonal (page shows k_αr, k_αf) ----
E(6,6) = P.kar;                             % k_αr
E(7,7) = P.kaf;                             % k_αf

% ---- Rate-definition rows for φ, δ, β  (unity on diagonal) ----
E(8,8) = 1;
E(9,9) = 1;
E(10,10) = 1;

% --- A matrix as function of Vx ---
Afun = @(Vx,ax) buildA(P,D,Vx,ax);

% --- Input matrix (page 257) ---
B = [0;0;0;1;0;0;0];
end

function A = buildA(P,D,Vx,ax)
A = zeros(10);

% Kinematic wheel spin speeds *** TODO: source? ***
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
A(2,6)  = P.kar * Nr;                                    % a_{2,6}
% a_{2,7}, a_{2,8}, a_{2,9}, a_{2,10}:  % TODO add if present on your page
A(2,7)  = P.kaf*omega*Nf + P.kaf*Nf;    % TODO: *** omega? ***

A(2,8)  = ktr*Nr + (ktf+omega*kgf)*Nf + h*m*ax + ha*Fad + Iwr*omega_r_dot ...
    - Xf*rho_f - Xr *rho_r
A(2,9)  = ktf*Nf*s) + (omega*c - rho_f*s + an)*Xf + mf*ef*ax + Iwf*omega_f_dot*s
A(2,10) = (lb - rho_f*c - omega*s)*Xf + Iwf*omega_f_dot*c - mb*zb*ax + ktf*Nf*c

% ---------------- Row 3: ω_φ (roll rate) ----------------
% Page shows gravity/tire-camber couplings into roll; fill once legible:
% A(3,*) = ... % TODO: transcribe from page 256–257
A(3,2)  = -m*h*Vx - Iwr*omega_r - Iwf*omega_f
A(3,4)  = -Iwf*omega_f*c
A(3,5)  = Iwf*omega_f*s
A(3,8)  = m*g*h - rho_f*Nf - rho_r*Nr
A(3,9)  = (an - rho_f*s)*Nf + mf*ef*g - Iwf*omega_f_dot*c
A(3,10) = (lb - rho_f*c)*Nf - mb*zb*g + Iwf*omega_f_dot*s

% ---------------- Row 4: ω_δ (steer rate) ----------------
A(4,2)  = -mf*ef*Vx - Iwf*omega_f*s
A(4,3)  = Iwf*omega_f*c
A(4,4)  = -P.c_delta;                                    % a_{4,4}
% Other steer couplings (trail/head-angle) usually: a_{4,7}, a_{4,8}, a_{4,9}, a_{4,10}
% A(4,7) = ... ; A(4,8) = ... ; A(4,9) = ... ; A(4,10) = ...  % TODO
A(4,5)  = Iwf*omega_f
A(4,7)  = (kaf*c - an*kaf)*Nf
A(4,8)  = (an*(1 - Kgf)*ktf*c - rho_f*s)*Nf - rho_f*Xf*c + mf*ef*g
A(4,9)  = kgf*an*Nf*s mf*ef*ax*c + A(4,8)*s
A(4,10) = (ktf*c^2 - rho_f*s*c + lb*s)*Nf - mb*zb*(g*s + ax*c) + ...
    (an*s - rho*s*c + lb*s)*Xf + Iwf*omega_f_dot*s

% ---------------- Row 5: ω_β (fork bending rate) ----------------
% Gyroscopic & geometric couplings seen on page (e.g., a_{5,2}, a_{5,3}, a_{5,4}, a_{5,7}, a_{5,8})
% A(5,2) =  P.m * P.b * D.zb * Vx  -  P.Iwf * omega_f * c;     % example seen
% A(5,3) = ... ;  A(5,4) = ... ;  A(5,7) = ... ;  A(5,8) = ...  % TODO
A(5,2)  = mb*zb*Vx - Iwf*omega_f*c
A(5,3)  = -Iwf*omega_f*s
A(5,4)  = -Iwf*omega_f
A(5,7)  = -(kaf*s + lb*kaf)*Nf
A(5,8)  = ((1 - kgf)*lb - ktf*s - rho_f*c)*Nf + rho_f*Xf*s - mb*zb*g
A(5,9)  = kgf*lb*Nf*s - mb*zb*ax*c + A(5,8)*s
A(5,10) = kgf*lb*Nf*c + mb*zb*s + A(5,8)*c - kb

% ---------------- Row 6: α_r (rear slip angle) ----------------
% Entries involving Vx, geometry, and possibly ω_ψ, V_y:
% A(6,1) = ... ; A(6,2) = ... ; etc.   % TODO (page gives explicit forms)
% Often: A(6,1) ~ - (kα_r/N_r) * (...) etc.  % leave TODO to avoid guessing
A(6,1)  = -kgr/Nr
A(6,3)  = 1 - kgr
A(6,6)  = -Vx*kgr/Nr

% ---------------- Row 7: α_f (front slip angle) ----------------
% Page shows dependences on V_y, ω_ψ, δ, β, ε:
% A(7,1) = ... ; A(7,2) = ... ; A(7,4) = ... ; A(7,5) = ... ; A(7,9) = ... ; A(7,10) = ...  % TODO
A(7,1)  = -kyf/Nf
A(7,2)  = -omega*kyf/Nf
A(7,3)  = 1 - kgf
A(7,4)  = (1 - kgf)*s + an*kyf/Nf
A(7,5)  = (1 - kgf)*c + lb*kyf/Nf
A(7,7)  = -Vx*kyf/Nf
A(7,9)  = Vx*c*kyf/Nf
A(7,10) = -Vx*s*kyf/Nf

% ---------------- Rows 8–10: rate definitions ----------------
A(8,3) = 1;                                              % φdot = ω_φ
A(9,4) = 1;                                              % δdot = ω_δ
A(10,5)= 1;                                              % βdot = ω_β

end
