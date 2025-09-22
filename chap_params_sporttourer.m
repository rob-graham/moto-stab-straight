function P = chap_params_sporttourer()
% Sport touring motorcycle parameters (Table 7.2, p.260)
% All SI units; caster angle converted to rad; kNm→N·m

P.g      = 9.80665;       % gravity [m/s^2]

% Geometry
P.w      = 1.448;         % wheelbase [m]
P.eps    = deg2rad(26.80);% caster angle ε [rad]
P.a_n    = 0.105;         % mechanical trail [m]

% Masses
P.m0     = 195.0;         % vehicle mass without rider [kg]
P.m      = 270.0;         % overall mass with rider [kg]
P.m_f    = 34.0;          % front assembly [kg]
P.m_beta = 18.0;          % fork bending lump [kg]

% CoM positions
P.b0     = 0.722;  P.h0 = 0.482; % vehicle only [m]
P.b      = 0.688;  P.h  = 0.640; % with rider [m]
P.e_f    = 0.025;  P.h_f= 0.600; % front assembly [m]
P.e_b    = 0.000;  P.h_b= 0.350; % bending lump [m]

% Inertias
P.I0xx   = 13.5;   P.I0xz= 3.0;   P.I0zz= 55.0;  % vehicle only
P.Ixx    = 35.5;   P.Ixz= -1.7;   P.Izz = 59.3;  % overall with rider
P.Ifzz   = 0.83;                                   % front assy yaw inertia
P.Iwf    = 0.600;  P.Iwr = 0.800;                  % wheel spin inertia
P.Ibxx   = 0.8;                                    % bending inertia

% Tyres
P.Rf     = 0.294;  P.Rr  = 0.299;  % rolling radii [m]
P.rhof   = 0.064;  P.rhor= 0.078;  % cross-section radii [m]

% Normalised stiffness
P.kaf    = 16.00;  P.kar = 14.50;  % cornering [1/rad]
P.kgf    = 0.850;  P.kgr = 0.950;  % camber [1/rad]
P.ksf    = -0.200; P.ksr = -0.200; % self-aligning [m/rad]
P.ktf    = 0.015;  P.ktr = 0.018;  % twist [m/rad]
P.klf    = 160000; P.klr = 140000; % transverse [N/m]

% Fork bending
P.l_beta = 0.670;          % fork bending axis pos [m]
P.k_beta = 38.0e3;         % bending stiffness [N·m/rad]

% Steering damper
P.c_delta= 1.0;            % [N·m·s/rad]

% Aero
P.CdA    = 0.467;  P.hA = 0.350;
end
