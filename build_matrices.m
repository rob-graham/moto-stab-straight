function mats = build_matrices(P)
%BUILD_MATRICES  Construct Whipple/Sharp small-angle model matrices.
%   MATS = BUILD_MATRICES(P) constructs the mass, damping, and stiffness
%   matrices (M, C1, K0, K2) of the classic Whipple/Sharp motorcycle model
%   using the parameter struct P returned by PARAMS_SPORT_TOURER. Steering
%   damper effects are added as a viscous contribution. All matrices are
%   2-by-2 and correspond to the lean and steer coordinates ordered as
%   q = [phi; delta].
%
%   The returned struct contains:
%     M  - mass matrix multiplying qddot
%     C1 - speed-proportional damping matrix
%     K0 - gravity-proportional stiffness matrix
%     K2 - speed-squared stiffness (geometric) matrix
%     D  - additional viscous damping (steering damper) matrix
%
%   Example:
%     P = params_sport_tourer();
%     mats = build_matrices(P);
%     disp(mats.M);
%
%   The implementation follows the notation in Meijaard et al. (2007) and
%   Sharp (2001), adapted for motorcycle parameters.

arguments
    P struct
end

lambda = P.steer_axis_tilt;
c = P.trail;
w = P.wheelbase;
rF = P.front_radius;
rR = P.rear_radius;

% Masses and centres of mass
mR = P.m_frame;
mF = P.m_fork + P.m_front_wheel;

xR = P.x_frame;
zR = P.z_frame;
xF = P.x_fork;
zF = P.z_fork;

mT = mR + mF;
xT = (mR * xR + mF * xF) / mT;
zT = (mR * zR + mF * zF) / mT;

% Inertias about component centres (augment fork with wheel)
IFr = P.I_frame;          % frame (rear)
IFo = P.I_fork;           % fork/handlebar assembly
Iw = P.I_front_wheel;     % wheel inertias (yy spin, xx radial)

IFork.xx = IFo.xx + Iw.xx;
IFork.yy = IFo.yy + Iw.yy;
IFork.zz = IFo.zz + Iw.xx;  % assume radial inertia about vertical
IFork.xz = IFo.xz;

% Total body inertias about the total CoM (parallel axis)
dxR = xR - xT;
dzR = zR - zT;
dxF = xF - xT;
dzF = zF - zT;

ITxx = IFr.xx + mR * dzR^2 + IFork.xx + mF * dzF^2;
ITzz = IFr.zz + mR * dxR^2 + IFork.zz + mF * dxF^2;
ITxz = IFr.xz + mR * dxR * dzR + IFork.xz + mF * dxF * dzF;

% Convenience combinations for the front assembly
IFxx = IFork.xx + mF * (zF - rF)^2; % include elevation of wheel contact
IFyy = IFork.yy;
IFzz = IFork.zz + mF * (xF - w)^2;
IFxz = IFork.xz + mF * (zF - rF) * (xF - w);

% Mass matrix
M11 = mT * zT^2 + ITxx;
M12 = mF * zF * zT + IFxx + mF * c * zF * cos(lambda);
M22 = IFxx + mF * c^2 + IFyy * sin(lambda)^2 + IFzz * cos(lambda)^2 + 2 * IFxz * sin(lambda) * cos(lambda);

M = [M11, M12; M12, M22];

% Damping matrix proportional to speed
C1 = zeros(2);
C1(1,2) = (mT * zT + ITxz) * sin(lambda) / w;
C1(2,1) = -(mT * zT + ITxz) * sin(lambda) / w; % gyroscopic coupling
C1(2,2) = (IFxz + IFyy * tan(lambda) + mF * c * zF * cos(lambda)) * sin(lambda) / w;

% Gravity-dependent stiffness
K0 = zeros(2);
K0(1,1) = -mT * zT;
K0(1,2) = -mF * (zF + c * cos(lambda));
K0(2,1) = K0(1,2);
K0(2,2) = -(IFxz * sin(lambda) + IFzz * cos(lambda)) * cos(lambda) / w - mF * c * (zF * sin(lambda) + c * cos(lambda)) / w;

% Speed-squared stiffness
K2 = zeros(2);
K2(1,2) = (mT * zT + ITxz) * sin(lambda) / w;
K2(2,1) = -K2(1,2);
K2(2,2) = (IFyy * sin(lambda) + IFxz * cos(lambda)) * sin(lambda) / w + mF * c * zF * sin(lambda) / w;

% Steering damper viscous contribution (acts on steer rate)
D = zeros(2);
D(2,2) = P.c_steer;

mats = struct('M', M, 'C1', C1, 'K0', K0, 'K2', K2, 'D', D);
end
