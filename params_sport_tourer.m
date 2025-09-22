function P = params_sport_tourer()
%PARAMS_SPORT_TOURER  Sport-tourer motorcycle parameters for stability study.
%   P = PARAMS_SPORT_TOURER() returns a struct containing the geometric,
%   inertial, and tyre properties of a representative sport-tourer
%   motorcycle. All quantities are given in SI units. The helper function
%   CHECK_UNITS(P) validates the supplied numbers against plausible ranges
%   to catch accidental unit or transcription mistakes.
%
%   The returned struct provides the fields:
%     g            - gravitational acceleration (m/s^2)
%     wheelbase    - wheelbase (m)
%     steer_axis_tilt - steering axis tilt from vertical (rad)
%     trail        - mechanical trail (m)
%     front_radius - effective rolling radius of front tyre (m)
%     rear_radius  - effective rolling radius of rear tyre (m)
%     m_frame      - mass of the main frame, rider, and rear wheel (kg)
%     m_fork       - mass of the front assembly (kg)
%     x_frame      - longitudinal CoM of the frame (m, w.r.t. rear contact)
%     z_frame      - vertical CoM of the frame (m, positive upward)
%     x_fork       - longitudinal CoM of the front assembly (m)
%     z_fork       - vertical CoM of the front assembly (m)
%     I_frame      - struct with inertia tensor components of frame (kg·m^2)
%                    fields: xx, yy, zz, xz (about frame CoM axes)
%     I_fork       - struct with inertia tensor components of front (kg·m^2)
%                    fields: xx, yy, zz, xz
%     m_front_wheel- mass of the front wheel (kg)
%     I_front_wheel- polar and radial inertia of front wheel (kg·m^2) with
%                    fields: yy (spin) and xx (radial)
%     tyre_stiff_lat - lateral tyre stiffness (N/rad)
%     tyre_stiff_align - aligning torque stiffness (N·m/rad)
%     c_steer      - viscous damping coefficient of the steering damper
%                    (N·m·s/rad)
%     speed_range  - vector [vmin vmax] in m/s suggested for analysis
%
%   Example:
%     P = params_sport_tourer();
%     disp(P.wheelbase);
%
%   See also BUILD_MATRICES, EIGEN_SCAN, PLOTS_MODES.

% Nominal sport-tourer numbers gathered from literature and manufacturer
% data sheets. The geometry matches a 25 deg head angle and 100 mm trail
% machine with a 1.48 m wheelbase.
P.g = 9.81;
P.wheelbase = 1.48;
P.steer_axis_tilt = deg2rad(25);      % rad, from vertical
P.trail = 0.10;
P.front_radius = 0.33;
P.rear_radius = 0.32;

% Mass and inertia properties
P.m_frame = 215;                      % kg (frame + rider + rear wheel)
P.m_fork = 35;                        % kg (steering head, fork, bars)
P.x_frame = 0.55;                     % m forward of rear contact
P.z_frame = 0.55;                     % m above ground
P.x_fork = 1.35;                      % m (near front axle)
P.z_fork = 0.70;                      % m above ground

% Frame inertia about its CoM (principal axes tilted: keep xz coupling)
P.I_frame.xx = 28.0;
P.I_frame.yy = 52.0;
P.I_frame.zz = 35.0;
P.I_frame.xz = 3.5;

% Front assembly inertia about its CoM
P.I_fork.xx = 4.5;
P.I_fork.yy = 5.0;
P.I_fork.zz = 2.5;
P.I_fork.xz = 0.90;

% Front wheel properties (polar and radial inertia)
P.m_front_wheel = 16;
P.I_front_wheel.xx = 0.60;            % about axle (radial)
P.I_front_wheel.yy = 1.10;            % spin inertia

% Tyre stiffness (linearised about straight running)
P.tyre_stiff_lat = 160000;            % N/rad
P.tyre_stiff_align = 16000;           % N*m/rad

% Steering damper coefficient (nominal)
P.c_steer = 6.0;                      % N*m*s/rad

% Suggested analysis speed range (m/s)
P.speed_range = [20 220] / 3.6;         % 2 to 70 km/h

check_units(P);
end

function check_units(P)
%CHECK_UNITS Validate basic magnitude of parameters.
%   CHECK_UNITS(P) throws an error when a parameter sits outside the
%   expected range, flagging likely unit mistakes. The thresholds are
%   intentionally generous so that reasonable variations pass.

arguments
    P struct
end

assert(P.g > 0 && P.g < 20, 'Gravity must be between 0 and 20 m/s^2.');
assert(P.wheelbase > 1.2 && P.wheelbase < 1.7, 'Wheelbase outside expected sport-tourer range.');
assert(P.steer_axis_tilt > deg2rad(15) && P.steer_axis_tilt < deg2rad(35), 'Steer axis tilt unrealistic.');
assert(P.trail > 0.05 && P.trail < 0.15, 'Trail should sit near 0.1 m.');
assert(P.front_radius > 0.25 && P.front_radius < 0.37, 'Front tyre radius invalid.');
assert(P.rear_radius > 0.25 && P.rear_radius < 0.37, 'Rear tyre radius invalid.');
assert(P.m_frame > 150 && P.m_frame < 300, 'Frame/rider mass unrealistic.');
assert(P.m_fork > 20 && P.m_fork < 60, 'Front assembly mass unrealistic.');
assert(P.x_frame > 0.4 && P.x_frame < 0.7, 'Frame CoM longitudinal coordinate out of bounds.');
assert(P.z_frame > 0.4 && P.z_frame < 0.7, 'Frame CoM height out of bounds.');
assert(P.x_fork > 1.2 && P.x_fork < 1.5, 'Fork CoM longitudinal coordinate out of bounds.');
assert(P.z_fork > 0.5 && P.z_fork < 0.9, 'Fork CoM height out of bounds.');
assert(all(structfun(@(x) x > 0, P.I_front_wheel)), 'Front wheel inertias must be positive.');
assert(P.tyre_stiff_lat > 1e5 && P.tyre_stiff_lat < 3e5, 'Lateral tyre stiffness outside expected range.');
assert(P.tyre_stiff_align > 5e3 && P.tyre_stiff_align < 4e4, 'Aligning stiffness outside expected range.');
assert(P.c_steer >= 0 && P.c_steer < 50, 'Steering damper coefficient unrealistic.');
assert(numel(P.speed_range) == 2 && all(P.speed_range > 0) && P.speed_range(2) > P.speed_range(1), 'Speed range must be a positive increasing pair.');
end
