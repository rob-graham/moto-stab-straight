function impulse_sim(P, mats, speed, t_final, steer_impulse)
%IMPULSE_SIM  Simulate lean/steer response to a steering impulse.
%   IMPULSE_SIM(P, MATS, SPEED, T_FINAL, STEER_IMPULSE) integrates the
%   linear state dynamics (constructed from MATS at the forward SPEED) and
%   plots the resulting lean angle, steer angle, and steer rate. The initial
%   condition corresponds to an impulse in steering rate with magnitude
%   STEER_IMPULSE (default 1 rad/s). T_FINAL specifies the simulation
%   horizon in seconds (default 3 s).

if nargin < 4 || isempty(t_final)
    t_final = 3;
end
if nargin < 5 || isempty(steer_impulse)
    steer_impulse = 1;
end

M = mats.M;
B = speed * mats.C1 + mats.D;
K = P.g * mats.K0 + speed^2 * mats.K2;
A = [zeros(2), eye(2); -M\K, -M\B];

x0 = [0; 0; 0; steer_impulse];

[t, x] = ode45(@(t, x) A * x, linspace(0, t_final, 400), x0);

phi = x(:, 1);
delta = x(:, 2);
delta_dot = x(:, 4);

figure('Name', 'Impulse response');
subplot(3,1,1); plot(t, phi, 'LineWidth', 1.4); grid on;
ylabel('\phi (rad)'); title(sprintf('Lean and steer response, v = %.1f km/h', speed*3.6));
subplot(3,1,2); plot(t, delta, 'LineWidth', 1.4); grid on;
ylabel('\delta (rad)');
subplot(3,1,3); plot(t, delta_dot, 'LineWidth', 1.4); grid on;
ylabel('\dot{\delta} (rad/s)'); xlabel('Time (s)');
end
