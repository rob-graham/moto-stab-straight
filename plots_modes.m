function plots_modes(results)
%PLOTS_MODES  Plot modal characteristics versus forward speed.
%   PLOTS_MODES(RESULTS) creates three figures summarising the eigenvalue
%   data produced by EIGEN_SCAN: modal frequency, damping ratio, and real
%   part of the eigenvalues. Colours and labels remain consistent across
%   the plots to ease interpretation.
%
%   The forward speed axis is displayed in km/h. Frequencies are given in
%   Hz, damping ratios are dimensionless, and the real parts are in 1/s.

arguments
    results struct
end

speeds_kmh = results.speeds * 3.6;
mode_names = results.mode_order;
colours = get_mode_colours();

% Frequency plot
figure('Name', 'Modal frequency'); hold on; grid on;
for k = 1:numel(mode_names)
    name = mode_names{k};
    plot(speeds_kmh, results.frequency_hz.(name), 'LineWidth', 1.5, ...
        'Color', colours.(name));
end
xlabel('Speed (km/h)'); ylabel('Frequency (Hz)');
title('Natural frequency of modes');
legend(mode_names, 'Location', 'best');

% Damping ratio plot
figure('Name', 'Modal damping'); hold on; grid on;
for k = 1:numel(mode_names)
    name = mode_names{k};
    plot(speeds_kmh, results.damping_ratio.(name), 'LineWidth', 1.5, ...
        'Color', colours.(name));
end
xlabel('Speed (km/h)'); ylabel('Damping ratio');
title('Modal damping ratio vs speed');
yline(0, '--k', 'Zero damping');
legend(mode_names, 'Location', 'best');

% Real part plot
figure('Name', 'Modal real part'); hold on; grid on;
for k = 1:numel(mode_names)
    name = mode_names{k};
    plot(speeds_kmh, real(results.modes.(name)), 'LineWidth', 1.5, ...
        'Color', colours.(name));
end
xlabel('Speed (km/h)'); ylabel('Real part (1/s)');
title('Eigenvalue real part vs speed');
yline(0, '--k', 'Neutral');
legend(mode_names, 'Location', 'best');
end

function colours = get_mode_colours()
colours.wobble = [0.8500, 0.3250, 0.0980];
colours.weave  = [0.0000, 0.4470, 0.7410];
colours.capsize = [0.4940, 0.1840, 0.5560];
colours.caster = [0.4660, 0.6740, 0.1880];
end
