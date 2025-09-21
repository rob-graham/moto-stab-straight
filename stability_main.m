function stability_main()
%STABILITY_MAIN  Main driver for motorcycle straight-running stability.
%   STABILITY_MAIN orchestrates the analysis: it loads a sport-tourer
%   parameter set, constructs the Whipple/Sharp matrices, scans the eigen
%   spectrum across forward speed, generates plots, and reports stability
%   crossings. A sensitivity sweep on the steering damper coefficient is
%   included to highlight its impact on the wobble and weave modes.
%
%   The script produces console output summarising the most critical
%   metrics and creates several figures (frequency, damping, real parts).
%
%   See also PARAMS_SPORT_TOURER, BUILD_MATRICES, EIGEN_SCAN, PLOTS_MODES.

P = params_sport_tourer();
fprintf('Loaded sport-tourer parameters (wheelbase %.2f m, head angle %.1f deg).\n', ...
    P.wheelbase, rad2deg(P.steer_axis_tilt));

mats = build_matrices(P);
speeds = linspace(P.speed_range(1), P.speed_range(2), 120);
results = eigen_scan(P, mats, speeds);
plots_modes(results);

summary = mode_summary(results);
disp('Nominal modal summary:');
disp(summary);

cross = stability_crossings(results);
if isempty(cross)
    disp('No unstable crossings detected within the analysed speed range.');
else
    disp('Neutral stability crossings (km/h):');
    disp(struct2table(cross));
end

% Steering damper sensitivity
c_values = [0, 2, 4, 6, 8, 10];
[scan_table, scans] = steering_sensitivity(P, c_values, speeds);
disp('Steering damper sensitivity (crossing speed in km/h, NaN if stable):');
disp(scan_table);

% Optionally plot damping sensitivity for wobble mode
figure('Name', 'Steering damper sensitivity'); hold on; grid on;
colours = lines(numel(c_values));
for k = 1:numel(c_values)
    name = sprintf('c=%.1f', c_values(k));
    wobble = scans{k}.modes.wobble;
    plot(speeds * 3.6, real(wobble), 'Color', colours(k, :), 'LineWidth', 1.3, ...
        'DisplayName', name);
end
xlabel('Speed (km/h)'); ylabel('Re(\lambda_{wobble}) (1/s)');
title('Effect of steering damper on wobble stability');
yline(0, '--k');
legend('Location', 'best');

end

function summary = mode_summary(results)
%MODE_SUMMARY Compute key figures for each classified mode.

names = results.mode_order;
num_modes = numel(names);
metrics = struct('Mode', strings(num_modes, 1), ...
    'MinDamping', zeros(num_modes, 1), ...
    'MaxRealPart', zeros(num_modes, 1), ...
    'CriticalSpeedKmH', zeros(num_modes, 1));

for k = 1:num_modes
    name = names{k};
    metrics.Mode(k) = string(name);
    damping = results.damping_ratio.(name);
    metrics.MinDamping(k) = min(damping);
    realpart = real(results.modes.(name));
    metrics.MaxRealPart(k) = max(realpart);
    idx = find(realpart >= 0, 1, 'first');
    if isempty(idx)
        metrics.CriticalSpeedKmH(k) = NaN;
    else
        metrics.CriticalSpeedKmH(k) = results.speeds(idx) * 3.6;
    end
end

summary = struct2table(metrics);
end

function cross = stability_crossings(results)
%STABILITY_CROSSINGS Identify speeds where modes cross the imaginary axis.

names = results.mode_order;
cross = struct('Mode', {}, 'SpeedKmH', {});
for k = 1:numel(names)
    name = names{k};
    val = real(results.modes.(name));
    speed = find_crossing(results.speeds, val);
    if ~isnan(speed)
        entry.Mode = string(name);
        entry.SpeedKmH = speed * 3.6;
        cross(end+1) = entry; %#ok<AGROW>
    end
end
end

function [table_out, scans] = steering_sensitivity(P, c_values, speeds)
%STEERING_SENSITIVITY Sweep steering damper settings.

scans = cell(numel(c_values), 1);
rows = struct('Damper', [], 'WobbleCross', [], 'WeaveCross', []);
rows = repmat(rows, numel(c_values), 1);
for k = 1:numel(c_values)
    Pk = P;
    Pk.c_steer = c_values(k);
    mats = build_matrices(Pk);
    scans{k} = eigen_scan(Pk, mats, speeds);
    wobble_cross = find_crossing(speeds, real(scans{k}.modes.wobble));
    weave_cross = find_crossing(speeds, real(scans{k}.modes.weave));
    rows(k).Damper = c_values(k);
    rows(k).WobbleCross = wobble_cross * 3.6;
    rows(k).WeaveCross = weave_cross * 3.6;
end

table_out = struct2table(rows);
end

function speed = find_crossing(speeds, values)
%FIND_CROSSING Locate a zero crossing via linear interpolation.

sign_change = values(1:end-1) .* values(2:end);
idx = find(sign_change <= 0, 1, 'first');
if isempty(idx)
    speed = NaN;
    return;
end
v1 = values(idx);
v2 = values(idx+1);
if abs(v2 - v1) < eps
    speed = speeds(idx);
else
    speed = interp1([v1, v2], [speeds(idx), speeds(idx+1)], 0, 'linear', 'extrap');
end
end
