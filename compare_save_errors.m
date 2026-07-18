%% compare_and_save_errors.m
% ------------------------------------------------------------------------
% Compares the area-weighted average of the six individual material runs
% ("expected") against one equivalent combined model, for a given areal
% ratio, depth and resolution. Plots the surface time series, the depth-time
% temperature profiles, and appends the error metrics to a database for
% later analysis.
% ------------------------------------------------------------------------
close all; clear; clc;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end

% where the model results live, the forcing file, and the error database
results_dir   = fullfile(script_dir, 'Results', '2026-07-11');
forcing_file  = fullfile(script_dir, 'environmental_forcing.mat');
error_db_file = fullfile(script_dir, 'Results', 'error_metrics.mat');

% areal-ratio matrix (must match the one used to run the model)
A = [ ...
   0.1500 0.1500 0.2000 0.1500 0.2000 0.1500
   0.1400 0.1400 0.2400 0.1400 0.2000 0.1400
   0.1300 0.1300 0.2100 0.1900 0.1700 0.1700
   0.2000 0.1500 0.1500 0.1000 0.2500 0.1500
        0 0.4000      0 0.3500 0.2500      0
   0.5000      0      0 0.1000 0.4000      0
        0      0 0.6000 0.4000      0      0
   0.3300 0.3300 0.3400      0      0      0
   0.1000 0.2000 0.2000      0 0.3000 0.2000
   0.1000 0.2000 0.3000 0.1000 0.2000 0.1000
   0.1200 0.2800 0.2000 0.1000 0.1500 0.1500
   0.1600 0.1400 0.2000 0.1600 0.2400 0.1000
   0.1700 0.1700 0.1800 0.1800 0.1500 0.1500
   0.3000 0.2500      0 0.2000 0.2500      0
        0 0.6000      0 0.4000      0      0
   0.4000      0      0      0 0.6000      0
   0.2700 0.3000 0.2500      0      0 0.1800
        0 0.3000 0.4000 0.3000      0      0
   0.2200 0.2200 0.1200 0.1100 0.2000 0.1300
   0.1100 0.2100 0.1900 0.1200 0.1700 0.2000
        0      0 0.3000 0.3000 0.4000      0
   0.2600      0 0.2400      0 0.3000 0.2000
   0.1800 0.1200 0.2000 0.1000 0.2500 0.1500
        0      0      0 0.7000 0.3000      0
   0.4000 0.3000      0      0 0.3000      0
   0.1000 0.3000 0.2000      0 0.4000      0
   0.1500 0.1500 0.1500 0.1500 0.2500 0.1500
   0.2000 0.1000 0.2000 0.2000 0.1000 0.2000
   0.3000      0 0.1000 0.3000 0.3000      0
        0 0.3300 0.3300      0 0.3400      0
        0      0 0.5000 0.2500      0 0.2500
   0.1900 0.1100 0.2000 0.2000 0.1500 0.1500
   0.2000 0.2000 0.2000 0.1000 0.1000 0.2000
   0.1000 0.1000 0.3000 0.1000 0.1000 0.3000
   0.2100 0.1900 0.1800 0.1200 0.1500 0.1500
   0.1000 0.2000 0.3000 0.2000 0.2000      0
        0 0.1500 0.2500 0.1500 0.2500 0.2000
   0.1000 0.1000 0.1000 0.1000 0.6000      0
   0.2000 0.3000      0 0.2000 0.3000      0
   0.2500      0 0.2500      0 0.2500 0.2500
        0 0.2000 0.5000      0 0.3000      0
        0      0 0.3300 0.3300 0.3400      0
   0.4500      0      0      0 0.4500 0.1000
        0 0.5000      0      0 0.5000      0
   0.7000      0      0 0.3000      0      0
        0      0 0.6000      0      0 0.4000
        0      0 0.5000      0      0 0.5000
   0.1000 0.1000 0.1000      0 0.7000      0
        0 0.2000 0.2000 0.3000 0.3000      0
   0.1000 0.1000 0.2000 0.3000 0.3000      0];

% scenario selection
jjj         = 10;        % areal-ratio row of A
aaaa        = 2;         % index into matrix_dx
depth       = 0.50;      % domain depth (m)
eq_material = 7;        % which combined model is the equivalent (7..16)

figures = 1;

matrix_dx = [0.005, 0.01, 0.02, 0.03, 0.05, 0.06, 0.1, 0.15];



base_names = {'Concrete Pavement','Asphalt Pavement','Stone', ...
              'Gravel with fines','Paving Clay Tiles','Rubberized Pavement'};

%% --- load the six individual material runs ------------------------------
 for eq_material=7:8

     dz = matrix_dx(aaaa);
     a  = A(jjj, :);
T_profile_collection = [];
Ts_collection        = [];
G_flux_collection    = [];
time_vector          = [];
for i = 1:6
    p = scenario_path(results_dir, i, a, depth, dz);
    assert(isfile(p), 'Missing individual run: %s', p);
    S = load(p, 'T_profile', 'G_flux', 'time_vector');
    if isempty(time_vector), time_vector = S.time_vector(:).'; end
    if isempty(T_profile_collection)
        T_profile_collection = zeros([size(S.T_profile), 6]);
    end
    T_profile_collection(:,:,i) = S.T_profile;
    Ts_collection(i,:)          = S.T_profile(1,:);
    G_flux_collection(i,:)      = S.G_flux(:).';
end

n_steps    = numel(time_vector);
time_hours = time_vector / 3600;

%% --- area-weighted average of the individuals ("expected") --------------
Ts_avg        = zeros(1, n_steps);
G_flux_avg    = zeros(1, n_steps);
T_profile_avg = zeros(size(T_profile_collection(:,:,1)));
for i = 1:6
    Ts_avg        = Ts_avg        + a(i) * Ts_collection(i,:);
    G_flux_avg    = G_flux_avg    + a(i) * G_flux_collection(i,:);
    T_profile_avg = T_profile_avg + a(i) * T_profile_collection(:,:,i);
end

%% --- load the equivalent combined model ---------------------------------
eq_path = scenario_path(results_dir, eq_material, a, depth, dz);
assert(isfile(eq_path), 'Missing equivalent run: %s', eq_path);
E = load(eq_path, 'T_profile', 'G_flux', 'meta');
T_profile_eq = E.T_profile;
Ts_eq        = E.T_profile(1,:);
G_flux_eq    = E.G_flux(:).';
eq_name      = E.meta.material_name;

T_difference = T_profile_eq - T_profile_avg;

%% --- day / night split from the shortwave forcing -----------------------
F = load(forcing_file, 'SW_dir_ground', 'SW_diff_ground');
SW_all_g = compute_SW(F.SW_dir_ground, F.SW_diff_ground, n_steps);
day = SW_all_g(:).' > 0;
SunriseIdx = find(day, 1, 'first');
SunsetIdx  = find(day, 1, 'last');
if isempty(SunriseIdx), SunriseIdx = 1; SunsetIdx = n_steps; end
day_idx   = SunriseIdx:SunsetIdx;
night_idx = [1:SunriseIdx, SunsetIdx:n_steps];

%% --- errors: equivalent vs expected average -----------------------------
temp_error = Ts_eq - Ts_avg;
flux_error = G_flux_eq - G_flux_avg;

% material list for the figure titles
entries = arrayfun(@(i) sprintf('%s (%.2f)', base_names{i}, a(i)), 1:6, 'UniformOutput', false);
line1 = strjoin(entries(1:3), ', ');
line2 = strjoin(entries(4:6), ', ');

if figures == 1
%% ===== Figure 1: surface temperature and ground heat flux ===============
figure('Color','w','Name','Average vs Equivalent');
colors = lines(2);
sgtitle(sprintf('Averaged individuals vs equivalent model\ndz = %.2f cm\n%s\n%s', ...
                dz*100, line1, line2));

ax1 = subplot(2,1,1); hold(ax1,'on'); grid(ax1,'on');
plot(ax1, time_hours, Ts_avg, '--', 'Color', colors(1,:), 'LineWidth', 2, ...
     'DisplayName', 'average of individual');
plot(ax1, time_hours, Ts_eq, '-', 'Color', colors(2,:), 'LineWidth', 1.5, ...
     'DisplayName', sprintf('%s (MAE %.2f C, RMSE %.2f C)', ...
                            eq_name, mean(abs(temp_error)), sqrt(mean(temp_error.^2))));
ylabel(ax1, 'Surface temperature (\circC)'); xlabel(ax1, 'Time (hours)');
legend(ax1, 'Location', 'best');

ax2 = subplot(2,1,2); hold(ax2,'on'); grid(ax2,'on');
plot(ax2, time_hours(2:end), G_flux_avg(2:end), '--', 'Color', colors(1,:), 'LineWidth', 2, ...
     'DisplayName', 'average of individual');
plot(ax2, time_hours(2:end), G_flux_eq(2:end), '-', 'Color', colors(2,:), 'LineWidth', 1.5, ...
     'DisplayName', sprintf('%s (MAE %.2f, RMSE %.2f W/m^2)', ...
                            eq_name, mean(abs(flux_error)), sqrt(mean(flux_error.^2))));
ylabel(ax2, 'Ground heat flux (W/m^2)'); xlabel(ax2, 'Time (hours)');
legend(ax2, 'Location', 'best');

%% ===== Figure 2: depth-time profiles ====================================
num_layers = size(T_profile_avg, 1) - 1;
z_nodes    = (0:num_layers) * dz * 100;     % cm, surface at 0

figure('Color','w','Name','Depth-Time Profiles','Position',[100 100 800 900]);
clim_shared = [min([T_profile_avg(:); T_profile_eq(:)]), ...
               max([T_profile_avg(:); T_profile_eq(:)])];

ax1 = subplot(3,1,1);
imagesc(ax1, time_hours, z_nodes, T_profile_avg); set(ax1,'YDir','reverse');
clim(ax1, clim_shared); colorbar(ax1);
title(ax1,'Expected (area-weighted average of individuals)');
ylabel(ax1,'Depth (cm)'); xlabel(ax1,'Time (hours)');

ax2 = subplot(3,1,2);
imagesc(ax2, time_hours, z_nodes, T_profile_eq); set(ax2,'YDir','reverse');
clim(ax2, clim_shared); colorbar(ax2);
title(ax2,'Equivalent system');
ylabel(ax2,'Depth (cm)'); xlabel(ax2,'Time (hours)');

ax3 = subplot(3,1,3);
imagesc(ax3, time_hours, z_nodes, T_difference); set(ax3,'YDir','reverse');
cb = colorbar(ax3); ylabel(cb,'Temperature difference (\circC)');
colormap(ax3, flipud(hot));
title(ax3,'Difference (equivalent - expected)');
ylabel(ax3,'Depth (cm)'); xlabel(ax3,'Time (hours)');

%% ===== Figure 3: individual base scenario profiles ======================
figure('Color','w','Name','Individual Base Scenario Profiles','Position',[950 100 900 800]);
clim_ind = [min(T_profile_collection(:)), max(T_profile_collection(:))];
for i = 1:6
    subplot(3,2,i);
    imagesc(time_hours, z_nodes, T_profile_collection(:,:,i)); set(gca,'YDir','reverse');
    clim(clim_ind); colorbar;
    title(base_names{i}); xlabel('Time (hours)'); ylabel('Depth (cm)');
end
sgtitle('Base scenario temperature profiles');
end
%% ===== error metrics -> append to database ==============================
S = struct();
S.mbe_temp  = mean(temp_error);
S.mae_temp  = mean(abs(temp_error));
S.rmse_temp = sqrt(mean(temp_error.^2));
S.mbe_flux  = mean(flux_error);
S.mae_flux  = mean(abs(flux_error));
S.rmse_flux = sqrt(mean(flux_error.^2));

S.mbe_temp_day  = mean(temp_error(day_idx));
S.mae_temp_day  = mean(abs(temp_error(day_idx)));
S.rmse_temp_day = sqrt(mean(temp_error(day_idx).^2));
S.mbe_flux_day  = mean(flux_error(day_idx));
S.mae_flux_day  = mean(abs(flux_error(day_idx)));
S.rmse_flux_day = sqrt(mean(flux_error(day_idx).^2));

S.mbe_temp_night  = mean(temp_error(night_idx));
S.mae_temp_night  = mean(abs(temp_error(night_idx)));
S.rmse_temp_night = sqrt(mean(temp_error(night_idx).^2));
S.mbe_flux_night  = mean(flux_error(night_idx));
S.mae_flux_night  = mean(abs(flux_error(night_idx)));
S.rmse_flux_night = sqrt(mean(flux_error(night_idx).^2));

S.surface_ratios = a;
S.resolution     = dz;
S.depth          = depth;      
S.n_dz           = n_steps;
S.eq_material    = eq_material;

if ~exist(fileparts(error_db_file),'dir'), mkdir(fileparts(error_db_file)); end
if isfile(error_db_file)
    L = load(error_db_file, 'Data');
    if isfield(L,'Data') && isstruct(L.Data), Data = L.Data; else, Data = struct([]); end
else
    Data = struct([]);
end

% overwrite a matching entry if present, otherwise append
tol = 1e-6; found = false;
for i = 1:numel(Data)
    if all(abs(Data(i).surface_ratios - S.surface_ratios) < tol) && ...
       abs(Data(i).resolution - S.resolution) < tol && ...
       Data(i).n_dz == S.n_dz && Data(i).eq_material == S.eq_material
        Data(i) = S; found = true; break;
    end
end
if ~found
    if isempty(Data), Data = S; else, Data(end+1) = S; end
end
save(error_db_file, 'Data');
fprintf('Saved error metrics to %s (%d entries).\n', error_db_file, numel(Data));
 end
 
%% ===== local functions ==================================================
function p = scenario_path(results_dir, material, a, depth, dz)
% rebuild the filename the model wrote
if material > 6
    ratio_tag = ['_a' strjoin(compose('%02d', round(a*100)), '-')];
else
    ratio_tag = '';
end
name = sprintf('m%02d%s_d%.2f_dz%.3f', material, ratio_tag, depth, dz);
p = fullfile(results_dir, [name '.mat']);
end

function SW_all_g = compute_SW(dir_ground_native, diff_ground_native, n_model)
idx   = linspace(1, numel(dir_ground_native), n_model).';
Dirrg = interp1(dir_ground_native,  idx, 'linear');
Diffg = interp1(diff_ground_native, idx, 'linear');
SW_all_g = Dirrg + Diffg;
end