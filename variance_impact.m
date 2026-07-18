%% variance_impact.m
% ------------------------------------------------------------------------
% Error vs area-weighted property variance, for ONE approach, ONE depth and
% ONE resolution, using every design (areal ratio) available under those.
% ------------------------------------------------------------------------
close all; clear; clc;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
error_db_file = fullfile(script_dir, 'Results', 'error_metrics.mat');

assert(isfile(error_db_file), 'Error database not found: %s', error_db_file);
load(error_db_file, 'Data');
assert(exist('Data','var') == 1 && isstruct(Data), 'Invalid Data structure.');
assert(isfield(Data,'eq_material'), 'Database has no approach field. Rebuild it with the updated compare_and_save_errors.m.');
assert(isfield(Data,'depth'),       'Database has no depth field. Rebuild it with the updated compare_and_save_errors.m.');

% base material properties (concrete, asphalt, stone, gravel, clay, rubber)
base.albedo = [0.4,     0.08,    0.3,     0.25,    0.3,     0.1];
base.K      = [1.6,     0.75,    2.5,     0.5,     1.0,     0.3];
base.C      = [2700000, 2000000, 2400000, 1500000, 2000000, 1300000];

% --- filters: pick one approach, one depth, one resolution ---
approach_target   = 7;     % combination model index (the approach); do not mix approaches
depth_target      = 0.50;   % m
resolution_target = 0.01;   % m

tol = 1e-6;
sel = ([Data.eq_material] == approach_target) ...
    & (abs([Data.depth]      - depth_target)      < tol) ...
    & (abs([Data.resolution] - resolution_target) < tol);
D = Data(sel);
assert(~isempty(D), 'No entries for approach %d at depth %.2f m, resolution %.4f m.', ...
       approach_target, depth_target, resolution_target);

approach_name = 'Present_study';  

% per-design property variance and MAE
N = numel(D);
var_albedo = zeros(N,1); var_K = zeros(N,1); var_C = zeros(N,1);
mae_temp   = zeros(N,1); mae_flux = zeros(N,1);
for i = 1:N
    a = D(i).surface_ratios;
    var_albedo(i) = sum(a .* base.albedo.^2) - sum(a .* base.albedo)^2;
    var_K(i)      = sum(a .* base.K.^2)      - sum(a .* base.K)^2;
    var_C(i)      = sum(a .* base.C.^2)      - sum(a .* base.C)^2;
    mae_temp(i)   = D(i).mae_temp;
    mae_flux(i)   = D(i).mae_flux;
end

var_data  = {var_albedo, var_K, var_C};
var_names = {'Variance of albedo', 'Variance of K', 'Variance of C'};
var_units = {'', '(W/m-K)^2', '(J/m^3-K)^2'};

figure('Color','w', ...
       'Name', sprintf('MAE vs Property Variance | %s | depth %.2f m | dz %.2f cm', ...
                       approach_name, depth_target, resolution_target*100), ...
       'Position',[50 50 1200 800]);
for i = 1:3
    ax = subplot(3,2,(i-1)*2 + 1); hold(ax,'on'); grid(ax,'on');
    scatter(ax, var_data{i}, mae_temp, 36, 'filled');
    r2 = plot_regression(ax, var_data{i}, mae_temp);
    title(ax, ['T MAE vs ' var_names{i}]);
    xlabel(ax, [var_names{i} ' ' var_units{i}]); ylabel(ax, 'MAE-T (\circC)');
    if ~isnan(r2)
        text(ax, 0.03, 0.97, sprintf('R^2 = %.3f', r2), 'Units','normalized', ...
             'VerticalAlignment','top', 'FontWeight','bold');
    end

    ax = subplot(3,2,(i-1)*2 + 2); hold(ax,'on'); grid(ax,'on');
    scatter(ax, var_data{i}, mae_flux, 36, 'filled');
    r2 = plot_regression(ax, var_data{i}, mae_flux);
    title(ax, ['G MAE vs ' var_names{i}]);
    xlabel(ax, [var_names{i} ' ' var_units{i}]); ylabel(ax, 'MAE-G (W/m^2)');
    if ~isnan(r2)
        text(ax, 0.03, 0.97, sprintf('R^2 = %.3f', r2), 'Units','normalized', ...
             'VerticalAlignment','top', 'FontWeight','bold');
    end
end

%% ===== local function ===================================================
function r_sq = plot_regression(ax, x, y)
r_sq = NaN;
ok = ~isnan(x(:)) & ~isnan(y(:));
x = x(ok); y = y(ok);
if numel(x) < 2, return; end
p  = polyfit(x, y, 1);
xl = linspace(min(x), max(x), 10);
plot(ax, xl, polyval(p, xl), 'LineWidth', 1.5, 'HandleVisibility','off');
yp = polyval(p, x);
ss_tot = sum((y - mean(y)).^2);
ss_res = sum((y - yp).^2);
if ss_tot == 0, r_sq = 1; else, r_sq = 1 - ss_res/ss_tot; end
end