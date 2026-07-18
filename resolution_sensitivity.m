%% resolution_sensitivity.m
% ------------------------------------------------------------------------
% Error vs vertical resolution for ONE approach, ONE depth and ONE design
% (areal ratio), using every resolution available under those.
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

% --- filters: pick one approach, one depth, one design ---
approach_target = 7;     % combination model index (the approach); do not mix approaches
depth_target    = 0.50;   % m
target_ratios   = [0.1000 0.2000 0.3000 0.1000 0.2000 0.1000];   % the one design

tol = 1e-6;
sel = ([Data.eq_material] == approach_target) ...
    & (abs([Data.depth] - depth_target) < tol) ...
    & arrayfun(@(d) all(abs(d.surface_ratios - target_ratios) < tol), Data);
D = Data(sel);
assert(~isempty(D), 'No entries for this approach / depth / design.');

approach_name = 'Present Study';   % readable label for titles

% group the matching entries by resolution
res_cm = round([D.resolution]*100 * 1e5) / 1e5;
[unique_res, ~, gidx] = unique(res_cm);
num_res = numel(unique_res);
if num_res < 2
    warning('Only one resolution found for this approach/depth/design; no trend to show.');
end

error_types = {'rmse','mbe','mae'};
variables   = {'temp','flux'};
time_suffix = {'', '_day', '_night'};
time_label  = {'whole day (24h)', 'daytime', 'nighttime'};

% mean error per resolution (handles duplicates; usually one entry each)
err = struct();
for t = 1:numel(time_suffix)
    for v = 1:numel(variables)
        for e = 1:numel(error_types)
            f = sprintf('%s_%s%s', error_types{e}, variables{v}, time_suffix{t});
            vals = nan(num_res,1);
            for r = 1:num_res
                vals(r) = mean(arrayfun(@(d) d.(f), D(gidx == r)), 'omitnan');
            end
            err.(f) = vals;
        end
    end
end

% one figure per time window, temperature on top and flux below
for t = 1:numel(time_suffix)
    figure('Color','w', ...
           'Name', sprintf('%s | %s | depth %.2f m', time_label{t}, approach_name, depth_target));
    tiledlayout(2,1);
    for v = 1:numel(variables)
        nexttile; hold on; grid on;
        for e = 1:numel(error_types)
            f = sprintf('%s_%s%s', error_types{e}, variables{v}, time_suffix{t});
            plot(unique_res, err.(f), '-o', 'LineWidth', 1.5, 'DisplayName', upper(error_types{e}));
        end
        title(sprintf('%s — %s', time_label{t}, variables{v}));
        xlabel('Resolution (cm)'); ylabel('Error'); legend('Location','best');
    end
end