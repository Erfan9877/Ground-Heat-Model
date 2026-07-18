%% ground_heat_model.m
% ------------------------------------------------------------------------
%
% Author: Erfan Hosseini (email:eh9561@princeton.edu), May 28, 2026
%
% ------------------------------------------------------------------------
% 1D vertical heat-transfer model for paved/urban surfaces.
%
% Solves the surface energy balance coupled to vertical heat conduction in
% the ground, for either a single material or a heterogeneous patch built
% from several materials via effective-medium combination rules (each with
% its own areal ratio).
%
% Workflow:
%   - environmental forcing is loaded from environmental_forcing.mat
%   - the data time step and the model time step are determined dynamically
%     and the forcing is interpolated onto the model grid
%   - thermal properties (k, C, albedo) are set per scenario
%   - the model is run (spin-up to convergence + explicit finite-difference
%     solver) for the requested material(s), areal ratio(s), domain depth(s)
%     and vertical resolution(s)
%   - results are written to ./Results/<date>/<scenario>.mat, with the
%     scenario name encoding material/combination, areal ratio, domain
%     depth and resolution
%
% Revision (Jul 2026): added material index 8 = Salamanca et al. (2009,
%   JAMC, 48, 1725-1732)
%
% Inputs : environmental_forcing.mat
% Outputs: per-scenario .mat files containing T_profile, G_flux, depth,
%          time vector and run metadata
% ------------------------------------------------------------------------
close all; clear; clc;

%% =================== 0. paths & load forcing ============================
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end

forcing_file = fullfile(script_dir, 'environmental_forcing.mat');
assert(isfile(forcing_file), ...
    'environmental_forcing.mat not found.');
F = load(forcing_file);

%% =================== 1. model time settings (dynamic) ==================
% --- pick the model time step here (in seconds).
model_dt_seconds = 10;                          % e.g. 10 s

total_time   = F.total_time_hours;              % hours, taken from the data
delta_t      = model_dt_seconds / 60;           % minutes 
dt           = model_dt_seconds;                % seconds 
n_model      = round(total_time * 3600 / dt);   % number of model steps
time         = (0:n_model-1).' * dt;            % model time vector

% native data step, reported for transparency
dt_data_env  = F.dt_env;                        % s
dt_data_sw   = F.dt_sw;                          % s
fprintf('Data dt (forcing) = %.2f s | Data dt (SW) = %.2f s | Model dt = %.2f s | %d steps\n', ...
        dt_data_env, dt_data_sw, dt, n_model);

% index-based linear interpolation onto the model grid (identical in
% behaviour to the original interp1(var, linspace(1,N,n_model)) call, but
% now n_model is derived from whatever model dt you chose).
to_model = @(v) interp1(v, linspace(1, numel(v), n_model).', 'linear');

%% =================== 2. interpolate forcing to model grid ==============
Ux            = to_model(F.Ux);
Uy            = to_model(F.Uy);
Uz            = to_model(F.Uz);
SWin_measured = to_model(F.SWin_measured);
LWin_measured = to_model(F.LWin_measured);
RH            = to_model(F.RH);
TPT           = to_model(F.TPT);
SBTC          = to_model(F.SBTC);
TMV           = to_model(F.TMV);
CO2           = to_model(F.CO2);
P             = to_model(F.P);
CO2_gm3       = to_model(F.CO2_gm3);
Tsonic        = to_model(F.Tsonic);
H2O           = to_model(F.H2O);
LWin=LWin_measured;

hw = ((Ux.^2) + (Uy.^2)).^0.5;

% shortwave on the model grid 
SW_all_g = compute_SW(F.SW_dir_ground, F.SW_diff_ground, n_model);

%% =================== 3. general constants ===============================
sigma   = 5.67037e-8;
k_vk    = 0.41;     % von Karman's constant 
Rv      = 8.31;     roh = 1.204;   Cp  = 1005;   Lv = 2450000;
kair    = 0.025;    kinvisc = 1.55e-5;   Pr = 1;
z_wind  = 6;        z0  = 0.01;
hw(:,:) = 1;        % (keep hw with 1, so that no impact from H)

% sensible-heat coefficient 
h_s = (roh * Cp * k_vk^2 / Pr) .* (hw ./ (log(z_wind ./ z0)).^2);

% air temperature, vapour pressure, incoming longwave 
qa = H2O * 0.000014694;
Ta = Tsonic ./ (1 + (0.51 * qa));
Ta = Ta + 10;
ea = 611 * (RH/100) .* exp((17.27*Ta)./(Ta+237.3));     % Pa


%% =================== 4. material / geometry tables ======================
matrix    = 1:16;
matrix_dx = [0.005, 0.01, 0.02, 0.03, 0.05, 0.06, 0.1, 0.15];

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

% base material properties
C_material      = [2700000, 2000000, 2400000, 1500000, 2000000, 1300000];
k_material      = [1.6,     0.75,    2.5,     0.5,     1,       0.3];
albedo_material = [0.4,     0.08,    0.3,     0.25,    0.3,     0.1];
es = 0.95;                                   % surface emissivity

plots = 1;

%% =================== 5. run requested scenarios =========================
% --- choose what to run here  -----------
depth_list = 0.5; % 0.2:0.2:1;        % domain depth(s)  [m]
A_rows     = 10; % 1:size(A,1);                 % which area-fraction row(s) of A
dz_idx     = 2;                  % which entries of matrix_dx
mat_idx    = 7; % 1:7;                % which materials/combinations

for domain_depth = depth_list
for jjj = A_rows
for aaaa = dz_idx
for aaa = mat_idx

    a        = A(jjj, :);
    material = matrix(aaa);

    % ---- scenario labels -------
    formatted_str   = ['resistance-weighted C (' strtrim(sprintf('%.2f ', a)) ')'];
    formatted_str2  = ['Carson 2022, k harmonic, c arithmic (' strtrim(sprintf('%.2f ', a)) ')'];
    formatted_str3  = ['mixture  theory (KBS) both arithmic (' strtrim(sprintf('%.2f ', a)) ')'];
    formatted_str4  = ['Urban Mosaic, Sugawara, k and kc arithmic (' strtrim(sprintf('%.2f ', a)) ')'];
    formatted_str5  = ['Multiphase Maxwell-Eucken, c arithmic, k inclusion (' strtrim(sprintf('%.2f ', a)) ')'];
    formatted_str6  = ['Multiphase Hashin-Shtrikman (' strtrim(sprintf('%.2f ', a)) ')'];
    formatted_str8  = ['Multiphase ETM (' strtrim(sprintf('%.2f ', a)) ')'];
    formatted_str9  = ['Arithmic C, matching effusivities (' strtrim(sprintf('%.2f ', a)) ')'];
    formatted_str10 = ['Matching effusivity and (diffusivity)-1 (' strtrim(sprintf('%.2f ', a)) ')'];
    material_name = {'Concrete Pavement','Asphalt Pavement','Stone', ...
        'Gravel with fines','Paving Clay Tiles','Rubberized Pavement', ...
        formatted_str, formatted_str2, formatted_str3, formatted_str4, ...
        formatted_str5, formatted_str6, '', formatted_str8, ...
        formatted_str9, formatted_str10};

    % =================== MATERIAL PROPERTY TREE ===============
    if material <= 6
        % Single-material cases
        k_pave  = k_material(material);
        C_pave  = C_material(material);
        albedos = albedo_material(material);
        intended_time_scale = 30; % minutes
        length_scale_pave = sqrt(k_pave*(intended_time_scale*60)/C_pave);
        sprintf('the thermal diffusive length scale for this pavement in %.1f minutes is %.2f cm', ...
                intended_time_scale, length_scale_pave*100);
        kapa_pave = k_pave / C_pave;
        alpha = kapa_pave;
        if exist('domain_depth','var')
            F_Tg_depth = domain_depth;
        else
            F_Tg_depth = (2.222*((43200*kapa_pave)^(0.5)));
        end

    elseif material == 7
        F_Tg_depth1 = (2.273*((43200*(k_material(1)/C_material(1))))^(0.5));
        F_Tg_depth2 = (2.273*((43200*(k_material(2)/C_material(2))))^(0.5));
        F_Tg_depth3 = (2.273*((43200*(k_material(3)/C_material(3))))^(0.5));
        F_Tg_depth4 = (2.273*((43200*(k_material(4)/C_material(4))))^(0.5));
        F_Tg_depth5 = (2.273*((43200*(k_material(5)/C_material(5))))^(0.5));
        F_Tg_depth6 = (2.273*((43200*(k_material(6)/C_material(6))))^(0.5));
        F_Tg_depth_T = mean([F_Tg_depth1, F_Tg_depth2, F_Tg_depth3, ...
                             F_Tg_depth4, F_Tg_depth5, F_Tg_depth6]);
        if exist('domain_depth','var')
            F_Tg_depth = domain_depth;
        else
            F_Tg_depth = F_Tg_depth_T;
        end
        albedos = (a(1)*albedo_material(1)) + (a(2)*albedo_material(2)) + ...
                  (a(3)*albedo_material(3)) + (a(4)*albedo_material(4)) + ...
                  (a(5)*albedo_material(5)) + (a(6)*albedo_material(6));
        F = [F_Tg_depth1, F_Tg_depth2, F_Tg_depth3, F_Tg_depth4, F_Tg_depth5, F_Tg_depth6]; 
        C_pave = sum(a.*C_material);
        R   = [F_Tg_depth1, F_Tg_depth2, F_Tg_depth3, F_Tg_depth4, F_Tg_depth5, F_Tg_depth6] ./ k_material;
        num = sum((a.*1) ./ (R .* C_material));
        k_pave = F_Tg_depth_T*C_pave*num;

    elseif material == 8
        k_pave = sum(a ./ (k_material));
        C_eq   = sum(a.*C_material);
        C_pave = C_eq;

    elseif material == 9
        k_pave = sum(a .* (k_material));
        C_eq   = sum(a.*C_material);
        C_pave = C_eq;

    elseif material == 10
        k_pave = sum(a .* (k_material));
        C_eq   = (sum(a.*C_material.*k_material))/k_pave;
        C_pave = C_eq;

    elseif material == 11
        C_pave = sum(C_material .* a);
        % Maxwell-Eucken iterative
        k_pave = k_material(1);
        phi_total = a(1);
        for i = 2:length(k_material)
            phi = a(i) / (phi_total + a(i));
            k_incl = k_material(i);
            k_matrix = k_pave;
            k_pave = k_matrix * ...
                (k_incl + 2*k_matrix - 2*phi*(k_matrix - k_incl)) / ...
                (k_incl + 2*k_matrix + phi*(k_matrix - k_incl));
            phi_total = phi_total + a(i);
        end

    elseif material == 12 || material == 13
        C_pave = sum(C_material .* a);
        present_idx = a > 0;
        if ~any(present_idx)
            warning('Hashin-Shtrikman: No materials with area ratio > 0.');
            if material == 12
                k_pave = 0;
            end
            k_eff_low = 0;
            k_eff_high = 0;
        else
            k_present = k_material(present_idx);
            a_present = a(present_idx);
            a_norm = a_present / sum(a_present);

            % lower bound
            [k_sorted, idx] = sort(k_present);
            a_sorted = a_norm(idx);
            k_eff_low = k_sorted(1);
            phi_total = a_sorted(1);
            for i = 2:length(k_sorted)
                if k_sorted(i) == 0, continue; end
                k_i = k_sorted(i);
                if (phi_total + a_sorted(i)) == 0
                    phi = 0;
                else
                    phi = a_sorted(i) / (phi_total + a_sorted(i));
                end
                if abs(k_eff_low - k_i) < 1e-9
                    % unchanged
                else
                    k_eff_low = k_i + (1 - phi) / ((1/(k_eff_low - k_i)) + phi/(3*k_i));
                end
                phi_total = phi_total + a_sorted(i);
            end

            % upper bound
            [k_sorted, idx] = sort(k_present, 'descend');
            a_sorted = a_norm(idx);
            k_eff_high = k_sorted(1);
            phi_total = a_sorted(1);
            for i = 2:length(k_sorted)
                k_i = k_sorted(i);
                if (phi_total + a_sorted(i)) == 0
                    phi = 0;
                else
                    phi = a_sorted(i) / (phi_total + a_sorted(i));
                end
                if abs(k_eff_high - k_i) < 1e-9
                    % unchanged
                else
                    k_eff_high = k_i + (1 - phi) / ((1/(k_eff_high - k_i)) + phi/(3*k_i));
                end
                phi_total = phi_total + a_sorted(i);
            end
        end
        if material == 12
            k_pave = (k_eff_low + k_eff_high)/2;
        end

    elseif material == 14   % EMT
        C_pave = sum(C_material .* a);
        present_idx = a > 0;
        if ~any(present_idx)
            warning('EMT: No materials with area ratio > 0.');
            k_pave = 0;
        else
            k_present = k_material(present_idx);
            a_present = a(present_idx);
            a_norm = a_present / sum(a_present);
            k_min = min(k_present);
            k_max = max(k_present);
            if length(k_present) == 1
                k_pave = k_present(1);
            elseif abs(k_min - k_max) < 1e-9
                k_pave = k_min;
            else
                emt_func = @(k_e) sum(a_norm .* (k_present - k_e) ./ (k_present + 2 * k_e));
                try
                    k_pave = fzero(emt_func, [k_min, k_max]);
                catch ME
                    warning('EMT calculation failed: %s', ME.message);
                    k_pave = NaN;
                end
            end
        end

    elseif material == 15   % arithmic C and matching effusivities
        C_pave = sum(C_material .* a);
        eff_eff = sum(((k_material.*C_material).^0.5) .* a);
        k_pave = (eff_eff^2)/C_pave;

    elseif material == 16   % matching inverse diffusivities & arithmic effusivities
        diff_eff = (sum(((k_material./C_material).^(-1)) .* a)).^(-1);
        eff_eff  = sum(((k_material.*C_material).^0.5) .* a);
        C_pave = eff_eff / (diff_eff^0.5);
        k_pave = eff_eff * (diff_eff^0.5);

        elseif material == 17
        % =============== Salamanca et al. (2009), approach c ==============
        % Ref: Salamanca, Krayenhoff & Martilli (2009), J. Appl. Meteor.
        %      Climatol., 48, 1725-1732, doi:10.1175/2009JAMC2176.1
        P_day  = 86400;                                  
        C_pave = sum(a .* C_material);                    

        if exist('domain_depth','var')
            d_sal = domain_depth;                         % common column depth
        else
            % fallback: mean diurnally active depth of the constituents
            d_sal = mean(2.273*sqrt(43200*k_material./C_material));
        end

        Z_Dj = sqrt(k_material*P_day./(C_material*pi));   % patch damping depths

        num_sal = sum(a .* exp(d_sal./Z_Dj) .* sin(d_sal./Z_Dj));
        den_sal = sum(a .* exp(d_sal./Z_Dj) .* cos(d_sal./Z_Dj));

        theta = atan2(num_sal, den_sal);                  % phase, (-pi, pi]
        if theta <= 0, theta = theta + 2*pi; end          % positive principal phase

        present_sal = a > 0;
        ZD_lo = min(Z_Dj(present_sal));
        ZD_hi = max(Z_Dj(present_sal));

        n_branch   = 0:ceil(d_sal/(2*pi*ZD_lo)) + 1;      % enough branches to cover ZD_lo
        phase_cand = theta + 2*pi*n_branch;
        ZD_cand    = d_sal ./ phase_cand;

        tol      = 1e-9;
        in_range = ZD_cand >= ZD_lo - tol & ZD_cand <= ZD_hi + tol;
        ZD_target = sum(a .* Z_Dj);                       % tiebreak: area-weighted Z_D
        if any(in_range)
            cand_ok = ZD_cand(in_range);
            [~, ibest] = min(abs(cand_ok - ZD_target));
            Z_D_sal = cand_ok(ibest);
        else
            [~, ibest] = min(abs(ZD_cand - ZD_target));
            Z_D_sal = ZD_cand(ibest);
            warning(['Salamanca 2009: no arctan branch places Z_D within ' ...
                     '[min,max] of constituent damping depths; using ' ...
                     'closest candidate Z_D = %.4f m.'], Z_D_sal);
        end

        k_pave = C_pave * pi * Z_D_sal^2 / P_day;         % Eq. (5c): k = C*pi*Z_D^2/P
        fprintf(['Salamanca 2009 approach c: d = %.3f m | Z_D = %.4f m | ' ...
                 'k_eq = %.4f W/m/K | C_eq = %.3e J/m3/K\n'], ...
                 d_sal, Z_D_sal, k_pave, C_pave);
    end

    kapa_pave = k_pave / C_pave;
    alpha = kapa_pave;
    if exist('domain_depth','var')
        F_Tg_depth = domain_depth;
    else
        F_Tg_depth = (2.222*((43200*kapa_pave)^(0.5)));
    end
    if material > 6
        albedos = (a(1)*albedo_material(1)) + (a(2)*albedo_material(2)) + ...
                  (a(3)*albedo_material(3)) + (a(4)*albedo_material(4)) + ...
                  (a(5)*albedo_material(5)) + (a(6)*albedo_material(6));
    end

    % ---- ground properties & discretization -----------------
    k = k_pave;            % thermal conductivity (W/mK)
    c = C_pave;            % volumetric heat capacity (J/m3/K)
    alpha = k / c;         % thermal diffusivity (m^2/s)

    total_depth = 1*ceil(F_Tg_depth*100)/100;   % total model depth (m)
    dz = matrix_dx(aaaa);                        % layer thickness (m)
    num_layers = floor(total_depth/dz);

    time_steps  = n_model;                       % dynamic, from data span & model dt
    time_vector = time;                          % seconds

    Fo = alpha * dt / dz^2;                       % Fourier number
    if Fo > 0.5
        fprintf('Warning: stability criterion Fo = %.3f > 0.5 (may be unstable).\n', Fo);
    end

    % ---- initialize state -----------------------------------------------
    T_profile = zeros(num_layers + 1, time_steps);
    Ts_open   = zeros(time_steps, 1);
    G_flux    = zeros(time_steps, 1);
    initial_T = mean(Ta);
    T_profile(:, 1) = initial_T;
    Ts_open(1,1)    = initial_T;

    % =================== SPIN-UP + EXPLICIT FD SOLVER =========
    max_runs = 100;
    convergence_threshold = 0.1;
    is_converged = false;
    run_count = 0;
    T_profile_previous_run = T_profile;
    while ~is_converged && run_count < max_runs
        run_count = run_count + 1;
        fprintf('Starting simulation run #%d...\n', run_count);
        T_profile(num_layers + 1, :) = mean(Ts_open);
        for i = 2:time_steps
            G_func = @(T_surf) (k / dz) * (T_surf - T_profile(2, i-1));
            y = @(xx) ( ((1-albedos)*SW_all_g(i)) + (LWin(i)) ...
                       - (es*sigma*((xx+273.15).^4)) ...
                       - (h_s(i)*(xx-Ta(i))) - G_func(xx) );
            xx0 = Ts_open(i-1,1);
            Ts_open(i,1) = fzero(y, xx0);
            T_profile(1, i) = Ts_open(i, 1);
            for j = 2:num_layers
                T_profile(j, i) = T_profile(j, i-1) + ...
                    Fo * (T_profile(j+1, i-1) - 2*T_profile(j, i-1) + T_profile(j-1, i-1));
            end
            G_flux(i,1) = (k / dz) * (T_profile(1, i) - T_profile(2, i));
        end
        max_temp_change = max(abs(T_profile(:, end) - T_profile_previous_run(:, end)));
        fprintf('Maximum temperature change between runs: %.4f C\n', max_temp_change);
        if max_temp_change < convergence_threshold
            is_converged = true;
            fprintf('Solution has converged!\n');
        else
            T_profile_previous_run = T_profile;
            T_profile(:, 1) = T_profile(:, end);
            Ts_open(1,1) = T_profile(1,1);
        end
    end
    if ~is_converged
        fprintf('Warning: Solution did not converge after %d runs.\n', max_runs);
    end

    % =================== SAVE -> ./Results/<date>/<scenario>.mat =========
    F_Tg_depth = total_depth;

    results_root = fullfile(script_dir, 'Results');
    date_folder  = fullfile(results_root, datestr(now, 'yyyy-mm-dd'));
    if ~exist(date_folder, 'dir'), mkdir(date_folder); end

    mat_label = material_name{material};
    
    if material > 6
        ratio_tag = ['_a' strjoin(compose('%02d', round(a*100)), '-')];  % e.g. _a10-20-30-10-20-10
    else
        ratio_tag = '';
    end
    scenario = sprintf('m%02d%s_d%.2f_dz%.3f', material, ratio_tag, total_depth, dz);
    out_path = fullfile(date_folder, [scenario '.mat']);

    meta = struct('material', material, 'material_name', mat_label, ...
        'a', a, 'domain_depth', domain_depth, 'total_depth', total_depth, ...
        'dz', dz, 'num_layers', num_layers, 'model_dt_seconds', dt, ...
        'data_dt_seconds', dt_data_env, 'n_model', n_model, ...
        'k_pave', k_pave, 'C_pave', C_pave, 'albedos', albedos);

    save(out_path, 'T_profile', 'G_flux', 'F_Tg_depth', ...
        'time_vector', 'meta', '-v7.3');
    fprintf('Saved: %s\n', out_path);

    % =================== optional figures =====================
    if plots == 1
        figure('Name','Surface Conditions Over Time (Coupled Model)');
        yyaxis left;
        plot(time_vector/3600, Ts_open, 'b-', 'LineWidth', 1.5);
        ylabel('Surface Temperature (\circC)'); ax = gca; ax.YColor = 'b';
        yyaxis right;
        plot(time_vector/3600, G_flux, 'r--', 'LineWidth', 1.5);
        ylabel('Ground Heat Flux (W/m^2)'); ax.YColor = 'r';
        xlabel('Time (hours)');
        title(sprintf('Surface T and ground heat flux: %s, depth %.0f cm', ...
            mat_label, total_depth*100));
        legend('Surface Temp','Ground Heat Flux','Location','northwest'); grid on;

        figure('Name','Ground Temperature Profiles (Coupled Model)');
        depth_vector = (0:num_layers) * dz;
        plot_indices = round(linspace(1, time_steps, 7));
        hold on; color_map = lines(length(plot_indices));
        for i = 1:length(plot_indices)
            idx = plot_indices(i);
            plot(T_profile(:, idx), -depth_vector, 'LineWidth', 1.5, ...
                'Color', color_map(i,:), ...
                'DisplayName', sprintf('time = %.1f', time_vector(idx)/3600));
        end
        hold off; xlabel('Temperature (\circC)'); ylabel('Depth (m)');
        title(sprintf('Ground T profiles: %s, depth %.0f cm', mat_label, total_depth*100));
        legend('show','Location','southeast'); grid on;

        figure('Name','Ground Temperature Heatmap');
        imagesc(time_vector/3600, -depth_vector, T_profile);
        set(gca,'YDir','normal'); colorbar;
        xlabel('Time (hours)'); ylabel('Depth (m)');
        title(sprintf('Ground T (depth vs time): %s, depth %.0f cm', mat_label, total_depth*100));
        clim([min(T_profile(:)), max(T_profile(:))]);
        ylabel(colorbar,'Temperature (\circC)');
    end

end
end
end
end

%% =================== local functions ====================================
function SW_all_g = compute_SW(dir_ground_native, diff_ground_native, n_model)
idx = linspace(1, numel(dir_ground_native), n_model).';
Dirrg = interp1(dir_ground_native,  idx, 'linear');
Diffg = interp1(diff_ground_native, idx, 'linear');

sunrise = find(diff(Dirrg > 0) ==  1) + 1;
sunset  = find(diff(Dirrg > 0) == -1) + 1;
z = zeros(size(Dirrg));
if ~isempty(sunrise) && ~isempty(sunset)
    z(sunrise:sunset) = linspace(180, 0, sunset - sunrise + 1);
end
SW_all_g = Dirrg + Diffg;                     % horizontal irradiance over ground
end

function s = sanitize_name(name)
s = char(name);
s = regexprep(s, '[^\w\.\-]+', '_');   
s = regexprep(s, '_+', '_');            
s = regexprep(s, '^_|_$', '');         
if isempty(s), s = 'scenario'; end
end