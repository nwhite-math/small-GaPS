% earth_sun_model
% 
% Two-body Sun-Earth model in 2D cylindrical coordinates
%
% This model is designed to replace the many earlier models for the paper. All
% the relevant studies - varying initial conditions, boundary size, and mesh
% fineness - are controlled by simplified variables.
% 
% Parameters which must be set:
% 
% log2au_bound      | Integer from 1 to 9 or -14 to -1. Log distance to
%                   |  computational boundary. Boundary set at 2^log2au_bound AU
% fineness_level    | Integer from 1 to 6, higher means finer mesh
%                   |
% ic                | Initial condition. Options same as vain_findiff:
%                   |  'sum_onebodies', 'combined_onebody', 'newton', etc.
%                   |
% sun_only          | Boolean. If true, simulate only the Sun without Earth
% earth_only        | Boolean. If true, simulate only Earth without the Sun
% earth_boxes       | Boolean. If false, perform simulations enclosing Earth
%                   |  and the Sun, in which case log2au_bound must be between
%                   |  1 and 9. If true, perform simulations around Earth with
%                   |  the single-body Sun solution as boundary condition, in
%                   |  which case log2au_bound must be between -14 and -1.

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-08-18
% Updated: 2020-01-06

% Iteration parameters
max_iter = 500;
min_nu_limit = 1E-10;
solution_method_switch_limit = 10;

% Physical parameters
% We nondimensionalize by the Sun radius (6.96 * 10^8 m) and density (1410 kg/m^3),
% and take r_c = c/H_0 = 1.499*10^26.

% The (1 AU) / (sun radius) = 215 (approximately)
AU = 215.032;

k = 0;                  % k = 3 / [ sqrt(8 pi G rho_sun) r_c ]
RSun = 1.0;             % Sun radius / Sun radius
rhoSun = 1;             % Sun density / Sun density
zSun = AU/2;            % Sun z location / Sun radius
REarth = 0.00915768;    % Earth radius / Sun radius
rhoEarth = 3.9169;      % Earth density / Sun density
zEarth = -AU/2;         % Earth z location / Sun radius

if (sun_only)
    body_str = 'sunonly';
    assert(earth_boxes == false)
    assert(earth_only == false)
    body_radii = [RSun];
    body_densities = [rhoSun];
    body_z_coordinates = [zSun];
elseif earth_only
    body_str = 'earthonly';
    assert(earth_boxes == false)
    body_radii = [REarth];
    body_densities = [rhoEarth];
    body_z_coordinates = [zEarth];
elseif earth_boxes
    body_str = 'earthbox';
    body_radii = [REarth];
    body_densities = [rhoEarth];
    body_z_coordinates = [zEarth];
else
    body_str = 'sunearth';
    body_radii = [RSun; REarth];
    body_densities = [rhoSun; rhoEarth];
    body_z_coordinates = [zSun; zEarth];
end


% Initial condition
% ic = sum_onebodies OR zero OR random
ic_str = ic;
if strcmpi(ic, 'sum_onebodies')
    initial_condition = 'sum_onebodies';
elseif strcmpi(ic, 'zero')
    initial_condition = 'custom';
    ic_func = @(r, z) 0*r;
elseif strcmpi(ic, 'random')
    % 270.2874
    phi_diff_max = abs(diff(vainshteinone(k, rhoSun, RSun, [0 AU*2^9])));

    assert(exist('seed','var')==1);
    assert(exist('noise','var')==1);

    % Initial condition
    initial_condition = 'custom';
    rng(seed);

    if strcmpi(noise, 'red')
        num_sins = 100; % add up 100 distinct sin functions
        noise_power = -1; % red noise
        % we want the random length scale to vary from about 10^(-3) (sub Earth radius) to 10^3 (about 4 AU)
        % So, we choose kr and kz uniformly at random on a log scale from 10^(-3) to 10^(3)
        kz = 10.^(2*(rand(1,num_sins)-0.5)*3);
        kr = 10.^(2*(rand(1,num_sins)-0.5)*3);
        rand_theta_z = rand(1,num_sins)*2*pi;
        rand_theta_r = rand(1,num_sins)*2*pi;
        rand_coeffs_a = randn(num_sins,1).* (( sqrt(kr.^2+kz.^2) )'.^noise_power);
        rand_coeffs_a = rand_coeffs_a*2*phi_diff_max/sqrt(sum(rand_coeffs_a.^2));
        ic_func = @(r, z) (sin(z*kz + (1+0*z)*rand_theta_z).*sin(r*kr + (1+0*r)*rand_theta_r))*rand_coeffs_a;
    elseif strcmpi(noise, 'white')
        ic_func = @(r, z) phi_diff_max*randn(size(z));
    else
        assert(false)
    end

    % NOTE: This is ok since ic_func is called only a single time
    ic_str = sprintf('%snoise_randseed%d', noise,seed);
else
    assert(false);
end


assert(log2au_bound == floor(log2au_bound))
if (~earth_boxes)
    assert((log2au_bound <= 9) && (log2au_bound >= 1))
else
    assert((log2au_bound <= -1))
    assert(log2au_bound >= -floor(log(AU/REarth)/log(2)));
end

% mesh settings
%box_width = (2^log2au_bound)*AU;

if (fineness_level == 1)
    n = 4;
elseif (fineness_level == 2)
    n = 8;
elseif (fineness_level == 3)
    n = 16;
elseif (fineness_level == 4) % "standard"
    n = 32;
elseif (fineness_level == 5)
    n = 64;
elseif (fineness_level == 6)
    n = 128;
else
    assert(false)
end

side_mesh_points = n;
assert(side_mesh_points/2==floor(side_mesh_points/2));
base_r_mesh = linspace(0,1,side_mesh_points+1);
base_z_mesh = linspace(-1,1,2*side_mesh_points+1);
base_r_finemesh = linspace(0,1,2*side_mesh_points+1);
base_z_finemesh = linspace(-1,1,4*side_mesh_points+1);

meshes = {};
% Start with wide meshes enclosing both the Earth and Sun
for j=1:log2au_bound
    this_r_mesh = base_r_mesh*AU*2^(j);
    this_z_mesh = base_z_mesh*AU*2^(j);
    meshes{end+1} = {this_r_mesh, this_z_mesh};
end
% The last combined mesh has extra fineness
this_r_mesh = base_r_finemesh*AU;
this_z_mesh = base_z_finemesh*AU;
meshes{end+1} = {this_r_mesh, this_z_mesh};

% Then, make nested meshes around the Earth and Sun individually
sun_bound_limit = floor(log(AU/RSun)/log(2));
for j=2:sun_bound_limit
    this_r_mesh = base_r_mesh*AU*2^(-j);
    this_z_mesh = base_z_mesh*AU*2^(-j) + zSun;
    meshes{end+1} = {this_r_mesh, this_z_mesh};
end

if (earth_boxes)
    meshes = {};
    earth_bound_upper_start = -log2au_bound;
else
    earth_bound_upper_start = 2;
end
earth_bound_limit = floor(log(AU/REarth)/log(2));
for j=earth_bound_upper_start:earth_bound_limit
    this_r_mesh = base_r_mesh*AU*2^(-j);
    this_z_mesh = base_z_mesh*AU*2^(-j) + zEarth;
    meshes{end+1} = {this_r_mesh, this_z_mesh};
end

% Boundary condition
boundary_condition = 'combined_onebody';
if earth_boxes
    boundary_condition = 'custom';
    bc_func = @(r,z) vainshteinone(0, rhoSun, RSun, sqrt(r.^2+(z-zSun).^2));
end

% Output parameters
data_dir = 'vfd_data';
output_filename = sprintf('paper_c_%s_log2aubound%02d_meshlevel%02d_ic_%s',body_str, log2au_bound, fineness_level, ic_str);
clobber = true; % true => overwrite

% Import mesh if possible
caller_m = dbstack(2, '-completenames');
caller_m = caller_m.file;
caller_dir = fileparts(caller_m);
full_data_dir = fullfile(caller_dir, data_dir);
similar_dirs = dir(sprintf('%saper_c_*_log2aubound%02d_meshlevel%02d_ic_*', fullfile(full_data_dir,'p'), log2au_bound, fineness_level));
if length(similar_dirs) > 0
    % there is another file with the same mesh
    import_mesh = true;
    similar_dir = similar_dirs(1).name;
    import_mesh_file = fullfile(data_dir, similar_dir, [similar_dir '.mat']);
else
    import_mesh = false;
end

plot_while_solving = true;
plot_update_steps = 1;
plots_to_show = {'residual_history'};
plots_to_save = plots_to_show;
