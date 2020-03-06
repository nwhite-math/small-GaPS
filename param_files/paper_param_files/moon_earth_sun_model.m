% moon_earth_sun_model
% 
% Three-body Sun-Earth model
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
% sun_only          | Boolean. If true, simulate only the Sun
% earth_only        | Boolean. If true, simulate only Earth
% moon_only         | Boolean. If true, simulate only the Moon
% earth_moon_only   | Boolean. If true, simulate only Earth and the Moon
% sun_earth_only    | Boolean. If true, simulate only the Sun and Earth
% earth_moon_boxes  | Boolean. If false, simulations enclose Earth, the Moon
%                   |  and the Sun, in which case log2au_bound must be between
%                   |  1 and 9. If true, simulations enclose Earth and the Moon
%                   |  with the single-body Sun solution as boundary condition,
%                   |  in which case log2au_bound must be between -14 and -1.

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-08-22
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

sun_r_meters = 695.7E6;

k = 0;                  % k = 3 / [ sqrt(8 pi G rho_sun) r_c ]
RSun = 1.0;             % Sun radius / Sun radius
rhoSun = 1;             % Sun density / Sun density
xSun = 0; ySun = 0;
zSun = AU/2;            % Sun z location / Sun radius
REarth = 0.00915768;    % Earth radius / Sun radius
rhoEarth = 3.9169;      % Earth density / Sun density
xEarth = 0; yEarth = 0;
zEarth = -AU/2;         % Earth z location / Sun radius

% Lunar inclination (from ecliptic plane)
moon_incl = 5.14 * pi/180;
moon_orbit_r = 385E6 / sun_r_meters; % ( approximating circular for now )


rotzx_matrix = @(theta) [ cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
rotxy_matrix = @(theta) [ cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
rotyz_matrix = @(theta) [ 1, 0, 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)];

%% Full and over-complicated moon position
%cur_moon_theta = 0 * pi/180; % (0 = away from sun; 180 * pi / 180 = pi = close to sun)
%incl_matrix = rotyz_matrix(-moon_incl);
%moon_angle_matrix = rotzx_matrix(cur_moon_theta);
%moon_xyz_rel = incl_matrix*moon_angle_matrix*[0; 0; -1]*moon_orbit_r;

% Direct angle with respect to sun
cur_moon_theta = -90 * pi/180; % (0 = away from sun; 180 * pi / 180 = pi = close to sun)
moon_angle_matrix = rotyz_matrix(cur_moon_theta);
moon_xyz_rel = moon_angle_matrix*[0; 0; -1]*moon_orbit_r;

RMoon = 2.497E-3;  % Moon radius / Sun radius 
rhoMoon = 2.375;   % Moon density / Sun density
xMoon = xEarth + moon_xyz_rel(1); % Moon x location / Sun radius
yMoon = yEarth + moon_xyz_rel(2); % Moon y location / Sun radius
zMoon = zEarth + moon_xyz_rel(3); % Moon z location / Sun radius

body_radii = [RSun; REarth; RMoon];
body_densities = [rhoSun; rhoEarth; rhoMoon];
body_coordinates = [ xSun   ySun   zSun;...
                     xEarth yEarth zEarth;...
                     xMoon  yMoon  zMoon];

if (sun_only)
    body_str = 'sunonly';
    assert(earth_moon_boxes == false)
    assert(earth_moon_only == false)
    assert(sun_earth_only == false)
    assert(earth_only == false)
    assert(moon_only == false)
    body_radii = [RSun];
    body_densities = [rhoSun];
    body_coordinates = [xSun ySun zSun];
elseif (earth_only)
    body_str = 'earthonly';
    assert(earth_moon_boxes == false)
    assert(earth_moon_only == false)
    assert(sun_earth_only == false)
    assert(moon_only == false)
    body_radii = [REarth];
    body_densities = [rhoEarth];
    body_coordinates = [xEarth yEarth zEarth];
elseif (moon_only)
    body_str = 'moononly';
    assert(earth_moon_boxes == false)
    assert(earth_moon_only == false)
    assert(sun_earth_only == false)
    body_radii = [RMoon];
    body_densities = [rhoMoon];
    body_coordinates = [xMoon yMoon zMoon];
elseif sun_earth_only
    body_str = 'sunearthonly';
    assert(earth_moon_boxes == false)
    assert(earth_moon_only == false)
    body_radii = [RSun; REarth];
    body_densities = [rhoSun; rhoEarth];
    body_coordinates = [ xSun   ySun   zSun;...
                         xEarth yEarth zEarth];
elseif earth_moon_only
    body_str = 'earthmoononly';
    assert(earth_moon_boxes == false)
    body_radii = [REarth; RMoon];
    body_densities = [rhoEarth; rhoMoon];
    body_coordinates = [ xEarth   yEarth   zEarth;...
                         xMoon yMoon zMoon];
elseif earth_moon_boxes
    body_str = 'earthmoonbox';
    body_radii = [REarth; RMoon];
    body_densities = [rhoEarth; rhoMoon];
    body_coordinates = [ xEarth   yEarth   zEarth;...
                         xMoon yMoon zMoon];
else
    body_str = 'sunearthmoon';
    body_radii = [RSun; REarth; RMoon];
    body_densities = [rhoSun; rhoEarth; rhoMoon];
    body_coordinates = [ xSun   ySun   zSun;...
                         xEarth yEarth zEarth;...
                         xMoon  yMoon  zMoon];
end


% Initial condition
% ic = sum_onebodies OR zero OR random
ic_str = ic;
if strcmpi(ic, 'sum_onebodies')
    initial_condition = 'sum_onebodies';
elseif strcmpi(ic, 'zero')
    initial_condition = 'custom';
    ic_func = @(x, y, z) 0*x;
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
        % So, we choose kx, ky, and kz uniformly at random on a log scale from 10^(-3) to 10^(3)
        kx = 10.^(2*(rand(1,num_sins)-0.5)*3);
        ky = 10.^(2*(rand(1,num_sins)-0.5)*3);
        kz = 10.^(2*(rand(1,num_sins)-0.5)*3);
        rand_theta_x = rand(1,num_sins)*2*pi;
        rand_theta_y = rand(1,num_sins)*2*pi;
        rand_theta_z = rand(1,num_sins)*2*pi;
        rand_coeffs_a = randn(num_sins,1) .* (( sqrt(kx.^2+ky.^2+kz.^2) )'.^noise_power);
        rand_coeffs_a = rand_coeffs_a*2*phi_diff_max/sqrt(sum(rand_coeffs_a.^2));
        ic_func = @(x, y, z) (sin(x*kx + (1+0*x)*rand_theta_x).*sin(y*ky + (1+0*y)*rand_theta_y).*sin(z*kz + (1+0*z)*rand_theta_z))*rand_coeffs_a;
    elseif strcmpi(noise, 'white')
        ic_func = @(x, y, z) phi_diff_max*randn(size(z));
    else
        assert(false)
    end

    % NOTE: This is ok since ic_func is called only a single time
    ic_str = sprintf('%snoise_randseed%d', noise,seed);
else
    assert(false);
end


assert(log2au_bound == floor(log2au_bound))
if (~earth_moon_boxes)
    assert((log2au_bound <= 9) && (log2au_bound >= 1))
else
    assert((log2au_bound <= -1))
    %assert(log2au_bound >= -floor(log(AU/REarth)/log(2)));
    assert(log2au_bound >= -floor(log(AU/moon_orbit_r)/log(2)));
end

% mesh settings
%box_width = (2^log2au_bound)*AU;

alignment_tol = 1E-8;
if (fineness_level == 1)
    n = 4;
elseif (fineness_level == 2) % "standard (?)"
    n = 6;
elseif (fineness_level == 3)
    n = 8;
elseif (fineness_level == 4)
    n = 10;
elseif (fineness_level == 5)
    n = 12;
elseif (fineness_level == 6)
    n = 16;
else
    assert(false)
end

side_mesh_points = n;
assert(side_mesh_points/2==floor(side_mesh_points/2));
base_x_mesh = linspace(-1,1,2*side_mesh_points+1);
base_y_mesh = linspace(-1,1,2*side_mesh_points+1);
base_z_mesh = linspace(-1,1,2*side_mesh_points+1);
base_x_finemesh = linspace(-1,1,4*side_mesh_points+1);
base_y_finemesh = linspace(-1,1,4*side_mesh_points+1);
base_z_finemesh = linspace(-1,1,4*side_mesh_points+1);

meshes = {};
% Start with wide meshes enclosing the Earth, Moon, and Sun
for j=1:log2au_bound
    this_x_mesh = base_x_mesh*AU*2^(j);
    this_y_mesh = base_y_mesh*AU*2^(j);
    this_z_mesh = base_z_mesh*AU*2^(j);
    meshes{end+1} = {this_x_mesh, this_y_mesh, this_z_mesh};
end
% The last combined mesh has extra fineness
this_x_mesh = base_x_finemesh*AU;
this_y_mesh = base_y_finemesh*AU;
this_z_mesh = base_z_finemesh*AU;
meshes{end+1} = {this_x_mesh, this_y_mesh, this_z_mesh};

% Then, make nested meshes around the Earth and Sun individually


% SUN BOXES
sun_bound_limit = floor(log(AU/RSun)/log(2));
for j=2:sun_bound_limit
    this_x_mesh = base_x_mesh*AU*2^(-j);
    this_y_mesh = base_y_mesh*AU*2^(-j);
    this_z_mesh = base_z_mesh*AU*2^(-j) + zSun;
    meshes{end+1} = {this_x_mesh, this_y_mesh, this_z_mesh};
end

if (earth_moon_boxes)
    meshes = {};
    earth_bound_upper_start = -log2au_bound;
else
    earth_bound_upper_start = 2;
end

earth_bound_limit = floor(log(AU/REarth)/log(2));
moon_bound_limit = floor(log(AU/RMoon)/log(2));
moon_orbit_limit = floor(log(AU/moon_orbit_r)/log(2));

% EARTH + MOON BOXES
for j=earth_bound_upper_start:moon_orbit_limit-1
    this_x_mesh = base_x_mesh*AU*2^(-j);
    this_y_mesh = base_y_mesh*AU*2^(-j);
    this_z_mesh = base_z_mesh*AU*2^(-j) + zEarth;
    meshes{end+1} = {this_x_mesh, this_y_mesh, this_z_mesh};
end

% FINAL EARTH + MOON BOX (extra fineness)
for j=moon_orbit_limit:moon_orbit_limit
    this_x_mesh = base_x_finemesh*AU*2^(-j);
    this_y_mesh = base_y_finemesh*AU*2^(-j);
    this_z_mesh = base_z_finemesh*AU*2^(-j) + zEarth;
    meshes{end+1} = {this_x_mesh, this_y_mesh, this_z_mesh};
end
final_earth_moon_mesh = meshes{end};

%EARTH BOXES
for j=moon_orbit_limit+2:earth_bound_limit
    this_x_mesh = base_x_mesh*AU*2^(-j);
    this_y_mesh = base_y_mesh*AU*2^(-j);
    this_z_mesh = base_z_mesh*AU*2^(-j) + zEarth;
    meshes{end+1} = {this_x_mesh, this_y_mesh, this_z_mesh};
end

%% MOON BOXES
%moon_bound_limit = moon_orbit_limit +2; % DEBUG
[~, xMoonApprox] = nearest_index(final_earth_moon_mesh{1}, xMoon);
[~, yMoonApprox] = nearest_index(final_earth_moon_mesh{2}, yMoon);
[~, zMoonApprox] = nearest_index(final_earth_moon_mesh{3}, zMoon);
for j=moon_orbit_limit+2:moon_bound_limit
    this_x_mesh = base_x_mesh*AU*2^(-j) + xMoonApprox;
    this_y_mesh = base_y_mesh*AU*2^(-j) + yMoonApprox;
    this_z_mesh = base_z_mesh*AU*2^(-j) + zMoonApprox;
    meshes{end+1} = {this_x_mesh, this_y_mesh, this_z_mesh};
    this_mesh = meshes{end};
    [~, xMoonApprox] = nearest_index(this_mesh{1}, xMoon);
    [~, yMoonApprox] = nearest_index(this_mesh{2}, yMoon);
    [~, zMoonApprox] = nearest_index(this_mesh{3}, zMoon);
end


% Boundary condition
boundary_condition = 'combined_onebody';
if earth_moon_boxes
    boundary_condition = 'custom';
    bc_func = @(x,y,z) vainshteinone(0, rhoSun, RSun, sqrt((x-xSun).^2+(y-ySun).^2+(z-zSun).^2));
end


% Output parameters
data_dir = 'vfd_data';
output_filename = sprintf('paper_c390_%s_log2aubound%02d_meshlevel%02d_ic_%s',body_str, log2au_bound, fineness_level, ic_str);
clobber = true; % true => overwrite

% Import mesh if possible
similar_dirs = {};
try
    caller_m = dbstack(2, '-completenames');
    caller_m = caller_m.file;
    caller_dir = fileparts(caller_m);
    full_data_dir = fullfile(caller_dir, data_dir);
    similar_dirs = dir(sprintf('%saper_c390_*_log2aubound%02d_meshlevel%02d_ic_*', fullfile(full_data_dir,'p'), log2au_bound, fineness_level));
end
if length(similar_dirs) > 0
    % there is another file with the same mesh
    import_mesh = true;
    similar_dir = similar_dirs(1).name;
    import_mesh_file = fullfile(data_dir, similar_dir, [similar_dir '.mat']);
else
    import_mesh = false;
end



% plot while solving
plot_while_solving = true;
plot_update_steps = 1;
plots_to_show = {'residual_history'};
plots_to_save ={'residual_history'};
