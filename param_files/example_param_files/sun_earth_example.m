% sun_earth_example.m
% 
% Compute 2-body solution for Sun and Earth system.
% 

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-11-03
% Updated: 2020-01-06

% Iteration parameters
max_iter = 500; % Actually, the simulation converges and exits long before this.
min_nu_limit = 1E-10;
solution_method_switch_limit = 4;

% We nondimensionalize by the Sun's radius and density
k = 0;                  % Really k is around 1E-15, but that is negligibly small

AU = 215.032;           % 1 AU / Sun radius
RSun = 1.0;             % Sun radius / Sun radius
rhoSun = 1;             % Sun density / Sun density
zSun = AU/2;            % Sun z location / Sun radius
REarth = 0.00915768;    % Earth radius / Sun radius
rhoEarth = 3.9169;      % Earth density / Sun density
zEarth = -AU/2;         % Earth z location / Sun radius

body_radii = [RSun; REarth];
body_densities = [rhoSun; rhoEarth];
body_z_coordinates = [zSun; zEarth];

% Initial condition
initial_condition = 'sum_onebodies';

% Mesh settings
% Let the computational domain be a cylinder with radius 2^x AU and z-extent
% 2^(x+1) AU, where x here is denoted by the variable log2au_bound.
log2au_bound = 6;

% Let nested meshes have (n+1) x (2n +1) points. In this example, we make
% subsequent meshes half the extent (one quarter the area) of prior meshes.
% Thus, n should be divible by 2.
n = 32; % => 33 x 65 point meshes
assert(n/2==floor(n/2));

base_r_mesh = linspace(0,1,n+1);
base_z_mesh = linspace(-1,1,2*n+1);
base_r_finemesh = linspace(0,1,2*n+1);
base_z_finemesh = linspace(-1,1,4*n+1);

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

% Make nested meshes around the Sun and Earth individually. We keep nesting
% smaller and smaller meshes until we get down to the size of each body.

sun_bound_limit = floor(log(AU/RSun)/log(2)); % Smallest mesh around Sun
for j=2:sun_bound_limit
    this_r_mesh = base_r_mesh*AU*2^(-j);
    this_z_mesh = base_z_mesh*AU*2^(-j) + zSun;
    meshes{end+1} = {this_r_mesh, this_z_mesh};
end

earth_bound_limit = floor(log(AU/REarth)/log(2)); % Smallest mesh around Earth
for j=2:earth_bound_limit
    this_r_mesh = base_r_mesh*AU*2^(-j);
    this_z_mesh = base_z_mesh*AU*2^(-j) + zEarth;
    meshes{end+1} = {this_r_mesh, this_z_mesh};
end

% Boundary condition
boundary_condition = 'combined_onebody';

% Output parameters
data_dir = 'example_output';
output_filename = 'sun_earth_example';
clobber = true; % true => overwrite

% Show the integrated residual while the simulation is running
plot_while_solving = true;
plot_update_steps = 1;
plots_to_show = {'residual_history'};
plots_to_save = plots_to_show;
