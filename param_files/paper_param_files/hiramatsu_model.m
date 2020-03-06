% hiramatsu model
% 
% Two-body Hiramatsu model (to show we don't need to go past rv)
%
% This model is designed to replace the many earlier models for the paper. The
% relevant studies - varying boundary size and mesh fineness - are controlled
% by simplified variables.
% 
% Parameters which must be set:
% 
% log2au_bound      | Integer from 1 to 10. Log distance to computational
%                   |  boundary. Boundary set at 2^log2au_bound.
% fineness_level    | Integer from 1 to 6, higher means finer mesh

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-09-11
% Updated: 2020-01-06

% Iteration parameters
max_iter = 500;
min_nu_limit = 1E-10;
solution_method_switch_limit = 10;

% Physical parameters

k = 5.934424261E-4; % k = 3 / [ sqrt(8 pi G rho_c) r_c ]
body_radii = [0.3; 0.1];
body_densities = [1; 0.3375];
body_z_coordinates = [1; 0];

zA = body_z_coordinates(1);
zB = body_z_coordinates(2);
RsA = body_radii(1);
RsB = body_radii(2);
rhoA = body_densities(1);
rhoB = body_densities(2);


% Initial condition
initial_condition = 'sum_onebodies';


assert(log2au_bound == floor(log2au_bound))
assert((log2au_bound <= 10) && (log2au_bound >= 1))

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

z_points = [zA + RsA, zA - RsA, zB + RsB, zB - RsB];
z_center = (max(z_points) + min(z_points))/2;

meshes = {};
for j=0:0
    this_r_mesh = base_r_finemesh*2^(j);
    this_z_mesh = base_z_finemesh*2^(j) + z_center;
    meshes{end+1} = {this_r_mesh, this_z_mesh};
end
for j=1:log2au_bound
    this_r_mesh = base_r_mesh*2^(j);
    this_z_mesh = base_z_mesh*2^(j) + z_center;
    meshes{end+1} = {this_r_mesh, this_z_mesh};
end

% Boundary condition
boundary_condition = 'combined_onebody';

% Output parameters
data_dir = 'vfd_data';
output_filename = sprintf('paper_hiramatsu_log2aubound%02d_meshlevel%02d',log2au_bound, fineness_level);
clobber = true; % true => overwrite

plot_while_solving = true;
plot_update_steps = 1;
plots_to_show = {'residual_history'};
plots_to_save = plots_to_show;
