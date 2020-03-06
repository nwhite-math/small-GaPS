% vain_findiff_3d.m
%
% A finite difference method to solve the cubic Galileon gravity (CGG) potential
% equation in axially-symmetric cylindrical coordinates. IMPORTANT: the user is
% responsible for nondimensionalizing the system.
%
% OVERVIEW OF ANALYTICAL SYSTEM AND NONDIMENSIONALIZATION
%
% The CGG potential equation is:
%    0 = k lap phi + [  (lap phi)^2 - d_ij phi d^ij phi ] - rho
%
% The user is responsible for nondimensionalizing the problem so that it can be
% described in this form. In particular, the user must choose a characteristic
% length scale d and a characteristic density rho_c. Then, all coordinates must
% be expressed in units of d, and densities in units of rho_0. The value of k
% can then be determined as:
%    k = 3 / [ sqrt(8 pi G rho_0) r_c ],
% where G is Newton's gravitational constant and r_c is the Vainshtein crossover
% scale. These choices automatically determine a characteristic potential
%    phi_0 = d^2 sqrt(8 pi G rho_0) / r_c
%
% Thus, the vain_findiff_3d solver does not have to know rho_0, d, or r_c. As
% long as the user inputs coordinates and body sizes in terms of d, densities in
% terms of rho_0, and the proper k, the result can be computed. Then, to
% interpret the result as a physical quantity, the user must multiply distances
% by d and potentials by phi_c to get back into physical units.
%
% OVERVIEW OF NUMERICAL SCHEME
%
% We define the residual as:
%    Res(phi) = lap phi + k/2 - sqrt( k^2/4 + rho + d_ij phi d^ij phi )
% so that for a true solution, Res(phi) = 0.
%
% Linearizing around a given function phi, we get the operator
%    Lop = lap - [(dij phi)/sqrt( k^2/4 + rho + d_ij phi d^ij phi )] d_ij
%
% After iteration phi_n, an update is computed by:
%    Lop gamma + Res(phi_n) = 0
%    phi_{n+1} = phi_n + nu*gamma
%
%
% USAGE
%
% WARNING: To avoid erroneous input, vain_findiff_3d clears all functions and
% all relevant input variables from memory. Be sure to save your work before
% running it. Also, in case there is a bug in the program, please clear the
% workspace (clearvars; clear functions) before running if possible.
%
% First, create a parameter file (e.g., my_params.m) defining the all relevant
% system and output parameters (see list below for options).
% Then, run the following:
%
% > param_file = 'my_params.m';
% > save_output = true; % (or false, for a trial run)
% > vain_findiff_3d
%
% Note: the save_output variable may be defined inside the parameter file.
%
%
% OVERVIEW OF INPUT PARAMETERS
%
% Physical Parameters
%
% k                  -> 3 / [ sqrt(8 pi G rho_0) r_c ]
%
% body_radii         -> List specifying the radius of each spherical body.
% body_masses        -> List specifying the mass of each spherical body.
% body_densities     -> List specifying the density of each spherical body.
% body_coordinates   -> Array specifying the (x,y,z) position of each spherical
%                       body.
%
% rho                -> Rather than provide a list of spherical bodies, an
%                       arbitrary density field rho may be provided instead.
%                       This should be passed as a function handle of (x,y,z) or
%                       as a vector of data values.
%
% boundary_condition -> Dirichlet boundary conditions on phi are applied, which
%                       may be specified as:
%                       'combined_onebody' => single body solution with the same
%                         total mass as the actual configuration
%                       'zero'          => 0
%                       'firstbody'     => single body solution of 1st body
%                       'sum_onebodies' => sum of single body solutions
%                       'custom'        => condition must be set by bc_func
% bc_func            -> Function handle of (x,y,z) specifying Dirichlet boundary
%                       condition on phi.
%
%
% Iteration Parameters
%
% max_iter           -> Maximum iterations to perform before exiting the
%                       computation. Note that if the residual ceases to
%                       decrease, the computation may exit before max_iter is
%                       reached.
%
% solution_method_switch_limit -> maximum times the solver should switch
%                       between least-squares and exact residual methods before
%                       ending.
%
% least_square_iterations -> when the method switches to least-squares, how many
%                       iterations should be performed before switching back.
%                       Default: 2
%
% nu                 -> Iterative update size (should be between 0 and 1), i.e.,
%                       step size in the gradient descent of the residual.
%                       nu only needs to be set manually if res_increase is set
%                       to 'abort' or 'ignore'; in general, it is better to let
%                       the code choose its own value of nu by setting
%                       res_increase='adaptive' and optionally specifying
%                       min_nu_limit and max_nu_limit.
% min_nu_limit       -> Lower limit on nu. If not set, the adaptive method is
%                       allowed to choose arbitrarily small nu values.
% max_nu_limit       -> Upper limit on nu. If not set, the adaptive method is
%                       allowed to choose arbitrarily large nu values.
%
% res_increase       -> 'adaptive': decrease nu as necessary (default)
%                       'abort': abort when the integrated residual increases
%                       'ignore': ignore residual increases (not recommended)
%
% linear_solver      -> 'bicg': solve linear problem iteratively with
%                         biconjugate gradient method. (default)
%                       'direct': solve linear problem exactly with mldivide
%                       For large problems (hence most 3D problems), bicg is
%                       much faster and is recommended.
%
% initial_condition  -> 'combined_onebody' => use a single body solution with
%                         the same total mass as the actual configuration
%                       'firstbody'     => use single body solution of 1st body
%                       'sum_onebodies' => use sum of single body solutions
%                       'newton'        => use Newtonian solution
%                       'continue'      => use existing workspace phiv
%                       'custom'        => condition must be set by ic_func
% ic_func            -> Function handle of (x,y,z) specifying initial condition.
% k_init             -> If initial_condition='newton', then the user may set
%                       k_init, to be used as: k_init lap phi = rho.
%                       If k_init is not set, the existing value of k will be
%                       used. If k=0 and initial_condition='newton', then
%                       k_init is required.
%
%
% Mesh Parameters
%
% meshes             -> A cell list of cell-tuples describing x, y, and z
%                       coordinates of meshes on which to solve. For a simple
%                       grid, meshes={{x_list, y_list, z_list}} will suffice,
%                       where x_list, y_list, z_list are lists of grid
%                       coordinates. For multi-scale meshes, a series of triples
%                       may be used as
%                       meshes={ {x_list1, y_list1, z_list1},
%                                {x_list2, y_list2, z_list2},...}.
%                       For multiscale meshes, each finer mesh must be fully
%                       contained within a larger mesh (no partial overlaps) and
%                       must have spacing which divides the larger mesh spacing.
% alignment_tol      -> Alignment tolerance, the distance below which two nearby
%                       points should be considered a single point.
%                       Default: 1E-12
% import_mesh        -> Boolean, whether mesh data should be imported from a
%                       separate mat file to avoid recalculation.
% import_mesh_file   -> Name of mat file from which to import mesh_data, if
%                       import_mesh is true.
%
%
% Ramp Parameters
% A computation may proceed by ramping the value of k. This may improve
% convergence in some cases, although it should not be needed in general.
%
% use_k_ramp         -> If true, k should be ramped through the list kvals.
% kvals              -> List of values k should take
% ramp_iter          -> Number of iterations after which to go to next ramp value
% post_ramp_iter     -> Number of iterations to perform after the ramping ends
%
%
% Data Output Parameters
%
% save_output        -> Whether or not to save data and figures upon completion
% clobber            -> True => overwrite existing data files. False => don't.
% data_dir           -> Directory for saved data and figures.
% output_filename    -> Filename to use for the saved data (should be specific
%                       to the run's parameters).
%
%
% Plotting Parameters
% vain_findiff_3d has a few basic plots defined to watch while the solver
% iterates and to save as output. These can cover the entire range or a
% zoomed-in range ("close region") which is automatically set to surround one of
% the bodies.
%
% plot_while_solving -> Whether or not to display plots while the solver runs.
%                       Drawing plots can get very slow on large simulations.
% plot_update_steps  -> If plot_while_solving=true and plot_update_steps=n, then
%                       update plots every nth iteration.
% plots_to_show      -> List of plots to show while solving.
% plots_to_save      -> List of plots to save when the computation ends.
%                       Note that plots_to_save does not have to contain or
%                       intersect plots_to_show.
% use_custom_close   -> If true, use a manually set close region instead of the
%                       automatically chosen one.
% xc_close           -> x-coordinate of center of close region.
% yc_close           -> y-coordinate of center of close region.
% zc_close           -> z-coordinate of center of close region.
% xl_close           -> Extent of close region in x.
% yl_close           -> Extent of close region in y.
% zl_close           -> Extent of close region in z.
%

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2018-09-14
% Updated: 2019-09-10

% Copyright 2020 Nicholas C. White
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

addpath('tools');

% Clear functions and input variables
clear functions
clear alignment_tol bc_func body_coordinates body_densities body_masses body_radii boundary_condition clobber data_dir ic_func initial_condition k k_init kvals linear_solver max_iter max_nu_limit import_mesh import_mesh_file meshes min_nu_limit nu output_filename plot_update_steps plot_while_solving plots_to_save plots_to_show post_ramp_iter ramp_iter res_increase rho solution_method_switch_limit least_square_iterations use_custom_close use_k_ramp xc_close xl_close yc_close yl_close zc_close zl_close

% Set input parameters from param_file
if ((~exist('param_file','var')) || ~ischar(param_file));
    disp('Specify a parameter file, e.g.');
    disp('   param_file = ''param_files/param01.m'';');
    disp(' ');
    if exist('param_files', 'dir')
        if length(dir('param_files/*.m'))>0
            fprintf('List of param_files/ contents:');
            dir('param_files/*.m')
        end
    end
    error('param_file not set.');
elseif (~exist(param_file,'file'));
    error(['Parameter file ''' param_file ''' does not exist.']);
end;

run(param_file);

if ((~exist('save_output','var')) || ~islogical(save_output));
    disp('Set save_output to true or false before running.');
    error('save_output not set.');
end

%%%%%%%%%% ITERATION AND MESH PARAMETERS %%%%%%%%%%


fprintf('\n');
disp('===Iteration settings===');

ramp_length = 1;

if (~exist('use_k_ramp','var'));
    use_k_ramp = false;
end
cadp_names = { 'use_k_ramp' };
cadp_attributes = { {'logical'}, {'scalar'} };
check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

if (use_k_ramp)
    % Check ramp iteration parameters
    cadp_names = {...
                    'ramp_iter',...
                    'post_ramp_iter',...
                  };
    cadp_attributes = {...
                        {'numeric'}, {'integer', 'scalar', 'nonnegative'};...
                        {'numeric'}, {'integer', 'scalar', 'nonnegative'};...
                      };
    check_and_disp_params
    clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
else
    % If we don't ramp, then max_iter should be set instead
    cadp_names = {'max_iter'};
    cadp_attributes = {{'numeric'}, {'integer', 'scalar', 'nonnegative'}};
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
end

if (~exist('solution_method_switch_limit', 'var'))
    solution_method_switch_limit = 4;
end
if (~exist('least_square_iterations', 'var'))
    least_square_iterations = 2;
end
cadp_names = {'solution_method_switch_limit', 'least_square_iterations'};
    cadp_attributes = {...
                        {'numeric'}, {'integer', 'scalar', 'nonnegative'};...
                        {'numeric'}, {'integer', 'scalar', 'nonnegative'};...
                      };
check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

% Check other iteration parameters

if (~exist('res_increase','var'));
    res_increase = 'best';
end
cadp_names = {...
                'res_increase',...
              };
cadp_attributes = {...
                    'validatestring', {'ignore', 'best', 'adaptive', 'abort'};...
                  };
check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

% If res_increase is adaptive, we need to define the smallest allowable step
if (strcmpi(res_increase, 'adaptive') || strcmpi(res_increase, 'best'))
    cadp_names = {...
                    'min_nu_limit',...
                  };
    cadp_attributes = {...
                        {'numeric'}, {'scalar', 'positive'};...
                      };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

    % Optionally set max_nu_limit also
    if ~exist('max_nu_limit', 'var')
        max_nu_limit = 1;
    else
        cadp_names = {...
                        'max_nu_limit',...
                      };
        cadp_attributes = {...
                            {'numeric'}, {'scalar', 'positive'};...
                          };
        check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
    end
else
    min_nu_limit = NaN;
    cadp_names = {...
                    'nu',...
                  };
    cadp_attributes = {...
                        {'numeric'}, {'scalar', 'nonnegative'};...
                      };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
end

if (~exist('linear_solver','var'));
    linear_solver = 'bicg';
end
cadp_names = {...
                'linear_solver',...
              };
cadp_attributes = {...
                    'validatestring', {'direct', 'bicg'};...
                  };
check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

fprintf('\n');
disp('===Mesh settings===');

if (~exist('import_mesh','var'));
    import_mesh = false;
end
if import_mesh
    assert(exist('import_mesh_file','var')==1, 'Because import_mesh = true, you must specify import_mesh_file, a mat file from which to import the mesh data.');
    cadp_names = { 'import_mesh_file' };
    cadp_descriptions = { 'mat file containing mesh data' };
    cadp_attributes = { {'char'}, {'vector', 'nonempty'} };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

    % If the mesh import file doesn't exist, just give a warning and then generate the mesh as usual
   if exist(import_mesh_file,'file')

        mesh_data = load(import_mesh_file);
        try
            meshes = mesh_data.meshes;
            idgrids = mesh_data.idgrids;
            nodedata = mesh_data.nodedata;
            dxv = mesh_data.dxv;
            d2xv = mesh_data.d2xv;
            dyv = mesh_data.dyv;
            d2yv = mesh_data.d2yv;
            dzv = mesh_data.dzv;
            d2zv = mesh_data.d2zv;
            int_weight = mesh_data.int_weight;
            if any(ismember(fields(mesh_data), 'alignment_tol'))
                alignment_tol = mesh_data.alignment_tol;
            end
        catch ME
            warning(['Error loading mesh data from ' import_mesh_file '. This file must contain the following variables: meshes, idgrids, nodedata, dxv, d2xv, dyv, d2yv, dzv, d2zv, int_weight.']);
            rethrow(ME);
        end
        clearvars mesh_data;
    else
        warning(sprintf('Mesh import file (%s) not found! Generating mesh instead.', import_mesh_file));
        import_mesh = false;
    end
end

if ~import_mesh
    cadp_names = { 'meshes' };
    cadp_descriptions = { 'input mesh' };
    cadp_attributes = { {'cell'}, {'vector'} };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
end
% TODO: more validation of mesh input

if (~exist('alignment_tol','var'));
    alignment_tol = 1E-12;
end
cadp_names = { 'alignment_tol' };
cadp_descriptions = { 'alignment tolerance' };
cadp_attributes = { {'numeric'}, {'scalar', 'positive'} };
check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% PHYSICAL PARAMETERS %%%%%%%%%%
fprintf('\n');
disp('===Physical paramaeters===');

if (use_k_ramp)
    % kvals, a list of k values must be set
    cadp_names = { 'kvals' };
    cadp_attributes = { {'numeric'}, {'vector', 'nonempty'} };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

    % set k access function
    kval = @(j) kvals(max(1,min(j, length(kvals))));
    ramp_length = length(kvals);

    k = kvals(end);
else
    % a single k must be set
    cadp_names = { 'k' };
    cadp_attributes = { {'numeric'}, {'scalar'} };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

    % set k access function
    kval = @(j) k;
end

% The user can either provide data for spherical bodies, or provide a custom forcing rho
hasBodies = (exist('body_radii', 'var') || exist('body_densities', 'var') || exist('body_masses', 'var') || exist('body_coordinates', 'var'));
% Make sure that there's no ambiguity
if (hasBodies)
    assert(exist('body_radii', 'var') && (exist('body_densities', 'var') || exist('body_masses', 'var')) && exist('body_coordinates', 'var'), 'To define spherical bodies, the variables body_radii, body_densities/body_masses, and body_coordinates must all be set.');
    assert(~exist('rho', 'var'), 'Ambiguous forcing. Either spherical bodies may be defined with body_radii, body_densities/body_masses, and body_coordinates, or the forcing rho may be specified; not both.');
    assert(~(exist('body_densities', 'var') && exist('body_masses', 'var')), 'Ambiguous forcing. EIther body_densities or body_masses may be defined; not both.');
else
    assert(exist('rho', 'var')==1, 'No forcing specified. Either spherical bodies may be defined with body_radii, body_densities/body_masses, and body_coordinates, or the forcing rho may be specified.');
end

if hasBodies
    if exist('body_densities', 'var')
        % body_densities is defined
        cadp_names = {...
                        'body_densities',...
                      };
        cadp_descriptions = {...
                        'Spherical body densities',...
                      };
        cadp_attributes = {...
                            {'numeric'}, {'vector'};...
                          };
        check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
        hasDensities = true;
    else
        % body_masses must be defined
        cadp_names = {...
                        'body_masses',...
                      };
        cadp_descriptions = {...
                        'Spherical body masses',...
                      };
        cadp_attributes = {...
                            {'numeric'}, {'vector'};...
                          };
        check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
        hasDensities = false;
    end

    cadp_names = {...
                    'body_radii',...
                    'body_coordinates',...
                  };
    cadp_descriptions = {...
                    'Spherical body physical radii',...
                    'Spherical body locations [x,y,z]',...
                  };
    cadp_attributes = {...
                        {'numeric'}, {'vector', 'positive'};...
                        {'numeric'}, {'size', [NaN 3]};...
                      };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

    nBodies = length(body_radii);

    if hasDensities
        assert((length(body_densities)==length(body_radii)) && (length(body_densities)==size(body_coordinates,1)), 'body_radii, body_densities, and body_coordinates must all have the same length (that length being the number of bodies).')
        body_masses = (4/3)*pi*body_densities.*body_radii.^3;
    else
        assert((length(body_masses)==length(body_radii)) && (length(body_masses)==size(body_coordinates,1)), 'body_radii, body_masses, and body_coordinates must all have the same length (that length being the number of bodies).')
        body_densities = body_masses./( (4/3)*pi*body_radii.^3 );
    end

    % Ensure certain variables are column vectors
    columnize = @(vec) reshape(vec, length(vec), 1);
    body_masses = columnize(body_masses);
    body_densities = columnize(body_densities);
    body_radii = columnize(body_radii);

    Rvs = body_radii.* (8*body_densities/(3*k^2)).^(1/3); % Vainshtein radius
    Rv_max = max(Rvs); % Max Vainshtein radius
    Rss = body_radii;
    rhos = body_densities;
    xbs = body_coordinates(:,1);
    ybs = body_coordinates(:,2);
    zbs = body_coordinates(:,3);
else
    Rss = [];
    Rvs = [];
    xbs = [];
    ybs = [];
    zbs = [];
    nBodies = 0;
end

custom_forcing = false;
if exist('rho', 'var')
    custom_forcing = true;
    % Custom forcing rho; we'll parse this later
end

% Read boundary condition
assert(exist('boundary_condition', 'var')==1, 'boundary_condition must be specified. Options are: ''combined_onebody'', ''zero'', ''firstbody'', ''sum_onebodies'', and ''custom''');
cadp_names = { 'boundary_condition' };
cadp_attributes = { 'validatestring', {'combined_onebody', 'zero', 'firstbody', 'sum_onebodies', 'custom'} };
check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

if strcmpi(boundary_condition, 'custom')
    assert(exist('bc_func', 'var')==1, 'To use a custom boundary condition, you must define bc_func(r,z).');
    cadp_names = { 'bc_func' };
    cadp_attributes = { {'function_handle'}, {'scalar'} };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
end

use_ramp = (use_k_ramp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% INITIAL CONDITION %%%%%%%%%%

fprintf('\n');
disp('===Initial condition settings===');
cadp_names = { 'initial_condition' };
cadp_attributes = { 'validatestring', {'Newton', 'combined_onebody', 'firstbody', 'sum_onebodies', 'continue', 'custom'} };
check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

if strcmpi(initial_condition, 'continue')
    assert(exist('phiv','var')==1, 'To continue from a previous simulation, the variable phiv must still exist.');
end

if strcmpi(initial_condition, 'combined_onebody')
    assert(hasBodies, 'To use ''combined_onebody'' as an initial condition, a spherical body must be defined (instead of a custom rho forcing)');
end
if strcmpi(initial_condition, 'firstbody')
    assert(hasBodies, 'To use ''firstbody'' as an initial condition, a spherical body must be defined (instead of a custom rho forcing)');
end
if strcmpi(initial_condition, 'sum_onebodies')
    assert(hasBodies, 'To use ''sum_onebodies'' as an initial condition, a spherical body must be defined (instead of a custom rho forcing)');
end
if strcmpi(initial_condition, 'custom')
    assert(exist('ic_func', 'var')==1, 'To use a custom initial condition, you must define ic_func(x,y,z).');
    cadp_names = { 'ic_func' };
    cadp_attributes = { {'function_handle'}, {'scalar'} };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
end

% The user may optionally set k_init for the initial condition (only Newton currently)
if strcmpi(initial_condition, 'Newton')
    if ~exist('k_init', 'var')
        k_init = k;
    end
    cadp_names = { 'k_init' };
    cadp_attributes = { {'numeric'}, {'scalar', 'positive'} };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% OUTPUT FILENAME AND PLOT SETTINGS %%%%%%%%%%

fprintf('\n');
disp('===Output settings===');
if (save_output)
    cadp_names = {...
                    'data_dir',...
                    'output_filename',...
                    'clobber',...
                  };
    cadp_descriptions = {...
                    'data directory',...
                    'output filename base',...
                    'whether existing files should be clobbered (overwritten)',...
                  };
    cadp_attributes = {...
                        {'char'}, {'vector', 'nonempty'};...
                        {'char'}, {'vector', 'nonempty'};...
                        {'logical'}, {'scalar'};...
                      };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
end

if ~exist('plot_while_solving', 'var')
    plot_while_solving = true;
end
if ~exist('plot_update_steps', 'var')
    plot_update_steps = 5;
end
if ~exist('plots_to_show', 'var')
    % plots to show while solving or upon completion (must be subset of valid_plot_types)
    plots_to_show = {'residual_history'};
end
if ~exist('plots_to_save', 'var')
    % plots to save upon completion (must be subset of valid_plot_types)
    plots_to_save = {'all'};
end
cadp_names = {...
                'plot_while_solving',...
                'plot_update_steps',...
                'plots_to_show',...
                'plots_to_save',...
              };
cadp_attributes = {...
                    {'logical'}, {'scalar'};...
                    {'numeric'}, {'scalar', 'positive'};...
                    {'cell'}, {};...
                    {'cell'}, {};...
                  };
check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');

% Custom close coordinates
if ~exist('use_custom_close', 'var')
    use_custom_close = false;
end
if use_custom_close
    assert(exist('xc_close', 'var') && exist('xl_close', 'var') && exist('yc_close', 'var') && exist('yl_close', 'var') && exist('zc_close', 'var') && exist('zl_close', 'var'), 'To use custom close coordinates, xc_close, xl_close, yc_close, yl_close, zc_close, and zl_close must be defined.');
    cadp_names = {...
                    'xc_close',...
                    'xl_close',...
                    'yc_close',...
                    'yl_close',...
                    'zc_close',...
                    'zl_close',...
                  };
    cadp_attributes = {...
                        {'numeric'}, {'scalar'};...
                        {'numeric'}, {'scalar'};...
                        {'numeric'}, {'scalar'};...
                        {'numeric'}, {'scalar'};...
                        {'numeric'}, {'scalar'};...
                        {'numeric'}, {'scalar'};...
                      };
    check_and_disp_params; clear('cadp_names', 'cadp_descriptions', 'cadp_attributes');
end

fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure max_iter is set
if use_ramp
    max_iter = ramp_length*ramp_iter + post_ramp_iter;
end

% Check that plots_to_show and plots_to_save are subsets of allowed plots
valid_plot_types = { 'residual_history' };
if (length(plots_to_save) == 1) && strcmpi(plots_to_save{1}, 'all')
    plots_to_save = valid_plot_types;
end

if (length(plots_to_show)~=length(intersect(valid_plot_types, plots_to_show)));
    warning('plots_to_show contains invalid entries, which will be ignored.');
    plots_to_show = intersect(valid_plot_types, plots_to_show);
end
if (length(plots_to_save)~=length(intersect(valid_plot_types, plots_to_save)));
    warning('plots_to_save contains invalid entries, which will be ignored.');
    plots_to_save = intersect(valid_plot_types, plots_to_save);
end

% construct list of extra plots (to save but not show)
plots_noshow = setdiff(plots_to_save, plots_to_show);

if (nBodies < 2) || (body_densities(1) == 0 || all(body_densities(2:end) == 0))
    disp('Note: since two bodies are not both defined, Hiramatsu''s screening statistic will not be plotted.');
    plots_to_save(strcmp(plots_to_save,'hm_screen_stat')) = [];
    plots_to_show(strcmp(plots_to_show,'hm_screen_stat')) = [];
    plots_noshow(strcmp(plots_noshow,'hm_screen_stat')) = [];
end

if (~hasBodies) || (k <= 0)
    % TODO implement Newton gravity for arbitrary rho
    disp('Note: since spherical bodies not defined, the Newtonian gravity solution will not be computed.');
    plots_to_save(strcmp(plots_to_save,'newton_frac')) = [];
    plots_to_show(strcmp(plots_to_show,'newton_frac')) = [];
    plots_noshow(strcmp(plots_noshow,'newton_frac')) = [];
    plots_to_save(strcmp(plots_to_save,'newton_frac_close')) = [];
    plots_to_show(strcmp(plots_to_show,'newton_frac_close')) = [];
    plots_noshow(strcmp(plots_noshow,'newton_frac_close')) = [];
end

%%%%%%%%%% Set up finite difference method %%%%%%%%%%

disp('Setting up mesh...');

% Start r half a unit below 0 instead of at 0 (despite it being the radial
% coordinate), since this makes the Neumann bc easier to implement.
% Another good idea from Hiramatsu.
% TODO: Ensure that x begins below 0

mesh_time = NaN;
if ~import_mesh
    tic
    [idgrids, nodedata, dxv, d2xv, dyv, d2yv, dzv, d2zv, int_weight] = multimesh_3d(meshes, 'alignmentTol', alignment_tol, 'edgeOrder', 2, 'progressBar', true);
    mesh_time = toc;
end
nodeN = length(nodedata.x);
nodes = [1:nodeN]';
xv = nodedata.x';
yv = nodedata.y';
zv = nodedata.z';
xmin = min(xv);
xmax = max(xv);
ymin = min(yv);
ymax = max(yv);
zmin = min(zv);
zmax = max(zv);

approx_eq = @(x1, x2) abs(x1 - x2) < alignment_tol;

% Define a shorthand function for constructing a sparse diagonal matrix
% from a vector.
sdiag = @(v) spdiags(v, 0, nodeN, nodeN);

% Define integration weight matrix, W.
% We will need this to compute adjoints of operators: adj(A) = inv(W)*transpose(A)*transpose(W)
W = sdiag(int_weight);

% Define integration operator
trapez_int = @(fieldv) sum(W*fieldv);


% Laplacian operator, since it is frequently used.
lapv = d2xv + d2yv + d2zv;

disp('Setting up BC matrices...');

% Just a field with value 1 where boundary conditions should be applied
% and value 0 everywhere else.
bcindv = nodes( approx_eq(xv, xmin) | approx_eq(xv, xmax) | approx_eq(yv, ymin) | approx_eq(yv, ymax) | approx_eq(zv, zmin) | approx_eq(zv, zmax)  );
nindv = setdiff(nodes, bcindv);

dir_bcindv = bcindv; % Dirichlet bc everywhere
dir_nindv = setdiff(nodes, dir_bcindv); % all nodes except those at the dirichlet boundaries

sdiagr = @(v) spdiags(v, 0, length(nindv), length(nindv));

% Homogeneous Dirichlet condition at box edges
dirv = zeros(nodeN, 1);
dirv( dir_bcindv ) = 1;
dirvop = sdiag(dirv);

% Combined boundary condition matrix.
%bcv = dirvop + neumvop;
bcv = dirvop;

% Reduced weight matrix
Wr = W(nindv, nindv);
%Winvr = inv(Wr);

% Reduced derivative matrices
dxvr = add_operator_bc(bcv, bcindv, dxv);
dyvr = add_operator_bc(bcv, bcindv, dyv);
dzvr = add_operator_bc(bcv, bcindv, dzv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Solve for CGG field %%%%%%%%%%

% Set up density field (forcing)
if custom_forcing
    if isa(rho, 'function_handle')
        rhov = rho(xv, yv, zv);
        assert( all(size(rhov) == size(xv)) );
    elseif all(size(rho) == size(xv))
        rhov = rho;
    elseif all(size(rho) == size(xv'))
        rhov = rho';
    else
        error('Input rho must be a function of x, y, and z, or a vector of length nodeN (%d).', nodeN);
    end
else
    assert(hasBodies);

    bodyDistvs = {};
    for bi = 1:nBodies
        bodyDistvs{bi} = sqrt((xv - xbs(bi)).^2 + (yv - ybs(bi)).^2 + (zv - zbs(bi)).^2);
    end

    rhovs = {};
    rho_origs = [];
    masses = [];
    rhov = 0;
    for bi = 1:nBodies
        this_rhov = rhos(bi)*(bodyDistvs{bi} <= Rss(bi));
        this_mass = (4/3)*pi*Rss(bi)^3*rhos(bi);
        rho_origs(bi) = rhos(bi);

        if (this_mass ~= 0)
            this_mass_num = trapez_int(this_rhov);
            if (this_mass_num == 0)
                warning(sprintf('Body %d does not overlap any mesh points. Continuing anyway.', bi));
            else
                if (this_mass_num/this_mass) <= 2/3
                    warning(sprintf('Body %d does not overlap many mesh points. Numerical volume is only %0.1f%% of true volume.', bi, 100*this_mass_num/this_mass));
                end
                rhos(bi) = rho_origs(bi)*(this_mass/this_mass_num);
                this_rhov = rhos(bi)*(bodyDistvs{bi} <= Rss(bi));
            end
        end

        masses(bi) = this_mass;
        rhovs{bi} = this_rhov;

        rhov = rhov + rhovs{bi};
    end
end

%assert(max(abs(rhov)) > 0, 'Density field is 0 everywhere.');

% A reference node at z=0 and r=rmax
%TODO: REDEFINE farthest_i, OR, BETTER YET, JUST DITCH IT ALTOGETHER
bodyDist_sumv = 0;
for bi = 1:nBodies
    bodyDist_sumv = bodyDist_sumv + bodyDistvs{bi};
end
[~,farthest_i] = max(bodyDist_sumv);

% Residual function
Resf = @(phiv,k) lapv*phiv + (k/2) - sqrt( (k/2)^2 + rhov + (d2xv*phiv).^2 + (d2yv*phiv).^2 + (d2zv*phiv).^2 + 2*(dxv*dyv*phiv).^2 + 2*(dyv*dzv*phiv).^2 + 2*(dzv*dxv*phiv).^2  );

% Objective function
Objf = @(phiv,k) Resf(phiv,k)'*W*Resf(phiv,k);

% Linear operator (function of phiv)
Lopf = @(phiv, k) lapv - sdiag(1./sqrt( (k/2)^2 + rhov + (d2xv*phiv).^2 + (d2yv*phiv).^2 + (d2zv*phiv).^2 + 2*(dxv*dyv*phiv).^2 + 2*(dyv*dzv*phiv).^2 + 2*(dzv*dxv*phiv).^2  )) * ...
        ( sdiag(d2xv*phiv)*d2xv...
        + sdiag(d2yv*phiv)*d2yv...
        + sdiag(d2zv*phiv)*d2zv...
        + 2*sdiag(dxv*dyv*phiv)*dxv*dyv...
        + 2*sdiag(dyv*dzv*phiv)*dyv*dzv...
        + 2*sdiag(dzv*dxv*phiv)*dzv*dxv...
        );

% Energy density function
% TODO: do we need boundary terms?????
%endensf = @(phiv) 0.5*(k*phiv - (drv*phiv).*(drv*phiv) + (dzv*phiv).*(dzv*phiv) ) .* (lapv*phiv) - rhov.*phiv;
%endensf = @(phiv) -0.5*k*((drv*phiv).*(drv*phiv) + (dzv*phiv).*(dzv*phiv) ) ...
            %+ (drv*phiv).*(drv*phiv).*(d2rv*phiv) + 2*(drv*phiv).*(dzv*phiv).*(dzv*drv*phiv) + (dzv*phiv).*(dzv*phiv).*(d2zv*phiv) ...
            %- rhov.*phiv;
endensf = @(phiv) ( (k/2)*phiv + lapv*phiv) .* ( (dxv*phiv).*(dxv*phiv) + (dyv*phiv).*(dyv*phiv) + (dzv*phiv).*(dzv*phiv) ) + rhov.*phiv;

disp('Setting up initial condition...');

% Set up one body solutions for each body
newton_psiv = zeros(nodeN, 1);
if hasBodies
    phi_onevs = {};
    lap_phi_onevs = {};
    phi_one_revvs = {};
    phi_sumv = 0;
    lap_phi_sumv = 0;
    for bi=1:nBodies
        phi_onevs{bi} = vainshteinone(k, rho_origs(bi), Rss(bi), bodyDistvs{bi}, false);
        lap_phi_onevs{bi} = vainshteinonelap(k, rho_origs(bi), Rss(bi), bodyDistvs{bi});
        phi_sumv = phi_sumv + phi_onevs{bi};
        lap_phi_sumv = lap_phi_sumv + lap_phi_onevs{bi};

        % Branch 2 solution
        phi_one_revvs{bi} = vainshteinone(k, rho_origs(bi), Rss(bi), bodyDistvs{bi}, false, true);
        phi_one_revvs{bi} = phi_one_revvs{bi} - phi_one_revvs{bi}(farthest_i) + phi_onevs{bi}(farthest_i);

        % Newtonian gravity
        if (k > 0)
            newton_psiv = newton_psiv + rho_origs(bi)/(4*k)*( (bodyDistvs{bi}.^2 - 3*Rss(bi)^2).*(bodyDistvs{bi} <= Rss(bi)) + (-2*Rss(bi)^3./bodyDistvs{bi}).*(bodyDistvs{bi} > Rss(bi)) );
        end
    end
end

% Set up single body solution with combined mass
if hasBodies
        mass_tot = sum(body_masses);

        % center of mass
        xyz_com = (body_masses'*body_coordinates)/mass_tot;
        x_com = xyz_com(1);
        y_com = xyz_com(2);
        z_com = xyz_com(3);

        R_combined = max(body_radii);
else
        % TODO: Warn about imperfect mass sum
        mass_tot = trapez_int(rhov);
        x_com = trapez_int(xv.*rhov)/mass_tot;
        y_com = trapez_int(yv.*rhov)/mass_tot;
        z_com = trapez_int(zv.*rhov)/mass_tot;

        x_com(isnan(x_com)) = 0;
        y_com(isnan(y_com)) = 0;
        z_com(isnan(z_com)) = 0;


        % TODO choose R_combined better
        R_rhomax = max(sqrt((xv(rhov>0)-x_com).^2 + (yv(rhov>0)-y_com).^2 + (zv(rhov>0)-z_com).^2));

        R_bounddist4 = min([min(abs(max(xv) - x_com), abs(min(xv) - x_com))/4, min(abs(max(yv) - y_com), abs(min(yv) - y_com))/4, min(abs(max(zv) - z_com), abs(min(zv) - z_com))/4]);

        R_combined = min(R_rhomax, R_bounddist4);
        if isempty(R_combined)
            R_combined = R_bounddist4;
        end
end
phi_combinedv = vainshteinone(k, mass_tot/((4/3)*pi*R_combined^3), R_combined, sqrt((xv-x_com).^2 + (yv-y_com).^2 + (zv-z_com).^2), false);

% Display a warning about the box being inside or outside the max Vainshtein radius
R_min_mass_edge_dist = min([min(abs(max(xv) - xv(rhov>0))), min(abs(xv(rhov>0) - min(xv))), min(abs(max(yv) - yv(rhov>0))), min(abs(yv(rhov>0) - min(yv))), min(abs(max(zv) - zv(rhov>0))), min(abs(zv(rhov>0) - min(zv)))]);
R_max_mass_edge_dist = max([max(abs(max(xv) - xv(rhov>0))), max(abs(xv(rhov>0) - min(xv))), max(abs(max(yv) - yv(rhov>0))), max(abs(yv(rhov>0) - min(yv))), max(abs(max(zv) - zv(rhov>0))), max(abs(zv(rhov>0) - min(zv)))]);
if exist('Rv_max', 'var')
    if (R_max_mass_edge_dist < 0.1*Rv_max)
        if strcmpi(boundary_condition, 'sum_onebodies')
            warning('Simulation bounds are below 10% of the max Vainshtein radius. In this nonlinear-dominated regime, ''sum_onebodies'' is not the ideal boundary condition. Consider ''combined_onebody'' instead.');
        else
            disp('Simulation bounds below 10% of the max Vainshtein radius.');
        end
    elseif (R_min_mass_edge_dist > 10*Rv_max)
        if strcmpi(boundary_condition, 'combined_onebody')
            warning('Simulation bounds are above 10x of the max Vainshtein radius. In this linear-dominated regime, ''combined_onebody'' is not the ideal boundary condition. Consider ''sum_onebodies'' instead.');
        else
            disp('Simulation bounds above 10 times the max Vainshtein radius.');
        end
    else
        if strcmpi(boundary_condition, 'combined_onebody') || strcmpi(boundary_condition, 'sum_onebodies')
            warning('Simulation bounds between 0.1 and 10 times the max Vainshtein radius. It is not clear whether superposition or combined boundary condition is better.');
        else
            warning('Simulation bounds between 0.1 and 10 times the max Vainshtein radius.');
        end
    end
end

% Define "reference solution" (phiA)
% TODO: allow custom reference solution ?
has_phi_ref = false;
phiv_refs = {};
phi_ref_titles = {};
if hasBodies && (nBodies == 1)
    phiv_refs = {phi_onevs{1}};
    phi_ref_titles = {'phiA'};
    has_phi_ref = true;
end
if ~has_phi_ref
    % remove phidiff from allowed plot strings
    disp('Note: since no reference plot is defined, phidiff plot will not be shown.');
    plots_to_save(strcmp(plots_to_save,'phidiff')) = [];
    plots_to_show(strcmp(plots_to_show,'phidiff')) = [];
    plots_noshow(strcmp(plots_noshow,'phidiff')) = [];
end

% set up initial conditions
if (strcmpi(initial_condition, 'Newton'))
    phiv = zeros(nodeN, 1);
    % Note: Newton condition will be computed after BCs
elseif (strcmpi(initial_condition, 'firstbody'))
    phiv = phi_onevs{1};
elseif (strcmpi(initial_condition, 'combined_onebody'))
    phiv = phi_combinedv;
elseif (strcmpi(initial_condition, 'sum_onebodies'))
    phiv = phi_sumv;
elseif (strcmpi(initial_condition, 'continue'))
    assert(exist('phiv','var')==1);
elseif (strcmpi(initial_condition, 'custom'))
    phiv = ic_func(xv, yv, zv);
else
    error('Invalid initial_condition string');
end

% We must enforce the boundary conditions on the initial conditions, since all corrections will be zero

% Enforce Dirichlet bcs at box edges
if strcmpi(boundary_condition, 'combined_onebody')
    phiv(dir_bcindv) = phi_combinedv(dir_bcindv);
elseif strcmpi(boundary_condition, 'sum_onebodies')
    phiv(dir_bcindv) = phi_sumv(dir_bcindv);
elseif strcmpi(boundary_condition, 'firstbody')
    phiv(dir_bcindv) = phi_onevs{1}(dir_bcindv);
elseif strcmpi(boundary_condition, 'zero')
    phiv(dir_bcindv) = 0;
elseif strcmpi(boundary_condition, 'custom')
    bc_datav = bc_func(xv, yv, zv);
    phiv(dir_bcindv) = bc_datav(dir_bcindv);
else
    error('boundary_condition incorrectly set');
end

% Set Newton initial condition
if strcmpi(initial_condition, 'Newton')

    newt_bcv = phiv;
    newt_bcv(dir_nindv) = 0;
    [klapvr, newt_Kvr, ~, add_rhs_factor, add_phi_factor] = add_operator_bc(bcv, bcindv, k_init*lapv, newt_bcv);
    rhovr = rhov(nindv);

    phiv = add_phi_factor + newt_Kvr*mldivide(klapvr, rhovr + add_rhs_factor);
end

% Save initial condition
phi_initv = phiv;

% Resv is the last residual; we need to initialize it before the loop
Resv = Resf(phiv, kval(1));

% intres_hist is a history of the integrated residual (so we can plot it)
intres_hist = [];

% energy_hist is a history of the integrated energy (so we can plot it)
energy_hist = [];

% Set up intres_hist and energy_hist with initial condition
endensv = endensf(phi_initv);
energy = trapez_int(endensv.^2);
intres = Objf(phi_initv, kval(1));
intres_hist(end+1) = intres;
energy_hist(end+1) = energy;

% Set up close domain and color scheme (for plotting)
if (~use_custom_close)
    xc_close = 0;
    xl_close = 2;
    yc_close = 0;
    yl_close = 2;
    zc_close = 0;
    zl_close = 2;
    if (hasBodies) && nBodies > 0
        xc_close = xbs(1);
        xl_close = Rss(1)*20;
        yc_close = ybs(1);
        yl_close = Rss(1)*20;
        zc_close = zbs(1);
        zl_close = Rss(1)*20;
        if nBodies > 1
            RsMax = max(Rss);
            xc_close = xbs(2);
            xl_close = Rss(2)*20;
            yc_close = ybs(2);
            yl_close = Rss(2)*20;
            zc_close = zbs(2);
            zl_close = Rss(2)*20;
        end
    end
end

xv_close = xv( (xc_close - 0.5*xl_close <= xv) & (xv <= xc_close + 0.5*xl_close) & (yc_close - 0.5*yl_close <= yv) & (yv <= yc_close + 0.5*yl_close) & (zc_close - 0.5*zl_close <= zv) & (zv <= zc_close + 0.5*zl_close) );
yv_close = yv( (xc_close - 0.5*xl_close <= xv) & (xv <= xc_close + 0.5*xl_close) & (yc_close - 0.5*yl_close <= yv) & (yv <= yc_close + 0.5*yl_close) & (zc_close - 0.5*zl_close <= zv) & (zv <= zc_close + 0.5*zl_close) );
zv_close = zv( (xc_close - 0.5*xl_close <= xv) & (xv <= xc_close + 0.5*xl_close) & (yc_close - 0.5*yl_close <= yv) & (yv <= yc_close + 0.5*yl_close) & (zc_close - 0.5*zl_close <= zv) & (zv <= zc_close + 0.5*zl_close) );
i_close = nodes( (xc_close - 0.5*xl_close <= xv) & (xv <= xc_close + 0.5*xl_close) & (yc_close - 0.5*yl_close <= yv) & (yv <= yc_close + 0.5*yl_close) & (zc_close - 0.5*zl_close <= zv) & (zv <= zc_close + 0.5*zl_close) );
xrange_close = [min(xv_close) max(xv_close)];
yrange_close = [min(yv_close) max(yv_close)];
zrange_close = [min(zv_close) max(zv_close)];
assert(length(xv_close)==length(zv_close))
assert(length(yv_close)==length(zv_close))

%if min(length(zv_close)) < 1
    %disp('Note: due to the coordinates, no close plots will be constructed.');
    %plots_to_save(strcmp(plots_to_save,'hm_screen_stat')) = [];
    %plots_to_show(strcmp(plots_to_show,'hm_screen_stat')) = [];
    %plots_noshow(strcmp(plots_noshow,'hm_screen_stat')) = [];
    %plots_to_save(strcmp(plots_to_save,'rad_close')) = [];
    %plots_to_show(strcmp(plots_to_show,'rad_close')) = [];
    %plots_noshow(strcmp(plots_noshow,'rad_close')) = [];
    %plots_to_save(strcmp(plots_to_save,'lap_phi_close')) = [];
    %plots_to_show(strcmp(plots_to_show,'lap_phi_close')) = [];
    %plots_noshow(strcmp(plots_noshow,'lap_phi_close')) = [];
%end

f_jet = figure; rjet = jet; rjet=rjet(end:-1:1,:); close(f_jet); % set colormap to match Hiramatsu

% Set up figures
fs = containers.Map;
axs = containers.Map;
legend_ps = containers.Map;
legend_strs = containers.Map;

for j = 1:length(plots_to_show)
    plot_str = plots_to_show{j};
    fs(plot_str) = figure;
    axs(plot_str) = axes('Parent', fs(plot_str));

    legend_ps(plot_str) = [];
    legend_strs(plot_str) = {};
end

for j = 1:length(plots_noshow)
    plot_str = plots_noshow{j};
    legend_ps(plot_str) = [];
    legend_strs(plot_str) = {};
end

% For MATLAB to make a nice contour plot, we block out the body edges
eblockv = 0*xv + 1;
for bi = 1:nBodies
    eblockv( abs(bodyDistvs{bi}.^2 - Rss(bi)^2) <= Rss(bi)^2*0.5) = NaN;
end

should_exit = false;

disp('Solving for CGG field...');

% number of _recorded_ iterations
% (since we sometimes have to step backwards when the residual increases)
true_iter = 0;

if (strcmpi(res_increase, 'adaptive') || strcmpi(res_increase, 'best'))
    cur_nu = max_nu_limit;
else
    cur_nu = nu;
end

solution_method_switches = 0;
solve_inside = true;
method_iter_j = 0;

rampj = 1;
ramp_iter_j = 0;
prev_step_ramp = false;
iter_times = [];
% Iteration loop
for (iter = 0:max_iter)
    disp(sprintf('\nIteration %d', iter));

    % Set up the linear operator (jacobian)
    % We update it on every iteration.

    % Full operator
    Lopv = Lopf(phiv, kval(rampj));

    % Reduced operator ("reduced" => taking into account BC)
    [Lopvr, Kvr, ~] = add_operator_bc(bcv, bcindv, Lopv);
    if (solve_inside)
        lhs_op_r = Lopvr;
        Resvr = Resv(nindv);
        rhs_op_r = -Resvr;
    else
        lhs_op_r = Lopv(:,nindv);
        rhs_op_r = -Resv;
    end

    tic
    % Solve reduced linear equation
    if strcmpi(linear_solver, 'direct') || ~solve_inside
        gammavr = mldivide(lhs_op_r,rhs_op_r);
    elseif strcmpi(linear_solver, 'bicg')
        [gammavr, bicg_flag, relres] = bicgstab(lhs_op_r,rhs_op_r, 1E-9, 1000, sdiagr(diag(lhs_op_r)));
        if bicg_flag ~= 0 && relres > 1E-6
            disp('Falling back to direct method for A*x = b');
            gammavr = mldivide(lhs_op_r,rhs_op_r);
        end
    else
        assert(false, 'Invalid linear solver method.');
    end
    iter_time = toc;
    iter_times(end+1) = iter_time;
    disp(sprintf('Elapsed time is %f seconds.',iter_time));

    % Rebuild full gamma (with BC)
    gammav = Kvr*gammavr;

    % Increase or decrease nu as needed
    if strcmpi(res_increase, 'adaptive') || strcmpi(res_increase, 'best')

        cur_nu = 1;
        new_intres = Objf(phiv + cur_nu*gammav, kval(rampj));
        bigstep_intres = Objf(phiv + 2*cur_nu*gammav, kval(rampj));
        smallstep_intres = Objf(phiv + 0.5*cur_nu*gammav, kval(rampj));

        while (cur_nu >= min_nu_limit) && (cur_nu <= max_nu_limit)
            if (smallstep_intres < new_intres) && (smallstep_intres < bigstep_intres) && (cur_nu*0.5 >= min_nu_limit)
                cur_nu = cur_nu*0.5;
                new_intres = smallstep_intres;
                smallstep_intres = Objf(phiv + 0.5*cur_nu*gammav, kval(rampj));
                bigstep_intres = Inf;
            elseif (bigstep_intres < new_intres) && (cur_nu*2 <= max_nu_limit)
                cur_nu = cur_nu*2;
                new_intres = bigstep_intres;
                bigstep_intres = Objf(phiv + 2*cur_nu*gammav, kval(rampj));
                smallstep_intres = Inf;
            else
                break
            end
        end

        disp(sprintf('nu = %g', cur_nu));
        disp(sprintf('res = %g', new_intres));

    end

    % Update phiv, and save phiv_prev (in case the residual goes up and we need to revert)
    phiv_prev = phiv;
    phiv = phiv + cur_nu * gammav;

    % Compute new residual
    Resv = Resf(phiv, kval(rampj));

    % Compute objective function ||R||^2
    intres = Objf(phiv, kval(rampj));

    % Compute integrated energy
    endensv = endensf(phiv);
    energy = trapez_int(endensv.^2);

    should_update_history = true;

    % What to do if the residual increases
    if ((length(intres_hist) > 1) && (intres > intres_hist(end)) && (~prev_step_ramp))

        if (method_iter_j == 0) && (solution_method_switches >= 1)
            % Both methods have failed, so don't keep switching back and forth.
            initiate_method_switch = false;
        elseif (solution_method_switches < solution_method_switch_limit)
            initiate_method_switch = true;
        else
            initiate_method_switch = false;
        end

        if initiate_method_switch
            % switch solution method before giving up
            solve_inside = ~solve_inside;
            solution_method_switches = solution_method_switches + 1;
            phiv = phiv_prev;
            should_update_history = false;
            method_iter_j = 0;
            disp('Residual did not decrease. Switching solution method');
        else
            if strcmpi(res_increase, 'abort') || strcmpi(res_increase, 'best') % stop
                phiv = phiv_prev;
                Resv = Resf(phiv, kval(rampj));
                endensv = endensf(phiv);
                should_exit = true;
                should_update_history = false;
                disp('Residual did not decrease (to machine precision).');

                % Instead of aborting, move on to the next k value if possible
                if use_ramp
                    if rampj < ramp_length
                        ramp_iter_j = ramp_iter - 1;
                        should_exit = false;
                        should_update_history = true;
                    end
                end
            else % ignore
            end
        end
    end

    % ramp_iter_j is the number of iterations done so far at this ramp value
    ramp_iter_j = ramp_iter_j + 1;

    % Update residual and energy history
    if (should_update_history)
        intres_hist(end+1) = intres;
        energy_hist(end+1) = energy;

        true_iter = length(intres_hist)-1;

        % method_iter_j is the number of iterations done so far with this solution method
        method_iter_j = method_iter_j + 1;

        if (~solve_inside) && (method_iter_j > least_square_iterations)
            solve_inside = true;
            solution_method_switches = solution_method_switches + 1;
            method_iter_j = 0;
            disp('Switching solution method');
        end
    end

    if use_ramp
        % Update the ramp value every ramp_iter steps
        prev_step_ramp = false;
        if (iter>0 && mod(ramp_iter_j,ramp_iter) == 0 && rampj < ramp_length)
            rampj = rampj + 1;
            prev_step_ramp = true;
            if use_k_ramp
                disp(sprintf('New k value: %g', kval(rampj)));
            end
        end
    end

    if (iter == max_iter)
        should_exit = true;
    end

    %%%%%%%%%% PLOT RESULTS %%%%%%%%%%

    if ( (plot_while_solving && mod(iter,plot_update_steps)==0) || should_exit)

        % if should_exit and save_output, then show all the extra plots
        show_extra_plots = (should_exit && save_output);

        if (show_extra_plots)
            for j = 1:length(plots_noshow)
                plot_str = plots_noshow{j};
                fs(plot_str) = figure;
                axs(plot_str) = axes('Parent', fs(plot_str));
            end
        end

        if (axs.isKey('residual_history') && ishandle(fs('residual_history')))
            % Plot residual history
            this_ax = axs('residual_history');
            semilogy(this_ax, [0:true_iter], intres_hist);
            xlabel(this_ax, 'Iteration');
            ylabel(this_ax, '||Res||^2');
            title(this_ax, sprintf('Integrated objective: %g (iteration %d)', intres, true_iter));
            xlim(this_ax, [0 true_iter+1]);
            drawnow
        end

        if (axs.isKey('energy_hist') && ishandle(fs('energy_hist')))
            % Plot energy history
            this_ax = axs('energy_hist');
            semilogy(this_ax, [0:true_iter], energy_hist);
            xlabel(this_ax, 'Iteration');
            ylabel(this_ax, 'E');
            title(this_ax, sprintf('Total energy: %g (iteration %d)', energy, true_iter));
            xlim(this_ax, [0 true_iter+1]);
            drawnow
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (should_exit)
        break
    end
end
mean_iter_time = mean(iter_times);
sprintf('Mean time per iteration: %f seconds.', mean_iter_time);

% Save results

if (save_output)
    % Create directories if necessary
    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end

    if ispc; % use backslash for directories on Windows :(
        dir_sep = '\';
    else;
        dir_sep = '/';
    end;
    if strcmp(data_dir(end), '/') || strcmp(data_dir(end), '\')
        data_dir = data_dir(1:end-1);
    end
    full_dir = [data_dir dir_sep output_filename dir_sep];
    if ~exist(full_dir, 'dir')
        mkdir(full_dir);
    end

    if ((~clobber) && exist([full_dir output_filename '.mat'],'file'))
        warning(['File ' full_dir output_filename '.mat' ' exists. Not overwriting.']);
    else
        % Save data
        disp(sprintf('Saving data to %s...',[full_dir output_filename '.mat']));

        % Rename some variables before saving
        if exist('Rss', 'var')
            body_radii = Rss;
        end
        if exist('xbs', 'var') && exist('ybs', 'var') && exist('zbs', 'var')
            body_coordinates = [xbs, ybs, zbs];
        end
        if exist('x_com', 'var') && exist('y_com', 'var') && exist('z_com', 'var')
            center_of_mass_coordinates = [x_com, y_com, z_com];
        end

        % Save variables (have to run in a loop since some variables don't always exist)
        vars_to_save = {
            'xv','yv','zv',...
            'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax',...
            'max_iter', 'initial_condition',...
            'linear_solver',...
            'k','phi_initv',...
            'rhov', 'phiv', 'Resv', 'intres_hist',...
            'energy_hist', 'endensv',...
            'idgrids', 'nodedata', 'int_weight','nodes','nodeN','meshes'...
            'dxv', 'd2xv', 'dyv', 'd2yv', 'dzv', 'd2zv',...
            'lapv', 'bcv', 'bcindv',...
            'param_file',...
            'body_masses',...
            'mass_tot', 'R_combined', 'center_of_mass_coordinates',...
            'body_radii', 'rhos', 'rho_origs', 'body_coordinates',...
            'phi_onevs', 'phi_sumv', 'phi_combinedv', 'newton_psiv',...
            'mean_iter_time', 'iter_times', 'mesh_time'...
            };
        save_append = false;
        for j = 1:length(vars_to_save)
            this_var = vars_to_save{j};
            if exist(this_var,'var')
                if save_append
                    save([full_dir output_filename '.mat'], this_var, '-append');
                else
                    save([full_dir output_filename '.mat'], this_var);
                end
                save_append = true;
            end
        end

        % Save figures
        disp('Saving figures...');
        plot_strs = keys(fs);
        pj = 0;
        for j = 1:length(fs)
            plot_str = plot_strs{j};
            if ~isempty(intersect({plot_str}, plots_to_save))
                pj = pj + 1;
                disp(sprintf('  [%2d/%d]: %s', pj, length(fs), plot_str));
                fname = plot_str;
                try
                    savefig(fs(plot_str), [full_dir output_filename '_' fname]); % MATLAB fig
                    print(fs(plot_str), [full_dir output_filename '_' fname], '-depsc'); % eps
                    print(fs(plot_str), [full_dir output_filename '_' fname], '-dpng'); % png
                catch ME
                    warning(['Error saving figure ' plot_str]);
                end
            end
        end

    end
else
    warning('Output not saved!')
end

% Close extra plots
for j = 1:length(plots_noshow)
    plot_str = plots_noshow{j};
    if fs.isKey(plot_str)
        close(fs(plot_str))
    end
end
