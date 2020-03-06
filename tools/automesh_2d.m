% AUTOMESH_2D Construct a set of nested meshes around a collection of
%  spherical bodies.
%
% {MESHES} = AUTOMESH_2D(BODY_Z_COORDINATES, BODY_RADII, BOX_BOUNDS,...
%                        'maxMeshSize', MAX_MESH_SIZE,...
%                        'minMeshSize', MIN_MESH_SIZE,...
%                        'stopRadii', STOP_RADII,...
%                        'sequenceBoxSizeFraction', SEQ_BOX_FRACTION,...
%                        'sequenceMeshSizeFraction', SEQ_MESH_FRACTION)
%  TODO: MAX_LEVELS or END_FINENESS or something
%
%
% Input:
%
%  BODY_Z_COORDINATES must be a list of z coordinates of the bodies (one row
%  per body).
%  BODY_RADII must be a list of the spherical body radii.
%  BOX_BOUNDS should be a list of the form:
%   [min_r, max_r; min_z, max_z]
%  MAX_MESH_SIZE should be a number indicating the length scale of the coarsest
%  mesh. Default: box length / 100
%  STOP_RADII is the minimum number of radii away from each body center at which
%  to stop producing smaller meshes (default 1). For example, if STOP_RADII is
%  3, then the smallest mesh around each body will have diameter approximately
%  6*r, r being that body's radius.
%  SEQ_BOX_FRACTION is a number between 0 and 1 (default 0.5) indicating the
%  length scale of the next nested box. For example, with 0.5, the first box
%  will have size L, the next approximately size 0.5*L, then 0.25*L, etc.
%  SEQ_MESH_FRACTION is a number between 0 and 1 (default 0.5) indicating the
%  length scale of the next nested mesh's spacing. For example, with 0.5, the
%  first mesh will have spacing dr, the next 0.5*dr, then 0.25*dr, etc.
%  SEQ_MESH_FRACTION must be of the form 1/N for some integer N; otherwise it
%  will be rounded to the nearest fraction of that form.
%  MIN_MESH_SIZE is the minimum allowed mesh spacing. If a nested mesh would be
%  finer, then nesting stops there. Default: 0
%
%
% Output:
%
%  MESHES is a cell list of grids, which can then be passed to multimesh_2d.
%
%
% Example:
%
%  TODO
%

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-08-03
% Updated: 2019-08-03

function meshes = automesh_2d(varargin)

    parser = inputParser;
    parser.addRequired('bodyZCoordinates', @(x) isnumeric(x) && isvector(x));
    parser.addRequired('bodyRadii', @(x) isnumeric(x) && isvector(x) && all(x>0));
    parser.addRequired('boxBounds', @(x) isnumeric(x) && size(x, 1)==2 && size(x,2)==2);
    parser.addParameter('maxMeshSize', NaN, @(x) isscalar(x) && isnumeric(x) && (x>0));
    parser.addParameter('minMeshSize', 0, @(x) isscalar(x) && isnumeric(x) && (x>=0));
    parser.addParameter('stopRadii', 1, @(x) isscalar(x) && isnumeric(x) && (x>0));
    parser.addParameter('sequenceBoxSizeFraction', 0.5, @(x) isscalar(x) && isnumeric(x) && (x>0) && (x < 1));
    parser.addParameter('sequenceMeshSizeFraction', 0.5, @(x) isscalar(x) && isnumeric(x) && (x>0) && (x < 1));
    parse(parser, varargin{:});

    body_z_coords = parser.Results.bodyZCoordinates;
    body_zs = body_z_coords;
    if size(body_zs,2) > 1
        body_zs = body_zs';
    end
    body_rs = body_zs * 0;
    body_radii = parser.Results.bodyRadii;
    box_bounds = parser.Results.boxBounds;
    max_dmesh = parser.Results.maxMeshSize;
    min_dmesh = parser.Results.minMeshSize;
    stop_radii = parser.Results.stopRadii;
    seq_box_frac = parser.Results.sequenceBoxSizeFraction;
    seq_mesh_frac = parser.Results.sequenceMeshSizeFraction;

    seq_mesh_int = max(round(1/seq_mesh_frac), 2);
    if abs((seq_mesh_frac - 1/seq_mesh_int)/seq_mesh_frac) > 0.01
        warning(sprintf('Changing SEQ_MESH_FRACTION from %f to %f (1/%d)',seq_mesh_frac, 1/seq_mesh_int, seq_mesh_int));
    end
    seq_mesh_frac = 1/seq_mesh_int;

    Nbody = length(body_radii);
    assert(length(body_z_coords) == Nbody, ...
        'body_z_coords and body_radii must have the same length (one entry for each body)');
    r_bounds = box_bounds(1,:);
    assert(r_bounds(1) == 0, 'The lower r_bound must be 0');
    z_bounds = box_bounds(2,:);
    assert(all(box_bounds(:,2) - box_bounds(:,1) > 0), ...
        'box_bounds must describe a postive-volume region');
    Lr = r_bounds(2) - r_bounds(1);
    Lz = z_bounds(2) - z_bounds(1);
    L = min([Lr Lz]);

    if isnan(max_dmesh)
        max_dmesh = L/100;
    end


    % Construct the first mesh

    % We want coarse_dmesh to divide Lr and Lz. For now, we use a dumb brute force approach
    coarse_dmesh = max_dmesh;
    coarse_dmesh_guesses = coarse_dmesh*[1:-0.001:0.9];
    best_res = Inf;
    for j = 1:length(coarse_dmesh_guesses)
        this_dmesh = coarse_dmesh_guesses(j);
        Lr_res = ceil(Lr/this_dmesh) - (Lr/this_dmesh);
        Lz_res = ceil(Lz/this_dmesh) - (Lz/this_dmesh);
        this_res = abs(Lr_res) + abs(Lz_res);
        if this_res < best_res*0.98
            best_res = this_res;
            coarse_dmesh = this_dmesh;
        end
    end
    coarse_dr = Lr/ceil(Lr/coarse_dmesh);
    coarse_dz = Lz/ceil(Lz/coarse_dmesh);

    coarse_mesh_r = linspace(r_bounds(1), r_bounds(2), 1 + Lr/coarse_dr);
    coarse_mesh_z = linspace(z_bounds(1), z_bounds(2), 1 + Lz/coarse_dz);

    % Set up initial rss, zss, and a mesh chain for each body
    rss = {coarse_mesh_r};
    zss = {coarse_mesh_z};
    drs = [coarse_dr];
    dzs = [coarse_dz];

    body_mesh_chain = {};
    for jbody = 1:Nbody
        body_mesh_chain{jbody} = [1];
    end

    unskipped_bodies = 1:Nbody;
    layers = 1;
    while ~isempty(unskipped_bodies)
        layers = layers + 1;
        % Determine the next mesh size around each body
        next_r_bounds = {};
        next_z_bounds = {};
        next_drs = [];
        next_dzs = [];
        new_skipped_bodies = [];
        for jbody = 1:Nbody
            this_level = length(body_mesh_chain{jbody});
            this_L = L*seq_box_frac^(this_level);
            this_dr = drs(body_mesh_chain{jbody}(end))*seq_mesh_frac;
            this_dz = dzs(body_mesh_chain{jbody}(end))*seq_mesh_frac;

            if (this_L*0.5 <= stop_radii*body_radii(jbody)) ||...
               (min([this_dr, this_dz]) < min_dmesh)
                %disp(sprintf('skipping body %d',jbody)) %DEBUG
                new_skipped_bodies(end+1) = jbody;
            else

                next_drs(jbody) = this_dr;
                next_dzs(jbody) = this_dz;

                prev_id = body_mesh_chain{jbody}(end);
                prev_rs = rss{prev_id};
                prev_zs = zss{prev_id};

                this_r_bounds = [0 body_rs(jbody) + this_L*0.5];
                this_z_bounds = body_zs(jbody) + this_L*[-0.5 0.5];

                %[~, this_z_bounds] = nearest_index(prev_zs, this_z_bounds);
                %[~, this_r_bounds(1)] = nearest_index(prev_rs, this_r_bounds(1), 'below');
                [~, this_r_bounds(2)] = nearest_index(prev_rs, this_r_bounds(2), 'above');
                [~, this_z_bounds(1)] = nearest_index(prev_zs, this_z_bounds(1), 'below');
                [~, this_z_bounds(2)] = nearest_index(prev_zs, this_z_bounds(2), 'above');
                %this_r_bounds(1) = max(this_r_bounds(1), prev_rs(2));
                this_z_bounds(1) = max(this_z_bounds(1), prev_zs(2));
                this_r_bounds(2) = min(this_r_bounds(2), prev_rs(end-1));
                this_z_bounds(2) = min(this_z_bounds(2), prev_zs(end-1));

                rz_relative_bounds = [ this_r_bounds(2) - 0,...
                                       this_z_bounds(1) - body_z_coords(jbody),...
                                       this_z_bounds(2) - body_z_coords(jbody) ];

                if min(min( abs( rz_relative_bounds ) - body_radii(jbody))) < 0
                    new_skipped_bodies(end+1) = jbody;
                end

                next_r_bounds{jbody} = this_r_bounds;
                next_z_bounds{jbody} = this_z_bounds;
                [this_r_bounds; this_z_bounds]; % DEBUG
            end
        end

        for jskippedbody = new_skipped_bodies
            unskipped_bodies(unskipped_bodies == jskippedbody) = [];
        end

        if ~isempty(unskipped_bodies)
            % cluster overlapping meshes

            clusters = num2cell(unskipped_bodies);
            cluster_r_bounds = next_r_bounds(unskipped_bodies);
            cluster_z_bounds = next_z_bounds(unskipped_bodies);

            cluster_parent_drs = [];
            cluster_parent_dzs = [];
            for jbody = unskipped_bodies
                cluster_parent_drs(end+1) = drs(body_mesh_chain{jbody}(end));
                cluster_parent_dzs(end+1) = dzs(body_mesh_chain{jbody}(end));
            end

            didmerge = true;
            while (didmerge)
                didmerge = false;
                for jcluster = 1:length(clusters)
                    for icluster = jcluster+1:length(clusters)
                        if ~( (cluster_r_bounds{jcluster}(2) < cluster_r_bounds{icluster}(1) - max(cluster_parent_drs(jcluster), cluster_parent_drs(icluster))) ||...
                                (cluster_r_bounds{icluster}(2) < cluster_r_bounds{jcluster}(1) - max(cluster_parent_drs(jcluster), cluster_parent_drs(icluster))) ||...
                                (cluster_z_bounds{jcluster}(2) < cluster_z_bounds{icluster}(1) - max(cluster_parent_dzs(jcluster), cluster_parent_dzs(icluster))) ||...
                                (cluster_z_bounds{icluster}(2) < cluster_z_bounds{jcluster}(1) - max(cluster_parent_dzs(jcluster), cluster_parent_dzs(icluster))) )
                            %  there's an overlap => merge the clusters
                            %disp('MERGE:')
                            %disp(clusters{jcluster})
                            %disp(clusters{icluster})

                            clusters{jcluster} = [clusters{jcluster} clusters{icluster}];
                            %cluster_r_bounds{jcluster}(1) = min(cluster_r_bounds{jcluster}(1), cluster_r_bounds{icluster}(1));
                            cluster_z_bounds{jcluster}(1) = min(cluster_z_bounds{jcluster}(1), cluster_z_bounds{icluster}(1));
                            cluster_r_bounds{jcluster}(2) = max(cluster_r_bounds{jcluster}(2), cluster_r_bounds{icluster}(2));
                            cluster_z_bounds{jcluster}(2) = max(cluster_z_bounds{jcluster}(2), cluster_z_bounds{icluster}(2));
                            cluster_parent_drs(jcluster) = max(cluster_parent_drs(jcluster), cluster_parent_drs(icluster));
                            cluster_parent_dzs(jcluster) = max(cluster_parent_dzs(jcluster), cluster_parent_dzs(icluster));
                            clusters(icluster) = [];
                            cluster_r_bounds(icluster) = [];
                            cluster_z_bounds(icluster) = [];
                            cluster_parent_drs(icluster) = [];
                            cluster_parent_dzs(icluster) = [];
                            didmerge = true;
                            break

                        end
                       if didmerge; break; end
                    end
                end


            end

            % create a new mesh for each cluster
            for jcluster = 1:length(clusters)
                if length(clusters{jcluster}) > 0

                    this_r_bounds = [0, 0];
                    this_z_bounds = [Inf, -Inf];
                    this_dr = Inf;
                    this_dz = Inf;
                    for jbody = clusters{jcluster}
                        %this_r_bounds(1) = min(this_r_bounds(1), next_r_bounds{jbody}(1));
                        this_z_bounds(1) = min(this_z_bounds(1), next_z_bounds{jbody}(1));
                        this_r_bounds(2) = max(this_r_bounds(2), next_r_bounds{jbody}(2));
                        this_z_bounds(2) = max(this_z_bounds(2), next_z_bounds{jbody}(2));
                        this_dr = min(this_dr, next_drs(jbody));
                        this_dz = min(this_dz, next_dzs(jbody));
                    end

                    % Ensure corners match up with parent nodes (probably should not be necessary)
                    parent_meshi = 1;
                    for jbody = clusters{jcluster}
                        this_parent_meshi = body_mesh_chain{jbody}(end);
                        if ( drs(this_parent_meshi) < drs(parent_meshi) ) % TODO: check dzs?
                            parent_meshi = this_parent_meshi;
                        end
                    end
                    %disp(sprintf('This mesh: %d    parent mesh: %d', length(rss)+1, parent_meshi))
                    parent_rs = rss{parent_meshi};
                    parent_zs = zss{parent_meshi};

                    %[~, this_r_bounds(1)] = nearest_index(parent_rs, this_r_bounds(1), 'below');
                    [~, this_r_bounds(2)] = nearest_index(parent_rs, this_r_bounds(2), 'above');
                    [~, this_z_bounds(1)] = nearest_index(parent_zs, this_z_bounds(1), 'below');
                    [~, this_z_bounds(2)] = nearest_index(parent_zs, this_z_bounds(2), 'above');
                    %this_r_bounds(1) = max(this_r_bounds(1), parent_rs(2));
                    this_z_bounds(1) = max(this_z_bounds(1), parent_zs(2));
                    this_r_bounds(2) = min(this_r_bounds(2), parent_rs(end-1));
                    this_z_bounds(2) = min(this_z_bounds(2), parent_zs(end-1));

                    % Construct mesh
                    this_rs = [this_r_bounds(1):this_dr:this_r_bounds(2)];
                    this_zs = [this_z_bounds(1):this_dz:this_z_bounds(2)];

                    % this sometimes skips the last point due to minimal errors
                    if (this_r_bounds(2) - this_rs(end)) > this_dr*0.1
                        assert(abs(this_r_bounds(2) - this_rs(end) - this_dr) < this_dr*1E-5);
                        this_rs(end+1) = this_r_bounds(2);
                    end
                    if (this_z_bounds(2) - this_zs(end)) > this_dz*0.1
                        assert(abs(this_z_bounds(2) - this_zs(end) - this_dz) < this_dz*1E-5);
                        this_zs(end+1) = this_z_bounds(2);
                    end

                    % only proceed if the mesh has at least two nodes per side
                    if ((length(this_rs) > 1) && (length(this_zs) > 1))
                        rss{end+1} = this_rs;
                        zss{end+1} = this_zs;
                        drs(end+1) = this_dr;
                        dzs(end+1) = this_dz;

                        this_meshi = length(rss);

                        for jbody = clusters{jcluster}
                            body_mesh_chain{jbody}(end+1) = this_meshi;
                        end
                    else
                        for jbody = clusters{jcluster}
                            % add to skipped bodies
                            unskipped_bodies(unskipped_bodies == jbody) = [];
                        end
                    end

                end
            end
        end
    end

    meshes = {};
    for j = 1:length(rss)
        meshes{j} = {rss{j}, zss{j}};
    end
