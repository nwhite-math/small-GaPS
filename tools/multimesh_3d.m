% MULTIMESH_3D Construct 3D multi-resolution mesh and derivative operators
%
% [idgrids, nodedata, dx, d2x, dy, d2y, dz, d2z, int_weight] = MULTIMESH_3D(MESHES)
% [...] = MULTIMESH_3D(MESHES, ..., 'alignmentTol', ALIGNMENTTOL)
% [...] = MULTIMESH_3D(MESHES, ..., 'edgeOrder', EDGEORDER)
%
%
% Input:
%
%  MESHES must be a cell-list of cell-tuples; each cell-tuple containing a list
%  of x values, a list of y values, and a list of z values to specify a mesh.
% For example, a valid input would be
%  >> MESHES = { {[0 1 2], [0 1 2 3], [0 5 10]},...
%                {[1 1.2 1.3 2], [0 0.1 0.5 1 1.5 2], [5 6 7 8 9 10]} }
%
%  Note that the corners of each submesh must align with nodes of a parent mesh.
%  This alignment is enforced up to ALIGNMENTTOL, default 1E-12.
%  The ordering of the input meshes does not matter. Meshes may be arbitrarily
%  nested, but they may not overlap without being fully contained by a parent.
%
%  Derivative matrices are computed to order 2 (i.e. central difference order).
%  The optional parameter EDGEORDER (default 2) may be 1 or 2, specifying the
%  derivative order at the outermost mesh edges.
%
%
% Output:
%
%  IDGRIDS is a cell containing one array for each input mesh. The entries of
%  each array are the ids of the corresponding nodes.
%
%  NODEDATA is a struct containing the following data for each node:
%   'grid': the main (parent) grid to which the node belongs
%   'subgrid': the sub (child) grid to which the node belongs, if applicable
%   'x': x value of node
%   'y': y value of node
%   'z': z value of node
%   'xi': x index in the main grid
%   'yi': y index in the main grid
%   'zi': z index in the main grid
%   'xi_sub': x index in the subgrid, if applicable
%   'yi_sub': y index in the subgrid, if applicable
%   'zi_sub': z index in the subgrid, if applicable
%   'type': NORMAL, SHARED,
%           INT_X_FACE, INT_Y_FACE, INT_Z_FACE,
%           INT_XY_EDGE, INT_YZ_EDGE, INT_ZX_EDGE,
%           EXT_X_FACE, EXT_Y_FACE, EXT_Z_FACE,
%           EXT_XY_EDGE, EXT_YZ_EDGE, EXT_ZX_EDGE,
%           EXT_CORNER
%
%  DX, D2X, DY, D2Y, DZ, and D2Z are differentiation matrices that act on
%  vectors indexed by node number.
%
%  INT_WEIGHT is a vector giving a 3D trapezoidal integration weight for each
%  node.
%
%
% Example:
%
%  Suppose we had a 4x3x3 coarse grid (marked 'O' below), with a 3x3x3 fine grid
%  (marked by *, with shared nodes marked by X). In this diagram, only two
%  z layers of the coarse grid are shown.
%
%  O     X  *  X     O   |         *  *  *         |   O     X  *  X     0
%        *  *  *         |         *  *  *         |         *  *  *
%  O     X  *  X     O   |         *  *  *         |   O     X  *  X     0
%                        |                         |
%  O     O     O     O   |                         |   O     O     O     O
%
%  The coarse grid might have x, y, and z vectors [0:3], [0:2], and [0:2], while
%  the fine grid might be [1 1.5 2], [1 1.5 2], [0 0.5 1]. Then, to construct
%  differentiation matrices, we could call
%
%    >> coarse_xs = [0 1 2 3];
%    >> coarse_ys = [0 1 2];
%    >> coarse_zs = [0 1 2];
%    >> fine_xs = [1 1.5 2];
%    >> fine_ys = [1 1.5 2];
%    >> fine_zs = [0 0.5 1];
%    >> [idgrids, nodedata, dx, d2x, dy, d2y, dz, d2z, int_weight] = ...
%          MULTIMESH_3D( { {coarse_xs, coarse_ys, coarse_zs}, ...
%                          {fine_xs, fine_ys, fine_zs} } );
%

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-04-15
% Updated: 2019-06-30

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


function [idgrids, nodedata, dx, d2x, dy, d2y, dz, d2z, int_weight] = multimesh_3d(varargin)

    debug = false;
    debug_verbose = false;

    % Define node types
    NORMAL      = 0;
    SHARED      = 1;
    INT_X_FACE  = 2;
    INT_Y_FACE  = 3;
    INT_Z_FACE  = 4;
    INT_XY_EDGE = 5;
    INT_YZ_EDGE = 6;
    INT_ZX_EDGE = 7;
    EXT_X_FACE  = 8;
    EXT_Y_FACE  = 9;
    EXT_Z_FACE  = 10;
    EXT_XY_EDGE = 11;
    EXT_YZ_EDGE = 12;
    EXT_ZX_EDGE = 13;
    EXT_CORNER  = 14;

    parser = inputParser;
    parser.addRequired('meshes', @(x) iscell(x));
    parser.addParameter('edgeOrder', 2, @(x) (x == 1 || x == 2));
    parser.addParameter('progressBar', false, @(x) isscalar(x));
    parser.addParameter('alignmentTol', 1E-12, @(x) (isscalar(x) && isnumeric(x) && x >= 0));
    parse(parser, varargin{:});
    meshes = parser.Results.meshes;

    edge_order = parser.Results.edgeOrder;

    show_waitbar= parser.Results.progressBar;

    alignment_tol = parser.Results.alignmentTol;

    % Number of meshes
    meshN = length(meshes);
    assert(meshN >= 1, 'At least one input mesh is required.');

    xss = {};
    yss = {};
    zss = {};
    xmaxs = [];
    xmins = [];
    ymaxs = [];
    ymins = [];
    zmaxs = [];
    zmins = [];

    % Extract x, y and z from each mesh
    for j = 1:meshN
        jmesh = meshes{j};
        assert(iscell(jmesh), 'Each mesh must be passed as a cell of form {[x...], [y...], [z...]}.');
        assert(length(jmesh)==3, 'Each mesh must be passed as a cell of form {[x...], [y...], [z...]}.');
        xss{j} = jmesh{1};
        yss{j} = jmesh{2};
        zss{j} = jmesh{3};
        assert(min([length(xss{j}) length(yss{j}) length(zss{j})]) > 1, 'Each mesh must contain at least two nodes per side.');
        xmaxs(j) = max(xss{j});
        xmins(j) = min(xss{j});
        ymaxs(j) = max(yss{j});
        ymins(j) = min(yss{j});
        zmaxs(j) = max(zss{j});
        zmins(j) = min(zss{j});
    end
    xmin_all = min(xmins);
    xmax_all = max(xmaxs);
    ymin_all = min(ymins);
    ymax_all = max(ymaxs);
    zmin_all = min(zmins);
    zmax_all = max(zmaxs);

    % Determine the dependency tree of the meshes
    % parentlist(j) gives the immediate parent of mesh j
    % ancestormat(j,k) returns 1 if k is an ancestor of j
    if (debug); disp('Constructing dependency tree..'); end
    parentlist = zeros(1,meshN);
    ancestormat = zeros(meshN,meshN);
    for j = 1:meshN
        ancestormat(j,j) = 1;
        for k = 1:meshN
            if (k ~= j)
                if ( (xmaxs(j) <= xmaxs(k)) && (xmins(j) >= xmins(k)) &&...
                     (ymaxs(j) <= ymaxs(k)) && (ymins(j) >= ymins(k)) &&...
                     (zmaxs(j) <= zmaxs(k)) && (zmins(j) >= zmins(k)) )
                    ancestormat(j,k) = 1;
                    if parentlist(j) == 0
                        parentlist(j) = k;
                    else
                        iparent = parentlist(j);
                        % mesh has multiple parents. Check which is smaller.
                        if ( (xmaxs(k) <= xmaxs(iparent)) && (xmins(k) >= xmins(iparent)) &&...
                             (ymaxs(k) <= ymaxs(iparent)) && (ymins(k) >= ymins(iparent)) &&...
                             (zmaxs(k) <= zmaxs(iparent)) && (zmins(k) >= zmins(iparent)) )
                            parentlist(j) = k;
                        end
                    end
                end
            end
        end
    end

    % Make sure there is only one root, and no overlaps of non-parent-child
    assert(sum(parentlist==0) == 1, 'Only one root mesh allowed (i.e., one mesh must contain all others).');
    root_gridindex = find(parentlist==0);
    for j = 1:meshN
        for k = 1:meshN
            if (ancestormat(j,k) == 0) && (ancestormat(k,j) == 0)
                % ensure no overlap
                if ~( (xmaxs(j) < xmins(k)) || (xmaxs(k) < xmins(j)) ||...
                      (ymaxs(j) < ymins(k)) || (ymaxs(k) < ymins(j)) ||...
                      (zmaxs(j) < zmins(k)) || (zmaxs(k) < zmins(j)) )
                    error('Submeshes may not overlap');
                end
            end
        end
    end

    % Ensure that the corners of each mesh line up with nodes of its parent mesh
    if (debug); disp('Checking alignment...'); end
    xmin_parent_indexs = zeros(1,meshN);
    xmax_parent_indexs = zeros(1,meshN);
    ymin_parent_indexs = zeros(1,meshN);
    ymax_parent_indexs = zeros(1,meshN);
    zmin_parent_indexs = zeros(1,meshN);
    zmax_parent_indexs = zeros(1,meshN);

    for j = 1:meshN

        iparent = parentlist(j);
        if iparent ~= 0 % if not the root mesh
            xs_child = xss{j};
            ys_child = yss{j};
            zs_child = zss{j};
            xs_parent = xss{iparent};
            ys_parent = yss{iparent};
            zs_parent = zss{iparent};

            [xmin_index, xmincheck] = nearest_index(xs_parent, xmins(j));
            [xmax_index, xmaxcheck] = nearest_index(xs_parent, xmaxs(j));
            [ymin_index, ymincheck] = nearest_index(ys_parent, ymins(j));
            [ymax_index, ymaxcheck] = nearest_index(ys_parent, ymaxs(j));
            [zmin_index, zmincheck] = nearest_index(zs_parent, zmins(j));
            [zmax_index, zmaxcheck] = nearest_index(zs_parent, zmaxs(j));
            xmin_parent_indexs(j) = xmin_index;
            xmax_parent_indexs(j) = xmax_index;
            ymin_parent_indexs(j) = ymin_index;
            ymax_parent_indexs(j) = ymax_index;
            zmin_parent_indexs(j) = zmin_index;
            zmax_parent_indexs(j) = zmax_index;

            % Ensure the corners of each mesh line up with nodes of the parent mesh
            errmsg_submesh_badalign = 'Each submesh must have corners aligning exactly with nodes of its parent mesh.';
            assert(abs(xmincheck - xmins(j)) < alignment_tol, errmsg_submesh_badalign)
            assert(abs(xmaxcheck - xmaxs(j)) < alignment_tol, errmsg_submesh_badalign)
            assert(abs(ymincheck - ymins(j)) < alignment_tol, errmsg_submesh_badalign)
            assert(abs(ymaxcheck - ymaxs(j)) < alignment_tol, errmsg_submesh_badalign)
            assert(abs(zmincheck - zmins(j)) < alignment_tol, errmsg_submesh_badalign)
            assert(abs(zmaxcheck - zmaxs(j)) < alignment_tol, errmsg_submesh_badalign)

            % Ensure that the mesh is not zero width or height
            errmsg_submesh_onespace = 'Each submesh must have nonzero width and height.';
            assert(max(xs_child) > min(xs_child), errmsg_submesh_onespace);
            assert(max(ys_child) > min(ys_child), errmsg_submesh_onespace);
            assert(max(zs_child) > min(zs_child), errmsg_submesh_onespace);
            assert(xmin_index < xmax_index, errmsg_submesh_onespace);
            assert(ymin_index < ymax_index, errmsg_submesh_onespace);
            assert(zmin_index < zmax_index, errmsg_submesh_onespace);

            % Ensure that the mesh spacing is above the tolerance
            errmsg_meshspacetol = sprintf('Mesh spacing must be greater than alignment tolerance (%g). Or, change the alignment tolerance.', alignment_tol) ;
            assert(min(abs(diff(xs_child)))>alignment_tol, errmsg_meshspacetol);
            assert(min(abs(diff(ys_child)))>alignment_tol, errmsg_meshspacetol);
            assert(min(abs(diff(zs_child)))>alignment_tol, errmsg_meshspacetol);

            % Ensure that the mesh has at least 3 nodes in each direction
            errmsg_submesh_threespace = 'Each submesh must have at least three nodes in width and height.';
            assert(length(xs_child) >= 3, errmsg_submesh_threespace);
            assert(length(ys_child) >= 3, errmsg_submesh_threespace);
            assert(length(zs_child) >= 3, errmsg_submesh_threespace);

            % Ensure that the submesh spacing divides the parent mesh spacing
            errmsg_submesh_divisible = 'The submesh spacing must divide the parent mesh spacing.';
            assert( ((xmax_index - xmin_index) <= length(xs_child)-1), errmsg_submesh_divisible );
            assert( ((ymax_index - ymin_index) <= length(ys_child)-1), errmsg_submesh_divisible );
            assert( ((zmax_index - zmin_index) <= length(zs_child)-1), errmsg_submesh_divisible );
            for xj = xmin_index+1:xmax_index-1
                [xcheck_index, xcheck] = nearest_index(xs_child, xs_parent(xj));
                assert(abs(xcheck - xs_child(xcheck_index)) < alignment_tol, errmsg_submesh_divisible)
            end
            for yj = ymin_index+1:ymax_index-1
                [ycheck_index, ycheck] = nearest_index(ys_child, ys_parent(yj));
                assert(abs(ycheck - ys_child(ycheck_index)) < alignment_tol, errmsg_submesh_divisible)
            end
            for zj = zmin_index+1:zmax_index-1
                [zcheck_index, zcheck] = nearest_index(zs_child, zs_parent(zj));
                assert(abs(zcheck - zs_child(zcheck_index)) < alignment_tol, errmsg_submesh_divisible)
            end

            % Ensure that the submesh is at least one full space inside the parent mesh
            errmsg_submesh_padding = 'Each submesh must be one full mesh space inside the edges of its parent mesh, unless it is at an edge.';
            assert( (xmin_index > 1) || (xs_child(1) == xmin_all), errmsg_submesh_padding );
            assert( (xmax_index < length(xss{iparent})) || (xs_child(end) == xmax_all), errmsg_submesh_padding);
            assert( (ymin_index > 1) || (ys_child(1) == ymin_all), errmsg_submesh_padding );
            assert( (ymax_index < length(yss{iparent})) || (ys_child(end) == ymax_all), errmsg_submesh_padding);
            assert( (zmin_index > 1) || (zs_child(1) == zmin_all), errmsg_submesh_padding );
            assert( (zmax_index < length(zss{iparent})) || (zs_child(end) == zmax_all), errmsg_submesh_padding);
        end
    end

    if (debug); disp('Constructing nodedata...'); end
    % Assign id numbers to each of the nodes. We start at the largest mesh and work ydownwards; shared nodes will take the value of their parent mesh, not the child.
    sort_mesh_i = sorttree(parentlist);

    % Nodedata contains data for each node. Nodedata has the following structure, each of which is indexed by node.
    %   'grid': the main (parent) grid to which the node belongs
    %   'subgrid': the sub (child) grid to which the node belongs, if applicable
    %   'x': x value of node
    %   'y': y value of node
    %   'z': z value of node
    %   'xi': x index in the main grid
    %   'yi': y index in the main grid
    %   'zi': z index in the main grid
    %   'xi_sub': x index in the subgrid, if applicable
    %   'yi_sub': y index in the subgrid, if applicable
    %   'zi_sub': z index in the subgrid, if applicable
    %   'type': NORMAL, SHARED,
    %           INT_X_FACE, INT_Y_FACE, INT_Z_FACE,
    %           INT_XY_EDGE, INT_YZ_EDGE, INT_ZX_EDGE,
    %           EXT_X_FACE, EXT_Y_FACE, EXT_Z_FACE,
    %           EXT_XY_EDGE, EXT_YZ_EDGE, EXT_ZX_EDGE,
    %           EXT_CORNER

    nodedata.grid = [];
    nodedata.subgrid = [];
    nodedata.x = [];
    nodedata.y = [];
    nodedata.z = [];
    nodedata.xi = [];
    nodedata.yi = [];
    nodedata.zi = [];
    nodedata.xi_sub = [];
    nodedata.yi_sub = [];
    nodedata.zi_sub = [];
    nodedata.type = [];

    cur_nodeid = 0;
    idgrids = {};
    for junsorted = 1:meshN
        j = sort_mesh_i(junsorted);

        iparent = parentlist(j);
        xs_self = xss{j};
        ys_self = yss{j};
        zs_self = zss{j};

        idgrid = zeros(length(xs_self), length(ys_self), length(zs_self));

        for zi = 1:length(zs_self)
            for yi = 1:length(ys_self)
                for xi = 1:length(xs_self)
                    x = xs_self(xi);
                    y = ys_self(yi);
                    z = zs_self(zi);

                    % determine if node is shared by a parent
                    node_is_shared = false;
                    if (iparent ~= 0)
                        xs_parent = xss{iparent};
                        ys_parent = yss{iparent};
                        zs_parent = zss{iparent};
                        [xi_parent, x_parent] = nearest_index(xs_parent, x);
                        [yi_parent, y_parent] = nearest_index(ys_parent, y);
                        [zi_parent, z_parent] = nearest_index(zs_parent, z);
                        if ( (abs(x_parent - x) < alignment_tol) &&...
                             (abs(y_parent - y) < alignment_tol) &&...
                             (abs(z_parent - z) < alignment_tol) )
                             node_is_shared = true;
                             idgrid_parent = idgrids{iparent};
                             parent_nodeid = idgrid_parent(xi_parent, yi_parent, zi_parent);
                        end
                    end


                    % if node is not shared, then a new entry is needed.
                    if ~node_is_shared
                        cur_nodeid = cur_nodeid + 1;
                        idgrid(xi, yi, zi) = cur_nodeid;
                        nodedata.grid(cur_nodeid) = j;
                        nodedata.subgrid(cur_nodeid) = 0;
                        nodedata.x(cur_nodeid) = x;
                        nodedata.y(cur_nodeid) = y;
                        nodedata.z(cur_nodeid) = z;
                        nodedata.xi(cur_nodeid) = xi;
                        nodedata.yi(cur_nodeid) = yi;
                        nodedata.zi(cur_nodeid) = zi;
                        nodedata.xi_sub(cur_nodeid) = 0;
                        nodedata.yi_sub(cur_nodeid) = 0;
                        nodedata.zi_sub(cur_nodeid) = 0;
                        if ( (x == xmax_all) || (x == xmin_all) ) &&...
                            ( (y == ymax_all) || (y == ymin_all) ) &&...
                            ( (z == zmax_all) || (z == zmin_all) )
                            nodedata.type(cur_nodeid) = EXT_CORNER;
                        elseif ( (x == xmax_all) || (x == xmin_all) ) &&...
                               ( (y == ymax_all) || (y == ymin_all) )
                               nodedata.type(cur_nodeid) = EXT_XY_EDGE;
                        elseif ( (y == ymax_all) || (y == ymin_all) ) &&...
                               ( (z == zmax_all) || (z == zmin_all) )
                               nodedata.type(cur_nodeid) = EXT_YZ_EDGE;
                        elseif ( (z == zmax_all) || (z == zmin_all) ) &&...
                               ( (x == xmax_all) || (x == xmin_all) )
                               nodedata.type(cur_nodeid) = EXT_ZX_EDGE;
                        elseif ( (x == xmax_all) || (x == xmin_all) )
                            nodedata.type(cur_nodeid) = EXT_X_FACE;
                        elseif ( (y == ymax_all) || (y == ymin_all) )
                            nodedata.type(cur_nodeid) = EXT_Y_FACE;
                        elseif ( (z == zmax_all) || (z == zmin_all) )
                            nodedata.type(cur_nodeid) = EXT_Z_FACE;
                        elseif ( (xi == length(xs_self)) || (xi == 1) ) &&...
                               ( (yi == length(ys_self)) || (yi == 1) ) &&...
                               ( (zi == length(zs_self)) || (zi == 1) )
                               % "Interior corner"
                               assert(false); % This should never occur, since it should always be a shared node.
                        elseif ( (xi == length(xs_self)) || (xi == 1) ) &&...
                               ( (yi == length(ys_self)) || (yi == 1) )
                               nodedata.type(cur_nodeid) = INT_XY_EDGE;
                        elseif ( (yi == length(ys_self)) || (yi == 1) ) &&...
                               ( (zi == length(zs_self)) || (zi == 1) )
                               nodedata.type(cur_nodeid) = INT_YZ_EDGE;
                        elseif ( (zi == length(zs_self)) || (zi == 1) ) &&...
                               ( (xi == length(xs_self)) || (xi == 1) )
                               nodedata.type(cur_nodeid) = INT_ZX_EDGE;
                        elseif ( (xi == length(xs_self)) || (xi == 1) )
                            nodedata.type(cur_nodeid) = INT_X_FACE;
                        elseif ( (yi == length(ys_self)) || (yi == 1))
                            nodedata.type(cur_nodeid) = INT_Y_FACE;
                        elseif ( (zi == length(zs_self)) || (zi == 1))
                            nodedata.type(cur_nodeid) = INT_Z_FACE;
                        else
                            nodedata.type(cur_nodeid) = NORMAL;
                        end
                    else
                        % if node is shared, then it has already been listed, but we
                        % still must repopulate it with the correct parent information.
                        assert(iparent ~= 0)
                        idgrid(xi, yi, zi) = parent_nodeid;
                        nodedata.grid(parent_nodeid) = iparent;
                        nodedata.x(parent_nodeid) = x_parent;
                        nodedata.y(parent_nodeid) = y_parent;
                        nodedata.z(parent_nodeid) = z_parent;
                        nodedata.xi(parent_nodeid) = xi_parent;
                        nodedata.yi(parent_nodeid) = yi_parent;
                        nodedata.zi(parent_nodeid) = zi_parent;
                        nodedata.subgrid(parent_nodeid) = j;
                        nodedata.xi_sub(parent_nodeid) = xi;
                        nodedata.yi_sub(parent_nodeid) = yi;
                        nodedata.zi_sub(parent_nodeid) = zi;
                        nodedata.type(parent_nodeid) = SHARED;

                        if ( (x == xmax_all) || (x == xmin_all) ) &&...
                            ( (y == ymax_all) || (y == ymin_all) ) &&...
                            ( (z == zmax_all) || (z == zmin_all) )
                            nodedata.type(parent_nodeid) = EXT_CORNER;
                        elseif ( (x == xmax_all) || (x == xmin_all) ) &&...
                               ( (y == ymax_all) || (y == ymin_all) )
                               nodedata.type(parent_nodeid) = EXT_XY_EDGE;
                        elseif ( (y == ymax_all) || (y == ymin_all) ) &&...
                               ( (z == zmax_all) || (z == zmin_all) )
                               nodedata.type(parent_nodeid) = EXT_YZ_EDGE;
                        elseif ( (z == zmax_all) || (z == zmin_all) ) &&...
                               ( (x == xmax_all) || (x == xmin_all) )
                               nodedata.type(parent_nodeid) = EXT_ZX_EDGE;
                        elseif ( (x == xmax_all) || (x == xmin_all) )
                            nodedata.type(parent_nodeid) = EXT_X_FACE;
                        elseif ( (y == ymax_all) || (y == ymin_all) )
                            nodedata.type(parent_nodeid) = EXT_Y_FACE;
                        elseif ( (z == zmax_all) || (z == zmin_all) )
                            nodedata.type(parent_nodeid) = EXT_Z_FACE;
                        end
                    end


                end
            end
        end

        idgrids{j} = idgrid;
    end
    nodeN = length(nodedata.x);
    assert(nodeN==cur_nodeid);

    % Set up shortcuts to get x, y, or z size of idgrid
    xlength = @(idgrid) size(idgrid, 1);
    ylength = @(idgrid) size(idgrid, 2);
    zlength = @(idgrid) size(idgrid, 3);

    % Set up neighbor fetching shortcuts
    get_xleft_node = @(node) get_neighbor_node(node, 'xleft', nodedata, idgrids);
    get_xright_node = @(node) get_neighbor_node(node, 'xright', nodedata, idgrids);
    get_ydown_node = @(node) get_neighbor_node(node, 'ydown', nodedata, idgrids);
    get_yup_node = @(node) get_neighbor_node(node, 'yup', nodedata, idgrids);
    get_zback_node = @(node) get_neighbor_node(node, 'zback', nodedata, idgrids);
    get_zfore_node = @(node) get_neighbor_node(node, 'zfore', nodedata, idgrids);

    if (debug);

        if (debug_verbose)
            nodedata
            for j = 1:meshN
                idgrids{j}
            end
        end


        disp('Checking neighbor integrity...');
        % check neighbor code
        for this_node = 1:nodeN
            this_xval = nodedata.x(this_node);
            this_yval = nodedata.y(this_node);
            this_zval = nodedata.z(this_node);

            this_xleft_node = get_xleft_node(this_node);
            this_xleft_xval = NaN;
            this_xleft_yval = NaN;
            this_xleft_zval = NaN;
            if (this_xleft_node > 0)
                this_xleft_xval = nodedata.x(this_xleft_node);
                this_xleft_yval = nodedata.y(this_xleft_node);
                this_xleft_zval = nodedata.z(this_xleft_node);
                assert(this_xleft_xval < this_xval);
                assert(abs(this_xleft_yval - this_yval) < alignment_tol);
                assert(abs(this_xleft_zval - this_zval) < alignment_tol);
            end

            this_xright_node = get_xright_node(this_node);
            this_xright_xval = NaN;
            this_xright_yval = NaN;
            this_xright_zval = NaN;
            if (this_xright_node > 0)
                this_xright_xval = nodedata.x(this_xright_node);
                this_xright_yval = nodedata.y(this_xright_node);
                this_xright_zval = nodedata.z(this_xright_node);
                assert(this_xright_xval > this_xval);
                assert(abs(this_xright_yval - this_yval) < alignment_tol);
                assert(abs(this_xright_zval - this_zval) < alignment_tol);
            end

            this_ydown_node = get_ydown_node(this_node);
            this_ydown_xval = NaN;
            this_ydown_yval = NaN;
            this_ydown_zval = NaN;
            if (this_ydown_node > 0)
                this_ydown_xval = nodedata.x(this_ydown_node);
                this_ydown_yval = nodedata.y(this_ydown_node);
                this_ydown_zval = nodedata.z(this_ydown_node);
                assert(abs(this_ydown_xval - this_xval) < alignment_tol);
                assert(this_ydown_yval < this_yval);
                assert(abs(this_ydown_zval - this_zval) < alignment_tol);
            end

            this_yup_node = get_yup_node(this_node);
            this_yup_xval = NaN;
            this_yup_yval = NaN;
            this_yup_zval = NaN;
            if (this_yup_node > 0)
                this_yup_xval = nodedata.x(this_yup_node);
                this_yup_yval = nodedata.y(this_yup_node);
                this_yup_zval = nodedata.z(this_yup_node);
                assert(abs(this_yup_xval - this_xval) < alignment_tol);
                assert(this_yup_yval > this_yval);
                assert(abs(this_yup_zval - this_zval) < alignment_tol);
            end

            this_zback_node = get_zback_node(this_node);
            this_zback_xval = NaN;
            this_zback_yval = NaN;
            this_zback_zval = NaN;
            if (this_zback_node > 0)
                this_zback_xval = nodedata.x(this_zback_node);
                this_zback_yval = nodedata.y(this_zback_node);
                this_zback_zval = nodedata.z(this_zback_node);
                assert(abs(this_zback_xval - this_xval) < alignment_tol);
                assert(abs(this_zback_yval - this_yval) < alignment_tol);
                assert(this_zback_zval < this_zval);
            end

            this_zfore_node = get_zfore_node(this_node);
            this_zfore_xval = NaN;
            this_zfore_yval = NaN;
            this_zfore_zval = NaN;
            if (this_zfore_node > 0)
                this_zfore_xval = nodedata.x(this_zfore_node);
                this_zfore_yval = nodedata.y(this_zfore_node);
                this_zfore_zval = nodedata.z(this_zfore_node);
                assert(abs(this_zfore_xval - this_xval) < alignment_tol);
                assert(abs(this_zfore_yval - this_yval) < alignment_tol);
                assert(this_zfore_zval > this_zval);
            end

        end
    end



    %%%%%%%%%% SET UP DIFFERENTIATION MATRICES %%%%%%%%%%

    % make sparse zero matrices to start with.
    dx = sparse(nodeN, nodeN, 0);
    d2x = sparse(nodeN, nodeN, 0);
    dy = sparse(nodeN, nodeN, 0);
    d2y = sparse(nodeN, nodeN, 0);
    dz = sparse(nodeN, nodeN, 0);
    d2z = sparse(nodeN, nodeN, 0);

    % First, efficiently set up differences on all the interior "NORMAL" nodes.
    % Next, go back and deal with the edge cases one by one.

    all_nodes = [1:nodeN];
    nontrivial_nodes_dx = all_nodes( ~ismember(nodedata.type, [NORMAL, INT_Y_FACE, INT_Z_FACE, INT_YZ_EDGE]) & ~( ismember(nodedata.type, [EXT_Y_FACE, EXT_Z_FACE, EXT_YZ_EDGE]) & (nodedata.subgrid == 0)) );
    nontrivial_nodes_dy = all_nodes( ~ismember(nodedata.type, [NORMAL, INT_Z_FACE, INT_X_FACE, INT_ZX_EDGE]) & ~( ismember(nodedata.type, [EXT_Z_FACE, EXT_X_FACE, EXT_ZX_EDGE]) & (nodedata.subgrid == 0)) );
    nontrivial_nodes_dz = all_nodes( ~ismember(nodedata.type, [NORMAL, INT_X_FACE, INT_Y_FACE, INT_XY_EDGE]) & ~( ismember(nodedata.type, [EXT_X_FACE, EXT_Y_FACE, EXT_XY_EDGE]) & (nodedata.subgrid == 0)) );

    % construct dx, dy, and dz for each mesh with deriv_findiff
    if (debug); disp('Constructing trivial dx, dy, and dz entries...'); end
    for j = 1:meshN
        this_xs = xss{j};
        this_ys = yss{j};
        this_zs = zss{j};
        this_idgrid = idgrids{j};

        [this_dx, this_d2x, ~, ~] = deriv_findiff(this_xs, 2);
        [this_dy, this_d2y, ~, ~] = deriv_findiff(this_ys, 2);
        [this_dz, this_d2z, ~, ~] = deriv_findiff(this_zs, 2);

        % permute results to match node ids...
        % TODO: is there a more efficient way to copy this data??
        for ky = 1:length(this_ys)
            for kz = 1:length(this_zs)
                dx(this_idgrid(:,ky,kz), this_idgrid(:,ky,kz)) = this_dx;
                d2x(this_idgrid(:,ky,kz), this_idgrid(:,ky,kz)) = this_d2x;
            end
        end
        for kz = 1:length(this_zs)
            for kx = 1:length(this_xs)
                dy(this_idgrid(kx,:,kz), this_idgrid(kx,:,kz)) = this_dy;
                d2y(this_idgrid(kx,:,kz), this_idgrid(kx,:,kz)) = this_d2y;
            end
        end
        for kx = 1:length(this_xs)
            for ky = 1:length(this_ys)
                dz(this_idgrid(kx,ky,:), this_idgrid(kx,ky,:)) = this_dz;
                d2z(this_idgrid(kx,ky,:), this_idgrid(kx,ky,:)) = this_d2z;
            end
        end

    end


    % Construct dx and d2x edge cases
    if (debug); disp('Constructing dx & d2x edge cases...'); end
    if (show_waitbar)
        f_wait = waitbar(0, sprintf('Constructing dx & d2x operators...'));
    end

    progress = 0;
    for this_node = nontrivial_nodes_dx
        progress = progress + 1;
        this_type         = nodedata.type(this_node);
        this_gridindex    = nodedata.grid(this_node);
        this_idgrid       = idgrids{this_gridindex};
        this_xi           = nodedata.xi(this_node);
        this_yi           = nodedata.yi(this_node);
        this_zi           = nodedata.zi(this_node);
        this_subgridindex = nodedata.subgrid(this_node);
        this_xi_sub       = nodedata.xi_sub(this_node);
        this_yi_sub       = nodedata.yi_sub(this_node);
        this_zi_sub       = nodedata.zi_sub(this_node);
        this_x            = nodedata.x(this_node);
        this_y            = nodedata.y(this_node);
        this_z            = nodedata.z(this_node);

        if (this_subgridindex == 0)
            this_subidgrid = [];
        else
            this_subidgrid = idgrids{this_subgridindex};
        end

        % clear this row of dx / d2x
        dx(this_node, :) = 0;
        d2x(this_node, :) = 0;

        if ismember(this_type, [NORMAL, SHARED, INT_Y_FACE, INT_Z_FACE, INT_YZ_EDGE, EXT_Y_FACE, EXT_Z_FACE, EXT_YZ_EDGE])
            % Easiest case, just do central difference

            xleft_node = get_xleft_node(this_node);
            xright_node = get_xright_node(this_node);
            assert((xleft_node ~= 0) && (xright_node ~= 0));
            xleft_x = nodedata.x(xleft_node);
            xright_x = nodedata.x(xright_node);

            h3 = (xright_x - this_x);
            h4 = (this_x - xleft_x);
            dx(this_node, xleft_node) = -(h3/h4)/(h3 + h4);
            dx(this_node, this_node) = 1/h4 - 1/h3;
            dx(this_node, xright_node) = (h4/h3)/(h3 + h4);

            d2x(this_node, xleft_node) = 2/(h4*(h3 + h4));
            d2x(this_node, this_node) = -2/(h3*h4);
            d2x(this_node, xright_node) = 2/(h3*(h3 + h4));
        elseif ismember(this_type, [EXT_X_FACE, EXT_XY_EDGE, EXT_ZX_EDGE, EXT_CORNER])
            % We're at the true edge of the whole mesh

            % use first or second order difference, depending on edge_order
            if (this_xi == 1)
                next_node = get_xright_node(this_node);
                far_node = get_xright_node(next_node);
            elseif (this_xi == xlength(this_idgrid));
                next_node = get_xleft_node(this_node);
                far_node = get_xleft_node(next_node);
            else
                assert(false)
            end
            assert((next_node ~= 0) && (far_node ~= 0));
            h3 = nodedata.x(next_node) - nodedata.x(this_node);
            h5 = nodedata.x(far_node) - nodedata.x(next_node);

            if (edge_order == 1)
                dx(this_node, this_node) = -1/h3;
                dx(this_node, next_node) = 1/h3;

                d2x(this_node, this_node) = 2/(h3*(h3+h5));
                d2x(this_node, next_node) = -2/(h3*h5);
                d2x(this_node, far_node) = 2/(h5*(h3+h5));
            else
                dx(this_node, this_node) = -(2*h3+h5)/(h3*(h3+h5));
                dx(this_node, next_node) = (1/h3 + 1/h5);
                dx(this_node, far_node) = -h3/(h5*(h3+h5));

                d2x(this_node, this_node) = 2/(h3*(h3+h5));
                d2x(this_node, next_node) = -2/(h3*h5);
                d2x(this_node, far_node) = 2/(h5*(h3+h5));
            end
        elseif ismember(this_type, [INT_XY_EDGE, INT_ZX_EDGE, INT_X_FACE])
            % We're at a mesh boundary within the system; we use 10 nearby points to construct
            % the correct second order difference (time for the ghost nodes).

            assert(this_subgridindex == 0);
            assert(this_gridindex ~= 0);

            % 1) get nearest neighbors
            bounding_nodes = get_bounding_shared_nodes(this_node, nodedata, idgrids, parentlist, xss, yss, zss, alignment_tol);
            assert(all(all(bounding_nodes(1,:,:) == bounding_nodes(2,:,:))));
            % get rid of x data in bounding_nodes
            bounding_nodes = permute(bounding_nodes, [2 3 1]);
            bounding_nodes = bounding_nodes(:,:,1);

            % 2) get side neighbors
            outside_bounding_nodes = zeros(2,2);
            if (this_xi == 1)
                inside_node = get_xright_node(this_node);
                for byi = 1:2; for bzi = 1:2;
                    outside_bounding_nodes(byi,bzi) = get_xleft_node(bounding_nodes(byi,bzi));
                end; end
            elseif (this_xi == xlength(this_idgrid));
                inside_node = get_xleft_node(this_node);
                for byi = 1:2; for bzi = 1:2;
                    outside_bounding_nodes(byi,bzi) = get_xright_node(bounding_nodes(byi,bzi));
                end; end
            end
            assert( all(all(outside_bounding_nodes ~= 0)) );
            assert( ( outside_bounding_nodes(1,1) ~= outside_bounding_nodes(1,2) ) || ( outside_bounding_nodes(1,1) ~= outside_bounding_nodes(2,1) ) );

            assert( all( nodedata.y(bounding_nodes(2,:)) >= nodedata.y(bounding_nodes(1,:)) ) );
            assert( all( nodedata.z(bounding_nodes(:,2)) >= nodedata.z(bounding_nodes(:,1)) ) );

            h1 = nodedata.y(bounding_nodes(2,1)) - nodedata.y(this_node);
            h2 = nodedata.y(this_node) - nodedata.y(bounding_nodes(1,1));

            h3 = nodedata.z(bounding_nodes(1,2)) - nodedata.z(this_node);
            h4 = nodedata.z(this_node) - nodedata.z(bounding_nodes(1,1));

            h5 = nodedata.x(outside_bounding_nodes(1,1)) - nodedata.x(this_node);
            assert(abs(h5-(nodedata.x(outside_bounding_nodes(2,1)) - nodedata.x(this_node)))<0.5*alignment_tol);
            assert(abs(h5-(nodedata.x(outside_bounding_nodes(1,2)) - nodedata.x(this_node)))<0.5*alignment_tol);
            assert(abs(h5-(nodedata.x(outside_bounding_nodes(2,2)) - nodedata.x(this_node)))<0.5*alignment_tol);

            h6 = nodedata.x(this_node) - nodedata.x(inside_node);

            % Check if aligned with one shared axis
            y_axis_aligned = false;
            z_axis_aligned = false;
            if (abs(h1) < alignment_tol) || (abs(h2) < alignment_tol)
                assert( (abs(h1) < alignment_tol) && (abs(h2) < alignment_tol))
                y_axis_aligned = true;
            end
            if (abs(h3) < alignment_tol) || (abs(h4) < alignment_tol)
                assert( (abs(h3) < alignment_tol) && (abs(h4) < alignment_tol))
                assert(~y_axis_aligned)
                z_axis_aligned = true;
            end

            base_wgt = -h6/(h5*(h5+h6));
            if ismember(this_type, [INT_X_FACE]) && (~y_axis_aligned) && (~z_axis_aligned)
                wgt_ydown_zback = base_wgt*h1*h3/( (h1+h2)*(h3+h4) );
                wgt_yup_zback   = base_wgt*h2*h3/( (h1+h2)*(h3+h4) );
                wgt_ydown_zfore = base_wgt*h1*h4/( (h1+h2)*(h3+h4) );
                wgt_yup_zfore   = base_wgt*h2*h4/( (h1+h2)*(h3+h4) );
                dx(this_node, bounding_nodes(1,1)) = wgt_ydown_zback;
                dx(this_node, bounding_nodes(2,1)) = wgt_yup_zback;
                dx(this_node, bounding_nodes(1,2)) = wgt_ydown_zfore;
                dx(this_node, bounding_nodes(2,2)) = wgt_yup_zfore;
                dx(this_node, outside_bounding_nodes(1,1)) = -dx(this_node, bounding_nodes(1,1));
                dx(this_node, outside_bounding_nodes(2,1)) = -dx(this_node, bounding_nodes(2,1));
                dx(this_node, outside_bounding_nodes(1,2)) = -dx(this_node, bounding_nodes(1,2));
                dx(this_node, outside_bounding_nodes(2,2)) = -dx(this_node, bounding_nodes(2,2));
            elseif ismember(this_type, [INT_XY_EDGE]) || (y_axis_aligned) % no variation in y
                assert(all(bounding_nodes(1,:)==bounding_nodes(2,:)));
                assert(all(outside_bounding_nodes(1,:)==outside_bounding_nodes(2,:)));
                wgt_zback = base_wgt*h3/( (h3+h4) );
                wgt_zfore = base_wgt*h4/( (h3+h4) );
                dx(this_node, bounding_nodes(1,1)) = wgt_zback;
                dx(this_node, bounding_nodes(1,2)) = wgt_zfore;
                dx(this_node, outside_bounding_nodes(1,1)) = -dx(this_node, bounding_nodes(1,1));
                dx(this_node, outside_bounding_nodes(1,2)) = -dx(this_node, bounding_nodes(1,2));
            elseif ismember(this_type, [INT_ZX_EDGE]) || (z_axis_aligned) % no variation in z
                assert(all(bounding_nodes(:,1)==bounding_nodes(:,2)));
                assert(all(outside_bounding_nodes(:,1)==outside_bounding_nodes(:,2)));
                wgt_ydown = base_wgt*h1/( (h1+h2) );
                wgt_yup = base_wgt*h2/( (h1+h2) );
                dx(this_node, bounding_nodes(1,1)) = wgt_ydown;
                dx(this_node, bounding_nodes(2,1)) = wgt_yup;
                dx(this_node, outside_bounding_nodes(1,1)) = -dx(this_node, bounding_nodes(1,1));
                dx(this_node, outside_bounding_nodes(2,1)) = -dx(this_node, bounding_nodes(2,1));
            end
            dx(this_node, this_node) = h5/(h6*(h5+h6));
            dx(this_node, inside_node) = -h5/(h6*(h5+h6));

            d2x(this_node, this_node) = dx(this_node, this_node)*-(2/h5);
            d2x(this_node, inside_node) = dx(this_node, inside_node)*-(2/h5);
            for byi = 1:2; for bzi = 1:2;
                d2x(this_node, bounding_nodes(byi,bzi)) = dx(this_node, bounding_nodes(byi,bzi))*(2/h6);
                d2x(this_node, outside_bounding_nodes(byi,bzi)) = dx(this_node, outside_bounding_nodes(byi,bzi))*(2/h6);
            end; end
        end

        if (show_waitbar)
            waitbar(progress/length(nontrivial_nodes_dx), f_wait, sprintf('Constructing dx & d2x operators...'));
        end
    end


    % Construct dy and d2y edge cases
    if (debug); disp('Constructing dy & d2y edge cases...'); end
    if (show_waitbar)
        waitbar(0, f_wait, sprintf('Constructing dy & d2y operators...'));
    end

    progress = 0;
    for this_node = nontrivial_nodes_dy
        progress = progress + 1;
        this_type         = nodedata.type(this_node);
        this_gridindex    = nodedata.grid(this_node);
        this_idgrid       = idgrids{this_gridindex};
        this_xi           = nodedata.xi(this_node);
        this_yi           = nodedata.yi(this_node);
        this_zi           = nodedata.zi(this_node);
        this_subgridindex = nodedata.subgrid(this_node);
        this_xi_sub       = nodedata.xi_sub(this_node);
        this_yi_sub       = nodedata.yi_sub(this_node);
        this_zi_sub       = nodedata.zi_sub(this_node);
        this_x            = nodedata.x(this_node);
        this_y            = nodedata.y(this_node);
        this_z            = nodedata.z(this_node);

        if (this_subgridindex == 0)
            this_subidgrid = [];
        else
            this_subidgrid = idgrids{this_subgridindex};
        end

        % clear this row of dy / d2y
        dy(this_node, :) = 0;
        d2y(this_node, :) = 0;

        if ismember(this_type, [NORMAL, SHARED, INT_Z_FACE, INT_X_FACE, INT_ZX_EDGE, EXT_Z_FACE, EXT_X_FACE, EXT_ZX_EDGE])
            % Easiest case, just do central difference

            ydown_node = get_ydown_node(this_node);
            yup_node = get_yup_node(this_node);
            assert((ydown_node ~= 0) && (yup_node ~= 0));
            ydown_y = nodedata.y(ydown_node);
            yup_y = nodedata.y(yup_node);

            h3 = (yup_y - this_y);
            h4 = (this_y - ydown_y);
            dy(this_node, ydown_node) = -(h3/h4)/(h3 + h4);
            dy(this_node, this_node) = 1/h4 - 1/h3;
            dy(this_node, yup_node) = (h4/h3)/(h3 + h4);

            d2y(this_node, ydown_node) = 2/(h4*(h3 + h4));
            d2y(this_node, this_node) = -2/(h3*h4);
            d2y(this_node, yup_node) = 2/(h3*(h3 + h4));
        elseif ismember(this_type, [EXT_Y_FACE, EXT_YZ_EDGE, EXT_XY_EDGE, EXT_CORNER])
            % We're at the true edge of the whole mesh

            % use first or second order difference, depending on edge_order
            if (this_yi == 1)
                next_node = get_yup_node(this_node);
                far_node = get_yup_node(next_node);
            elseif (this_yi == ylength(this_idgrid));
                next_node = get_ydown_node(this_node);
                far_node = get_ydown_node(next_node);
            else
                assert(false)
            end
            assert((next_node ~= 0) && (far_node ~= 0));
            h3 = nodedata.y(next_node) - nodedata.y(this_node);
            h5 = nodedata.y(far_node) - nodedata.y(next_node);

            if (edge_order == 1)
                dy(this_node, this_node) = -1/h3;
                dy(this_node, next_node) = 1/h3;

                d2y(this_node, this_node) = 2/(h3*(h3+h5));
                d2y(this_node, next_node) = -2/(h3*h5);
                d2y(this_node, far_node) = 2/(h5*(h3+h5));
            else
                dy(this_node, this_node) = -(2*h3+h5)/(h3*(h3+h5));
                dy(this_node, next_node) = (1/h3 + 1/h5);
                dy(this_node, far_node) = -h3/(h5*(h3+h5));

                d2y(this_node, this_node) = 2/(h3*(h3+h5));
                d2y(this_node, next_node) = -2/(h3*h5);
                d2y(this_node, far_node) = 2/(h5*(h3+h5));
            end
        elseif ismember(this_type, [INT_YZ_EDGE, INT_XY_EDGE, INT_Y_FACE])
            % We're at a mesh boundary within the system; we use 10 nearby points to construct
            % the correct second order difference (time for the ghost nodes).

            assert(this_subgridindex == 0);
            assert(this_gridindex ~= 0);

            % 1) get nearest neighbors
            bounding_nodes = get_bounding_shared_nodes(this_node, nodedata, idgrids, parentlist, xss, yss, zss, alignment_tol);
            assert(all(all(bounding_nodes(:,1,:) == bounding_nodes(:,2,:))));
            % get rid of y data in bounding_nodes
            bounding_nodes = permute(bounding_nodes, [3 1 2]);
            bounding_nodes = bounding_nodes(:,:,1);

            % 2) get side neighbors
            outside_bounding_nodes = zeros(2,2);
            if (this_yi == 1)
                inside_node = get_yup_node(this_node);
                for bzi = 1:2; for bxi = 1:2;
                    outside_bounding_nodes(bzi,bxi) = get_ydown_node(bounding_nodes(bzi,bxi));
                end; end
            elseif (this_yi == ylength(this_idgrid));
                inside_node = get_ydown_node(this_node);
                for bzi = 1:2; for bxi = 1:2;
                    outside_bounding_nodes(bzi,bxi) = get_yup_node(bounding_nodes(bzi,bxi));
                end; end
            end
            assert( all(all(outside_bounding_nodes ~= 0)) );
            assert( ( outside_bounding_nodes(1,1) ~= outside_bounding_nodes(1,2) ) || ( outside_bounding_nodes(1,1) ~= outside_bounding_nodes(2,1) ) );

            assert( all( nodedata.z(bounding_nodes(2,:)) >= nodedata.z(bounding_nodes(1,:)) ) );
            assert( all( nodedata.x(bounding_nodes(:,2)) >= nodedata.x(bounding_nodes(:,1)) ) );

            h1 = nodedata.z(bounding_nodes(2,1)) - nodedata.z(this_node);
            h2 = nodedata.z(this_node) - nodedata.z(bounding_nodes(1,1));

            h3 = nodedata.x(bounding_nodes(1,2)) - nodedata.x(this_node);
            h4 = nodedata.x(this_node) - nodedata.x(bounding_nodes(1,1));

            h5 = nodedata.y(outside_bounding_nodes(1,1)) - nodedata.y(this_node);
            assert(abs(h5-(nodedata.y(outside_bounding_nodes(2,1)) - nodedata.y(this_node)))<0.5*alignment_tol);
            assert(abs(h5-(nodedata.y(outside_bounding_nodes(1,2)) - nodedata.y(this_node)))<0.5*alignment_tol);
            assert(abs(h5-(nodedata.y(outside_bounding_nodes(2,2)) - nodedata.y(this_node)))<0.5*alignment_tol);

            h6 = nodedata.y(this_node) - nodedata.y(inside_node);

            % Check if aligned with one shared axis
            z_axis_aligned = false;
            x_axis_aligned = false;
            if (abs(h1) < alignment_tol) || (abs(h2) < alignment_tol)
                assert( (abs(h1) < alignment_tol) && (abs(h2) < alignment_tol))
                z_axis_aligned = true;
            end
            if (abs(h3) < alignment_tol) || (abs(h4) < alignment_tol)
                assert( (abs(h3) < alignment_tol) && (abs(h4) < alignment_tol))
                assert(~y_axis_aligned)
                x_axis_aligned = true;
            end

            base_wgt = -h6/(h5*(h5+h6));
            if ismember(this_type, [INT_Y_FACE]) && (~z_axis_aligned) && (~x_axis_aligned)
                wgt_zback_xleft = base_wgt*h1*h3/( (h1+h2)*(h3+h4) );
                wgt_zfore_xleft   = base_wgt*h2*h3/( (h1+h2)*(h3+h4) );
                wgt_zback_xright = base_wgt*h1*h4/( (h1+h2)*(h3+h4) );
                wgt_zfore_xright   = base_wgt*h2*h4/( (h1+h2)*(h3+h4) );
                dy(this_node, bounding_nodes(1,1)) = wgt_zback_xleft;
                dy(this_node, bounding_nodes(2,1)) = wgt_zfore_xleft;
                dy(this_node, bounding_nodes(1,2)) = wgt_zback_xright;
                dy(this_node, bounding_nodes(2,2)) = wgt_zfore_xright;
                dy(this_node, outside_bounding_nodes(1,1)) = -dy(this_node, bounding_nodes(1,1));
                dy(this_node, outside_bounding_nodes(2,1)) = -dy(this_node, bounding_nodes(2,1));
                dy(this_node, outside_bounding_nodes(1,2)) = -dy(this_node, bounding_nodes(1,2));
                dy(this_node, outside_bounding_nodes(2,2)) = -dy(this_node, bounding_nodes(2,2));
            elseif ismember(this_type, [INT_YZ_EDGE]) || (z_axis_aligned) % no variation in z
                assert(all(bounding_nodes(1,:)==bounding_nodes(2,:)));
                assert(all(outside_bounding_nodes(1,:)==outside_bounding_nodes(2,:)));
                wgt_xleft = base_wgt*h3/( (h3+h4) );
                wgt_xright = base_wgt*h4/( (h3+h4) );
                dy(this_node, bounding_nodes(1,1)) = wgt_xleft;
                dy(this_node, bounding_nodes(1,2)) = wgt_xright;
                dy(this_node, outside_bounding_nodes(1,1)) = -dy(this_node, bounding_nodes(1,1));
                dy(this_node, outside_bounding_nodes(1,2)) = -dy(this_node, bounding_nodes(1,2));
            elseif ismember(this_type, [INT_XY_EDGE]) || (x_axis_aligned) % no variation in x
                assert(all(bounding_nodes(:,1)==bounding_nodes(:,2)));
                assert(all(outside_bounding_nodes(:,1)==outside_bounding_nodes(:,2)));
                wgt_zback = base_wgt*h1/( (h1+h2) );
                wgt_zfore = base_wgt*h2/( (h1+h2) );
                dy(this_node, bounding_nodes(1,1)) = wgt_zback;
                dy(this_node, bounding_nodes(2,1)) = wgt_zfore;
                dy(this_node, outside_bounding_nodes(1,1)) = -dy(this_node, bounding_nodes(1,1));
                dy(this_node, outside_bounding_nodes(2,1)) = -dy(this_node, bounding_nodes(2,1));
            end
            dy(this_node, this_node) = h5/(h6*(h5+h6));
            dy(this_node, inside_node) = -h5/(h6*(h5+h6));

            d2y(this_node, this_node) = dy(this_node, this_node)*-(2/h5);
            d2y(this_node, inside_node) = dy(this_node, inside_node)*-(2/h5);
            for bzi = 1:2; for bxi = 1:2;
                d2y(this_node, bounding_nodes(bzi,bxi)) = dy(this_node, bounding_nodes(bzi,bxi))*(2/h6);
                d2y(this_node, outside_bounding_nodes(bzi,bxi)) = dy(this_node, outside_bounding_nodes(bzi,bxi))*(2/h6);
            end; end
        end

        if (show_waitbar)
            waitbar(progress/length(nontrivial_nodes_dy), f_wait, sprintf('Constructing dy & d2y operators...'));
        end
    end


    % Construct dz and d2z edge cases
    if (debug); disp('Constructing dz & d2z edge cases...'); end
    if (show_waitbar)
        waitbar(0, f_wait, sprintf('Constructing dz & d2z operators...'));
    end

    progress = 0;
    for this_node = nontrivial_nodes_dz
        progress = progress + 1;
        this_type         = nodedata.type(this_node);
        this_gridindex    = nodedata.grid(this_node);
        this_idgrid       = idgrids{this_gridindex};
        this_xi           = nodedata.xi(this_node);
        this_yi           = nodedata.yi(this_node);
        this_zi           = nodedata.zi(this_node);
        this_subgridindex = nodedata.subgrid(this_node);
        this_xi_sub       = nodedata.xi_sub(this_node);
        this_yi_sub       = nodedata.yi_sub(this_node);
        this_zi_sub       = nodedata.zi_sub(this_node);
        this_x            = nodedata.x(this_node);
        this_y            = nodedata.y(this_node);
        this_z            = nodedata.z(this_node);

        if (this_subgridindex == 0)
            this_subidgrid = [];
        else
            this_subidgrid = idgrids{this_subgridindex};
        end

        % clear this row of dz / d2z
        dz(this_node, :) = 0;
        d2z(this_node, :) = 0;

        if ismember(this_type, [NORMAL, SHARED, INT_X_FACE, INT_Y_FACE, INT_XY_EDGE, EXT_X_FACE, EXT_Y_FACE, EXT_XY_EDGE])
            % Easiest case, just do central difference

            zback_node = get_zback_node(this_node);
            zfore_node = get_zfore_node(this_node);
            assert((zback_node ~= 0) && (zfore_node ~= 0));
            zback_z = nodedata.z(zback_node);
            zfore_z = nodedata.z(zfore_node);

            h3 = (zfore_z - this_z);
            h4 = (this_z - zback_z);
            dz(this_node, zback_node) = -(h3/h4)/(h3 + h4);
            dz(this_node, this_node) = 1/h4 - 1/h3;
            dz(this_node, zfore_node) = (h4/h3)/(h3 + h4);

            d2z(this_node, zback_node) = 2/(h4*(h3 + h4));
            d2z(this_node, this_node) = -2/(h3*h4);
            d2z(this_node, zfore_node) = 2/(h3*(h3 + h4));
        elseif ismember(this_type, [EXT_Z_FACE, EXT_ZX_EDGE, EXT_YZ_EDGE, EXT_CORNER])
            % We're at the true edge of the whole mesh

            % use first or second order difference, depending on edge_order
            if (this_zi == 1)
                next_node = get_zfore_node(this_node);
                far_node = get_zfore_node(next_node);
            elseif (this_zi == zlength(this_idgrid));
                next_node = get_zback_node(this_node);
                far_node = get_zback_node(next_node);
            else
                assert(false)
            end
            assert((next_node ~= 0) && (far_node ~= 0));
            h3 = nodedata.z(next_node) - nodedata.z(this_node);
            h5 = nodedata.z(far_node) - nodedata.z(next_node);

            if (edge_order == 1)
                dz(this_node, this_node) = -1/h3;
                dz(this_node, next_node) = 1/h3;

                d2z(this_node, this_node) = 2/(h3*(h3+h5));
                d2z(this_node, next_node) = -2/(h3*h5);
                d2z(this_node, far_node) = 2/(h5*(h3+h5));
            else
                dz(this_node, this_node) = -(2*h3+h5)/(h3*(h3+h5));
                dz(this_node, next_node) = (1/h3 + 1/h5);
                dz(this_node, far_node) = -h3/(h5*(h3+h5));

                d2z(this_node, this_node) = 2/(h3*(h3+h5));
                d2z(this_node, next_node) = -2/(h3*h5);
                d2z(this_node, far_node) = 2/(h5*(h3+h5));
            end
        elseif ismember(this_type, [INT_ZX_EDGE, INT_YZ_EDGE, INT_Z_FACE])
            % We're at a mesh boundary within the system; we use 10 nearby points to construct
            % the correct second order difference (time for the ghost nodes).

            assert(this_subgridindex == 0);
            assert(this_gridindex ~= 0);

            % 1) get nearest neighbors
            bounding_nodes = get_bounding_shared_nodes(this_node, nodedata, idgrids, parentlist, xss, yss, zss, alignment_tol);
            assert(all(all(bounding_nodes(:,:,1) == bounding_nodes(:,:,2))));
            % get rid of z data in bounding_nodes
            %bounding_nodes = permute(bounding_nodes, [1 2 3]);
            bounding_nodes = bounding_nodes(:,:,1);

            % 2) get side neighbors
            outside_bounding_nodes = zeros(2,2);
            if (this_zi == 1)
                inside_node = get_zfore_node(this_node);
                for bxi = 1:2; for byi = 1:2;
                    outside_bounding_nodes(bxi,byi) = get_zback_node(bounding_nodes(bxi,byi));
                end; end
            elseif (this_zi == zlength(this_idgrid));
                inside_node = get_zback_node(this_node);
                for bxi = 1:2; for byi = 1:2;
                    outside_bounding_nodes(bxi,byi) = get_zfore_node(bounding_nodes(bxi,byi));
                end; end
            end
            assert( all(all(outside_bounding_nodes ~= 0)) );
            assert( ( outside_bounding_nodes(1,1) ~= outside_bounding_nodes(1,2) ) || ( outside_bounding_nodes(1,1) ~= outside_bounding_nodes(2,1) ) );

            assert( all( nodedata.x(bounding_nodes(2,:)) >= nodedata.x(bounding_nodes(1,:)) ) );
            assert( all( nodedata.y(bounding_nodes(:,2)) >= nodedata.y(bounding_nodes(:,1)) ) );

            h1 = nodedata.x(bounding_nodes(2,1)) - nodedata.x(this_node);
            h2 = nodedata.x(this_node) - nodedata.x(bounding_nodes(1,1));

            h3 = nodedata.y(bounding_nodes(1,2)) - nodedata.y(this_node);
            h4 = nodedata.y(this_node) - nodedata.y(bounding_nodes(1,1));

            h5 = nodedata.z(outside_bounding_nodes(1,1)) - nodedata.z(this_node);
            assert(abs(h5 - (nodedata.z(outside_bounding_nodes(2,1)) - nodedata.z(this_node))) < alignment_tol*0.5);
            assert(abs(h5 - (nodedata.z(outside_bounding_nodes(1,2)) - nodedata.z(this_node))) < alignment_tol*0.5);
            assert(abs(h5 - (nodedata.z(outside_bounding_nodes(2,2)) - nodedata.z(this_node))) < alignment_tol*0.5);

            h6 = nodedata.z(this_node) - nodedata.z(inside_node);

            % Check if aligned with one shared axis
            x_axis_aligned = false;
            y_axis_aligned = false;
            if (abs(h1) < alignment_tol) || (abs(h2) < alignment_tol)
                assert( (abs(h1) < alignment_tol) && (abs(h2) < alignment_tol))
                x_axis_aligned = true;
            end
            if (abs(h3) < alignment_tol) || (abs(h4) < alignment_tol)
                assert( (abs(h3) < alignment_tol) && (abs(h4) < alignment_tol))
                assert(~y_axis_aligned)
                y_axis_aligned = true;
            end

            base_wgt = -h6/(h5*(h5+h6));
            if ismember(this_type, [INT_Z_FACE]) && (~x_axis_aligned) && (~y_axis_aligned)
                wgt_xleft_ydown = base_wgt*h1*h3/( (h1+h2)*(h3+h4) );
                wgt_xright_ydown   = base_wgt*h2*h3/( (h1+h2)*(h3+h4) );
                wgt_xleft_yup = base_wgt*h1*h4/( (h1+h2)*(h3+h4) );
                wgt_xright_yup   = base_wgt*h2*h4/( (h1+h2)*(h3+h4) );
                dz(this_node, bounding_nodes(1,1)) = wgt_xleft_ydown;
                dz(this_node, bounding_nodes(2,1)) = wgt_xright_ydown;
                dz(this_node, bounding_nodes(1,2)) = wgt_xleft_yup;
                dz(this_node, bounding_nodes(2,2)) = wgt_xright_yup;
                dz(this_node, outside_bounding_nodes(1,1)) = -dz(this_node, bounding_nodes(1,1));
                dz(this_node, outside_bounding_nodes(2,1)) = -dz(this_node, bounding_nodes(2,1));
                dz(this_node, outside_bounding_nodes(1,2)) = -dz(this_node, bounding_nodes(1,2));
                dz(this_node, outside_bounding_nodes(2,2)) = -dz(this_node, bounding_nodes(2,2));
            elseif ismember(this_type, [INT_ZX_EDGE]) || (x_axis_aligned) % no variation in x
                assert(all(bounding_nodes(1,:)==bounding_nodes(2,:)));
                assert(all(outside_bounding_nodes(1,:)==outside_bounding_nodes(2,:)));
                wgt_ydown = base_wgt*h3/( (h3+h4) );
                wgt_yup = base_wgt*h4/( (h3+h4) );
                dz(this_node, bounding_nodes(1,1)) = wgt_ydown;
                dz(this_node, bounding_nodes(1,2)) = wgt_yup;
                dz(this_node, outside_bounding_nodes(1,1)) = -dz(this_node, bounding_nodes(1,1));
                dz(this_node, outside_bounding_nodes(1,2)) = -dz(this_node, bounding_nodes(1,2));
            elseif ismember(this_type, [INT_YZ_EDGE]) || (y_axis_aligned) % no variation in y
                assert(all(bounding_nodes(:,1)==bounding_nodes(:,2)));
                assert(all(outside_bounding_nodes(:,1)==outside_bounding_nodes(:,2)));
                wgt_xleft = base_wgt*h1/( (h1+h2) );
                wgt_xright = base_wgt*h2/( (h1+h2) );
                dz(this_node, bounding_nodes(1,1)) = wgt_xleft;
                dz(this_node, bounding_nodes(2,1)) = wgt_xright;
                dz(this_node, outside_bounding_nodes(1,1)) = -dz(this_node, bounding_nodes(1,1));
                dz(this_node, outside_bounding_nodes(2,1)) = -dz(this_node, bounding_nodes(2,1));
            end
            dz(this_node, this_node) = h5/(h6*(h5+h6));
            dz(this_node, inside_node) = -h5/(h6*(h5+h6));

            d2z(this_node, this_node) = dz(this_node, this_node)*-(2/h5);
            d2z(this_node, inside_node) = dz(this_node, inside_node)*-(2/h5);
            for bxi = 1:2; for byi = 1:2;
                d2z(this_node, bounding_nodes(bxi,byi)) = dz(this_node, bounding_nodes(bxi,byi))*(2/h6);
                d2z(this_node, outside_bounding_nodes(bxi,byi)) = dz(this_node, outside_bounding_nodes(bxi,byi))*(2/h6);
            end; end
        end

        if (show_waitbar)
            waitbar(progress/length(nontrivial_nodes_dz), f_wait, sprintf('Constructing dz & d2z operators...'));
        end
    end

    %%%%%%%%%% CONSTRUCT INTEGRATION WEIGHT VECTOR %%%%%%%%%%

    % I.e., a vector giving the area of a small rectangle around each node.
    if (debug); disp('Calculating integration weights...'); end

    % TODO: Construct with a proven nth order accuracy
    int_weight = zeros(nodeN, 1);

    nontrivial_nodes_intwgt = all_nodes((nodedata.type ~= NORMAL));

    % construct integration weight for each mesh...
    for j = 1:meshN
        this_xs = xss{j};
        this_ys = yss{j};
        this_zs = zss{j};
        this_idgrid = idgrids{j};

        this_x_lengths = 0.5*([0 diff(this_xs)] + [diff(this_xs) 0]);
        this_y_lengths = 0.5*([0 diff(this_ys)] + [diff(this_ys) 0]);
        this_z_lengths = 0.5*([0 diff(this_zs)] + [diff(this_zs) 0]);

        for xi = 1:length(this_xs); for yi = 1:length(this_ys); for zi = 1:length(this_zs)
                    this_node = this_idgrid(xi, yi, zi);
                    int_weight(this_node) = this_x_lengths(xi)*this_y_lengths(yi)*this_z_lengths(zi);
        end; end; end

    end


    progress = 0;
    if (show_waitbar)
        waitbar(progress/length(nontrivial_nodes_intwgt), f_wait, sprintf('Calculating integration weights...'));
    end

    for this_node = nontrivial_nodes_intwgt
        this_type         = nodedata.type(this_node);
        this_gridindex    = nodedata.grid(this_node);
        this_idgrid       = idgrids{this_gridindex};
        this_xi           = nodedata.xi(this_node);
        this_yi           = nodedata.yi(this_node);
        this_zi           = nodedata.zi(this_node);
        this_subgridindex = nodedata.subgrid(this_node);
        this_xi_sub       = nodedata.xi_sub(this_node);
        this_yi_sub       = nodedata.yi_sub(this_node);
        this_zi_sub       = nodedata.zi_sub(this_node);
        this_x            = nodedata.x(this_node);
        this_y            = nodedata.y(this_node);
        this_z            = nodedata.z(this_node);

        progress = progress + 1;

        xleft_node = get_xleft_node(this_node);
        if (xleft_node == 0)
            if abs(this_x - xmin_all) < alignment_tol
                xleft_node = 0; % we're on the total mesh edge
            else
                % get bounding nodes and ensure we're on an x-facing face/edge
                bounding_nodes = get_bounding_shared_nodes(this_node, nodedata, idgrids, parentlist, xss, yss, zss, alignment_tol);
                assert(all(all(bounding_nodes(1,:,:) == bounding_nodes(2,:,:))));
                assert(~all(all(bounding_nodes(1,:,:) == 0)));
                shared_node = 0;
                for i1=1:2;
                    for i2=1:2;
                        shared_node = bounding_nodes(1,i1,i2);
                        if (shared_node ~= 0)
                            break
                        end
                    end
                    if (shared_node ~= 0)
                        break
                    end
                end

                xleft_node = get_xleft_node(shared_node);
            end
        end

        xright_node = get_xright_node(this_node);
        if (xright_node == 0)
            if abs(this_x - xmax_all) < alignment_tol
                xright_node = 0; % we're on the total mesh edge
            else
                % get bounding nodes and ensure we're on an x-facing face/edge
                bounding_nodes = get_bounding_shared_nodes(this_node, nodedata, idgrids, parentlist, xss, yss, zss, alignment_tol);
                assert(all(all(bounding_nodes(1,:,:) == bounding_nodes(2,:,:))));
                assert(~all(all(bounding_nodes(1,:,:) == 0)));
                shared_node = 0;
                for i1=1:2;
                    for i2=1:2;
                        shared_node = bounding_nodes(1,i1,i2);
                        if (shared_node ~= 0)
                            break
                        end
                    end
                    if (shared_node ~= 0)
                        break
                    end
                end

                xright_node = get_xright_node(shared_node);
            end
        end

        ydown_node = get_ydown_node(this_node);
        if (ydown_node == 0)
            if abs(this_y - ymin_all) < alignment_tol
                ydown_node = 0; % we're on the total mesh edge
            else
                % get bounding nodes and ensure we're on an y-facing face/edge
                bounding_nodes = get_bounding_shared_nodes(this_node, nodedata, idgrids, parentlist, xss, yss, zss, alignment_tol);
                assert(all(all(bounding_nodes(:,1,:) == bounding_nodes(:,2,:))));
                assert(~all(all(bounding_nodes(:,1,:) == 0)));
                shared_node = 0;
                for i1=1:2;
                    for i2=1:2;
                        shared_node = bounding_nodes(i1,1,i2);
                        if (shared_node ~= 0)
                            break
                        end
                    end
                    if (shared_node ~= 0)
                        break
                    end
                end

                ydown_node = get_ydown_node(shared_node);
            end
        end

        yup_node = get_yup_node(this_node);
        if (yup_node == 0)
            if abs(this_y - ymax_all) < alignment_tol
                yup_node = 0; % we're on the total mesh edge
            else
                % get bounding nodes and ensure we're on an y-facing face/edge
                bounding_nodes = get_bounding_shared_nodes(this_node, nodedata, idgrids, parentlist, xss, yss, zss, alignment_tol);
                assert(all(all(bounding_nodes(:,1,:) == bounding_nodes(:,2,:))));
                assert(~all(all(bounding_nodes(:,1,:) == 0)));
                shared_node = 0;
                for i1=1:2;
                    for i2=1:2;
                        shared_node = bounding_nodes(i1,1,i2);
                        if (shared_node ~= 0)
                            break
                        end
                    end
                    if (shared_node ~= 0)
                        break
                    end
                end

                yup_node = get_yup_node(shared_node);
            end
        end

        zback_node = get_zback_node(this_node);
        if (zback_node == 0)
            if abs(this_z - zmin_all) < alignment_tol
                zback_node = 0; % we're on the total mesh edge
            else
                % get bounding nodes and ensure we're on an z-facing face/edge
                assert(all(all(bounding_nodes(:,:,1) == bounding_nodes(:,:,2))));
                assert(~all(all(bounding_nodes(:,:,1) == 0)));
                shared_node = 0;
                for i1=1:2;
                    for i2=1:2;
                        shared_node = bounding_nodes(i1,i2,1);
                        if (shared_node ~= 0)
                            break
                        end
                    end
                    if (shared_node ~= 0)
                        break
                    end
                end

                zback_node = get_zback_node(shared_node);
            end
        end

        zfore_node = get_zfore_node(this_node);
        if (zfore_node == 0)
            if abs(this_z - zmax_all) < alignment_tol
                zfore_node = 0; % we're on the total mesh edge
            else
                % get bounding nodes and ensure we're on an z-facing face/edge
                assert(all(all(bounding_nodes(:,:,1) == bounding_nodes(:,:,2))));
                assert(~all(all(bounding_nodes(:,:,1) == 0)));
                shared_node = 0;
                for i1=1:2;
                    for i2=1:2;
                        shared_node = bounding_nodes(i1,i2,1);
                        if (shared_node ~= 0)
                            break
                        end
                    end
                    if (shared_node ~= 0)
                        break
                    end
                end

                zfore_node = get_zfore_node(shared_node);
            end
        end


        dist_xleft = 0;
        dist_xright = 0;
        if (xleft_node ~= 0)
            dist_xleft = (nodedata.x(this_node) - nodedata.x(xleft_node));
        end
        if (xright_node ~= 0)
            dist_xright = (nodedata.x(xright_node) - nodedata.x(this_node));
        end
        this_x_length = 0.5*(dist_xright + dist_xleft);
        assert(this_x_length > 0);

        dist_ydown = 0;
        dist_yup = 0;
        if (ydown_node ~= 0)
            dist_ydown = (nodedata.y(this_node) - nodedata.y(ydown_node));
        end
        if (yup_node ~= 0)
            dist_yup = (nodedata.y(yup_node) - nodedata.y(this_node));
        end
        this_y_length = 0.5*(dist_yup + dist_ydown);
        assert(this_y_length > 0);

        dist_zback = 0;
        dist_zfore = 0;
        if (zback_node ~= 0)
            dist_zback = (nodedata.z(this_node) - nodedata.z(zback_node));
        end
        if (zfore_node ~= 0)
            dist_zfore = (nodedata.z(zfore_node) - nodedata.z(this_node));
        end
        this_z_length = 0.5*(dist_zfore + dist_zback);
        assert(this_z_length > 0);


        int_weight(this_node) = this_x_length*this_y_length*this_z_length;

        if (show_waitbar)
            waitbar(progress/length(nontrivial_nodes_intwgt), f_wait, sprintf('Calculating integration weights...'));
        end
    end

    if (show_waitbar)
        delete(f_wait);
    end

end

%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%

function next_node = get_neighbor_node(this_node, direction, nodedata, idgrids)
    % Direction: 'xleft', 'xright', 'ydown', 'yup', 'zback', 'zfore'

    % Define node types
    NORMAL      = 0;
    SHARED      = 1;
    INT_X_FACE  = 2;
    INT_Y_FACE  = 3;
    INT_Z_FACE  = 4;
    INT_XY_EDGE = 5;
    INT_YZ_EDGE = 6;
    INT_ZX_EDGE = 7;
    EXT_X_FACE  = 8;
    EXT_Y_FACE  = 9;
    EXT_Z_FACE  = 10;
    EXT_XY_EDGE = 11;
    EXT_YZ_EDGE = 12;
    EXT_ZX_EDGE = 13;
    EXT_CORNER  = 14;
    xlength = @(idgrid) size(idgrid, 1);
    ylength = @(idgrid) size(idgrid, 2);
    zlength = @(idgrid) size(idgrid, 3);

    this_type         = nodedata.type(this_node);
    this_gridindex    = nodedata.grid(this_node);
    this_idgrid       = idgrids{this_gridindex};
    this_xi           = nodedata.xi(this_node);
    this_yi           = nodedata.yi(this_node);
    this_zi           = nodedata.zi(this_node);
    this_subgridindex = nodedata.subgrid(this_node);
    this_xi_sub       = nodedata.xi_sub(this_node);
    this_yi_sub       = nodedata.yi_sub(this_node);
    this_zi_sub       = nodedata.zi_sub(this_node);
    this_x            = nodedata.x(this_node);
    this_y            = nodedata.y(this_node);
    this_z            = nodedata.z(this_node);

    if strcmpi(direction,'xleft')
        not_on_edge = (this_xi > 1);
        sub_not_on_edge = (this_xi_sub > 1);
        xistep = -1;
        yistep = 0;
        zistep = 0;
    elseif strcmpi(direction,'xright')
        not_on_edge = (this_xi < xlength(this_idgrid));
        if (this_subgridindex ~= 0)
            sub_not_on_edge = (this_xi_sub < xlength(idgrids{this_subgridindex}));
        end
        xistep = 1;
        yistep = 0;
        zistep = 0;
    elseif strcmpi(direction,'ydown')
        not_on_edge = (this_yi > 1);
        sub_not_on_edge = (this_yi_sub > 1);
        xistep = 0;
        yistep = -1;
        zistep = 0;
    elseif strcmpi(direction,'yup')
        not_on_edge = (this_yi < ylength(this_idgrid));
        if (this_subgridindex ~= 0)
            sub_not_on_edge = (this_yi_sub < ylength(idgrids{this_subgridindex}));
        end
        xistep = 0;
        yistep = 1;
        zistep = 0;
    elseif strcmpi(direction,'zback')
        not_on_edge = (this_zi > 1);
        sub_not_on_edge = (this_zi_sub > 1);
        xistep = 0;
        yistep = 0;
        zistep = -1;
    elseif strcmpi(direction,'zfore')
        not_on_edge = (this_zi < zlength(this_idgrid));
        if (this_subgridindex ~= 0)
            sub_not_on_edge = (this_zi_sub < zlength(idgrids{this_subgridindex}));
        end
        xistep = 0;
        yistep = 0;
        zistep = 1;
    else
        assert('direction must be ''xleft'', ''xright'', ''ydown'', ''yup'', ''zback'', or ''zfore''.');
    end

    if this_subgridindex == 0 % no subgrid
        if not_on_edge
            next_node = this_idgrid(this_xi + xistep, this_yi + yistep, this_zi + zistep);
        else
            next_node = 0; % true edge
        end
    else % yes subgrid
        if sub_not_on_edge
            next_node = idgrids{this_subgridindex}(this_xi_sub + xistep, this_yi_sub + yistep, this_zi_sub + zistep);
        else
            % we're at an interior boundary; we can look at the parent mesh
            if not_on_edge
                next_node = this_idgrid(this_xi + xistep, this_yi + yistep, this_zi + zistep);
            else
                next_node = 0; % edge of parent mesh too
            end
        end
    end

end

function bounding_nodes = get_bounding_shared_nodes(this_node, nodedata, idgrids, parentlist, xss, yss, zss, alignment_tol)
    % Returned in the following order:
    %   [ (-X, -Y, -Z)  (-X, Y, -Z) ; (X, -Y, -Z)  (X, Y, -Z) ]
    %   [ (-X, -Y, Z)   (-X, Y, Z)  ; (X, -Y, Z)   (X, Y, Z)  ]

    % Define node types
    NORMAL      = 0;
    SHARED      = 1;
    INT_X_FACE  = 2;
    INT_Y_FACE  = 3;
    INT_Z_FACE  = 4;
    INT_XY_EDGE = 5;
    INT_YZ_EDGE = 6;
    INT_ZX_EDGE = 7;
    EXT_X_FACE  = 8;
    EXT_Y_FACE  = 9;
    EXT_Z_FACE  = 10;
    EXT_XY_EDGE = 11;
    EXT_YZ_EDGE = 12;
    EXT_ZX_EDGE = 13;
    EXT_CORNER  = 14;
    xlength = @(idgrid) size(idgrid, 1);
    ylength = @(idgrid) size(idgrid, 2);
    zlength = @(idgrid) size(idgrid, 3);

    this_type         = nodedata.type(this_node);
    this_gridindex    = nodedata.grid(this_node);
    this_idgrid       = idgrids{this_gridindex};
    this_xi           = nodedata.xi(this_node);
    this_yi           = nodedata.yi(this_node);
    this_zi           = nodedata.zi(this_node);
    this_subgridindex = nodedata.subgrid(this_node);
    this_x            = nodedata.x(this_node);
    this_y            = nodedata.y(this_node);
    this_z            = nodedata.z(this_node);

    bounding_nodes = zeros(2,2,2);

    % Fail if we start on a shared node
    if(this_subgridindex ~= 0);
        assert(false);
        %bounding_nodes(:,:,:) = this_node;
        return
    end

    root_gridindex = find(parentlist==0);
    if (this_gridindex == root_gridindex)% && (this_subgridindex == 0)
        % Fail if this is the root mesh and not a shared node
        assert(false);
        return
    end

    parent_gridindex = parentlist(this_gridindex);
    parent_idgrid = idgrids{parent_gridindex};
    parent_xs = xss{parent_gridindex};
    parent_ys = yss{parent_gridindex};
    parent_zs = zss{parent_gridindex};

    [parent_xi0, parent_x0] = nearest_index(parent_xs, this_x);
    if abs(parent_x0 - this_x) < alignment_tol
        parent_xi1 = parent_xi0;
        parent_xi2 = parent_xi0;
    elseif parent_x0 < this_x
        parent_xi1 = parent_xi0;
        parent_xi2 = parent_xi0 + 1;
    else % parent_x0 > this_x
        parent_xi2 = parent_xi0;
        parent_xi1 = parent_xi0 - 1;
    end

    [parent_yi0, parent_y0] = nearest_index(parent_ys, this_y);
    if abs(parent_y0 - this_y) < alignment_tol
        parent_yi1 = parent_yi0;
        parent_yi2 = parent_yi0;
    elseif parent_y0 < this_y
        parent_yi1 = parent_yi0;
        parent_yi2 = parent_yi0 + 1;
    else % parent_y0 > this_y
        parent_yi2 = parent_yi0;
        parent_yi1 = parent_yi0 - 1;
    end

    [parent_zi0, parent_z0] = nearest_index(parent_zs, this_z);
    if abs(parent_z0 - this_z) < alignment_tol
        parent_zi1 = parent_zi0;
        parent_zi2 = parent_zi0;
    elseif parent_z0 < this_z
        parent_zi1 = parent_zi0;
        parent_zi2 = parent_zi0 + 1;
    else % parent_z0 > this_z
        parent_zi2 = parent_zi0;
        parent_zi1 = parent_zi0 - 1;
    end


    bounding_nodes = parent_idgrid([parent_xi1 parent_xi2], [parent_yi1 parent_yi2], [parent_zi1 parent_zi2]);

end
