% MULTIMESH_2D Construct 2D multi-resolution mesh and derivative operators
%
% [idgrids, nodedata, dx, d2x, dy, d2y, int_weight] = MULTIMESH_2D(MESHES)
% [...] = MULTIMESH_2D(MESHES, ..., 'alignmentTol', ALIGNMENTTOL)
% [...] = MULTIMESH_2D(MESHES, ..., 'edgeOrder', EDGEORDER)
% [...] = MULTIMESH_2D(MESHES, ..., 'x0EdgeOrder', X0EDGEORDER)
% [...] = MULTIMESH_2D(MESHES, ..., 'x0Cylindrical', X0CYLINDRICAL)
% [...] = MULTIMESH_2D(MESHES, ..., 'cylindricalWeight', CYLINDRICALWEIGHT)
%
%
% Input:
%
%  MESHES must be a cell-list of cell-tuples; each cell-tuple containing a list
%  of x values and a list of y values to specify a mesh. For example, a valid
%  input would be
%  >> MESHES = { {[0 1 2], [0 1 2 3]}, {[1 1.2 1.3 2], [0 0.1 0.5 1 1.5 2]} }
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
%  Two special options are available for the minimal x edge, for use with
%  cylindrical coordinates: X0EDGEORDER and X0CYLINDRICAL. These may not be
%  used together. X0EDGEORDER may be either 1 or 2, and simply uses derivatives
%  at that specified order at the minimal x (which may be < 0, for example).
%  The option X0CYLINDRICAL may be true or false. If enabled, the minimal x
%  must be exactly 0; in this case the derivative will be set to 0 there and
%  the second derivative will be computed with central difference by assuming
%  symmetry about the y axis.
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
%   'xi': x index in the main grid
%   'yi': y index in the main grid
%   'xi_sub': x index in the subgrid, if applicable
%   'yi_sub': y index in the subgrid, if applicable
%   'type': NORMAL, SHARED, INT_X_EDGE, INT_Y_EDGE, EXT_X_EDGE, EXT_Y_EDGE, EXT_CORNER
%
%  DX, D2X, DY, and D2Y are differentiation matrices that act on vectors indexed
%  by node number.
%
%  INT_WEIGHT is a vector giving a 2D trapezoidal integration weight for each
%  node. If CYLINDRICALWEIGHT is true, then integration weights will be computed
%  assuming that x is the r coordinate of a 2D cylindrical system.
%
%
% Example:
%
%  Suppose we had a 4x3 coarse grid (marked by O below), with a 3x3 fine grid
%  (marked by *, with shared nodes marked by X).
%
%  O     X  *  X     O
%        *  *  *
%  O     X  *  X     O
%
%  O     O     O     O
%
%  The coarse grid might have x and y vectors [0:3] and [0:2], while the fine
%  grid might be [1 1.5 2] by [1 1.5 2]. Then, to construct differentiation
%  matrices, we could call
%
%    >> coarse_xs = [0 1 2 3];
%    >> coarse_ys = [0 1 2];
%    >> fine_xs = [1 1.5 2];
%    >> fine_ys = [1 1.5 2];
%    >> [idgrids, nodedata, dx, d2x, dy, d2y, int_weight] = ...
%               MULTIMESH_2D( {{coarse_xs,coarse_ys}, {fine_xs, fine_ys}} );
%

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-03-13
% Updated: 2019-05-10

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


function [idgrids, nodedata, dx, d2x, dy, d2y, int_weight] = multimesh_2d(varargin)

    debug = false;
    debug_verbose = false;

    % Define node types
    NORMAL = 0;
    SHARED = 1;
    INT_X_EDGE = 2;
    INT_Y_EDGE = 3;
    EXT_X_EDGE = 5;
    EXT_Y_EDGE = 6;
    EXT_CORNER = 7;

    parser = inputParser;
    parser.addRequired('meshes', @(x) iscell(x));
    parser.addParameter('edgeOrder', 2, @(x) (x == 1 || x == 2));
    parser.addParameter('x0EdgeOrder', NaN, @(x) (x == 1 || x == 2));
    parser.addParameter('x0Cylindrical', 0, @(x) (x == 1 || x == 2));
    parser.addParameter('cylindricalWeight', false, @(x) isscalar(x));
    parser.addParameter('progressBar', false, @(x) isscalar(x));
    parser.addParameter('alignmentTol', 1E-12, @(x) (isscalar(x) && isnumeric(x) && x >= 0));
    parse(parser, varargin{:});
    meshes = parser.Results.meshes;

    edge_order = parser.Results.edgeOrder;
    x0_edge_order = parser.Results.x0EdgeOrder;
    x0_cylindrical = parser.Results.x0Cylindrical;
    if x0_cylindrical
        assert(isnan(x0_edge_order), 'If x0Cylindrical is true, x0EdgeOrder must not be set.');
    end
    if isnan(x0_edge_order); x0_edge_order = edge_order; end;

    cylindrical_weight = parser.Results.cylindricalWeight;

    show_waitbar= parser.Results.progressBar;

    alignment_tol = parser.Results.alignmentTol;

    % Number of meshes
    meshN = length(meshes);
    assert(meshN >= 1, 'At least one input mesh is required.');

    xss = {};
    yss = {};
    xmaxs = [];
    xmins = [];
    ymaxs = [];
    ymins = [];

    % Extract x and y from each mesh
    for j = 1:meshN
        jmesh = meshes{j};
        assert(iscell(jmesh), 'Each mesh must be passed as a cell of form {[x...], [y...]}.');
        assert(length(jmesh)==2, 'Each mesh must be passed as a cell of form {[x...], [y...]}.');
        xss{j} = jmesh{1};
        yss{j} = jmesh{2};
        assert(min([length(xss{j}) length(yss{j})]) > 1, 'Each mesh must contain at least two nodes per side.');
        xmaxs(j) = max(xss{j});
        xmins(j) = min(xss{j});
        ymaxs(j) = max(yss{j});
        ymins(j) = min(yss{j});
    end
    xmin_all = min(xmins);
    if (x0_cylindrical)
        assert(xmin_all == 0, 'To use x0Cylindrical, the grid must have one edge at x=0.');
    end
    xmax_all = max(xmaxs);
    ymin_all = min(ymins);
    ymax_all = max(ymaxs);

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
                        (ymaxs(j) <= ymaxs(k)) && (ymins(j) >= ymins(k)))
                    ancestormat(j,k) = 1;
                    if parentlist(j) == 0
                        parentlist(j) = k;
                    else
                        iparent = parentlist(j);
                        % mesh has multiple parents. Check which is smaller.
                        if ( (xmaxs(k) <= xmaxs(iparent)) && (xmins(k) >= xmins(iparent)) &&...
                                (ymaxs(k) <= ymaxs(iparent)) && (ymins(k) >= ymins(iparent)))
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
                        (ymaxs(j) < ymins(k)) || (ymaxs(k) < ymins(j)) )
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

    for j = 1:meshN

        iparent = parentlist(j);
        if iparent ~= 0 % if not the root mesh
            xs_child = xss{j};
            ys_child = yss{j};
            xs_parent = xss{iparent};
            ys_parent = yss{iparent};

            [xmin_index, xmincheck] = nearest_index(xs_parent, xmins(j));
            [xmax_index, xmaxcheck] = nearest_index(xs_parent, xmaxs(j));
            [ymin_index, ymincheck] = nearest_index(ys_parent, ymins(j));
            [ymax_index, ymaxcheck] = nearest_index(ys_parent, ymaxs(j));
            xmin_parent_indexs(j) = xmin_index;
            xmax_parent_indexs(j) = xmax_index;
            ymin_parent_indexs(j) = ymin_index;
            ymax_parent_indexs(j) = ymax_index;

            % Ensure the corners of each mesh line up with nodes of the parent mesh
            errmsg_submesh_badalign = 'Each submesh must have corners aligning exactly with nodes of its parent mesh.';
            assert(abs(xmincheck - xmins(j)) < alignment_tol, errmsg_submesh_badalign)
            assert(abs(xmaxcheck - xmaxs(j)) < alignment_tol, errmsg_submesh_badalign)
            assert(abs(ymincheck - ymins(j)) < alignment_tol, errmsg_submesh_badalign)
            assert(abs(ymaxcheck - ymaxs(j)) < alignment_tol, errmsg_submesh_badalign)

            % Ensure that the mesh is not zero width or height
            errmsg_submesh_onespace = 'Each submesh must have nonzero width and height.';
            assert(max(xs_child) > min(xs_child), errmsg_submesh_onespace);
            assert(max(ys_child) > min(ys_child), errmsg_submesh_onespace);
            assert(xmin_index < xmax_index, errmsg_submesh_onespace);
            assert(ymin_index < ymax_index, errmsg_submesh_onespace);

            % Ensure that the mesh spacing is above the tolerance
            errmsg_meshspacetol = sprintf('Mesh spacing must be greater than alignment tolerance (%g). Or, change the alignment tolerance.', alignment_tol) ;
            assert(min(abs(diff(xs_child)))>alignment_tol, errmsg_meshspacetol);
            assert(min(abs(diff(ys_child)))>alignment_tol, errmsg_meshspacetol);

            % Ensure that the mesh has at least 3 nodes in each direction
            errmsg_submesh_threespace = 'Each submesh must have at least three nodes in width and height.';
            assert(length(xs_child) >= 3, errmsg_submesh_threespace);
            assert(length(ys_child) >= 3, errmsg_submesh_threespace);

            % Ensure that the submesh spacing divides the parent mesh spacing
            errmsg_submesh_divisible = 'The submesh spacing must divide the parent mesh spacing.';
            assert( ((xmax_index - xmin_index) < length(xs_child)-1), errmsg_submesh_divisible );
            assert( ((ymax_index - ymin_index) < length(ys_child)-1), errmsg_submesh_divisible );
            for xj = xmin_index+1:xmax_index-1
                [xcheck_index, xcheck] = nearest_index(xs_child, xs_parent(xj));
                assert(abs(xcheck - xs_child(xcheck_index)) < alignment_tol, errmsg_submesh_divisible)
            end
            for yj = ymin_index+1:ymax_index-1
                [ycheck_index, ycheck] = nearest_index(ys_child, ys_parent(yj));
                assert(abs(ycheck - ys_child(ycheck_index)) < alignment_tol, errmsg_submesh_divisible)
            end

            % Ensure that the submesh is at least one full space inside the parent mesh
            errmsg_submesh_padding = 'Each submesh must be one full mesh space inside the edges of its parent mesh, unless it is at an edge.';
            assert( (xmin_index > 1) || (xs_child(1) == xmin_all), errmsg_submesh_padding );
            assert( (xmax_index < length(xss{iparent})) || (xs_child(end) == xmax_all), errmsg_submesh_padding);
            assert( (ymin_index > 1) || (ys_child(1) == ymin_all), errmsg_submesh_padding );
            assert( (ymax_index < length(yss{iparent})) || (ys_child(end) == ymax_all), errmsg_submesh_padding);
        end
    end

    if (debug); disp('Constructing nodedata...'); end
    % Assign id numbers to each of the nodes. We start at the largest mesh and work downwards; shared nodes will take the value of their parent mesh, not the child.
    sort_mesh_i = sorttree(parentlist);

    % Nodedata contains data for each node. Nodedata has the following structure, each of which is indexed by node.
    %   'grid': the main (parent) grid to which the node belongs
    %   'subgrid': the sub (child) grid to which the node belongs, if applicable
    %   'x': x value of node
    %   'y': y value of node
    %   'xi': x index in the main grid
    %   'yi': y index in the main grid
    %   'xi_sub': x index in the subgrid, if applicable
    %   'yi_sub': y index in the subgrid, if applicable
    %   'type': NORMAL, SHARED, INT_X_EDGE, INT_Y_EDGE, EXT_X_EDGE, EXT_Y_EDGE, EXT_CORNER
    nodedata.grid = [];
    nodedata.subgrid = [];
    nodedata.x = [];
    nodedata.y = [];
    nodedata.xi = [];
    nodedata.yi = [];
    nodedata.xi_sub = [];
    nodedata.yi_sub = [];
    nodedata.type = [];

    cur_nodeid = 0;
    idgrids = {};
    for junsorted = 1:meshN
        j = sort_mesh_i(junsorted);

        iparent = parentlist(j);
        xs_self = xss{j};
        ys_self = yss{j};

        idgrid = zeros(length(xs_self), length(ys_self));

        for yi = 1:length(ys_self)
            for xi = 1:length(xs_self)
                x = xs_self(xi);
                y = ys_self(yi);

                % determine if node is shared by a parent
                node_is_shared = false;
                if (iparent ~= 0)
                    xs_parent = xss{iparent};
                    ys_parent = yss{iparent};
                    [xi_parent, x_parent] = nearest_index(xs_parent, x);
                    [yi_parent, y_parent] = nearest_index(ys_parent, y);
                    if ( (abs(x_parent- x) < alignment_tol) && (abs(y_parent - y) < alignment_tol))
                        node_is_shared = true;
                        idgrid_parent = idgrids{iparent};
                        parent_nodeid = idgrid_parent(xi_parent, yi_parent);
                    end
                end


                % if node is not shared, then a new entry is needed.
                if ~node_is_shared
                    cur_nodeid = cur_nodeid + 1;
                    idgrid(xi, yi) = cur_nodeid;
                    nodedata.grid(cur_nodeid) = j;
                    nodedata.subgrid(cur_nodeid) = 0;
                    nodedata.x(cur_nodeid) = x;
                    nodedata.y(cur_nodeid) = y;
                    nodedata.xi(cur_nodeid) = xi;
                    nodedata.yi(cur_nodeid) = yi;
                    nodedata.xi_sub(cur_nodeid) = 0;
                    nodedata.yi_sub(cur_nodeid) = 0;
                    if ( (x == xmax_all) || (x == xmin_all) ) && ( (y == ymax_all) || (y == ymin_all))
                        nodedata.type(cur_nodeid) = EXT_CORNER;
                    elseif ( (x == xmax_all) || (x == xmin_all) )
                        nodedata.type(cur_nodeid) = EXT_X_EDGE;
                    elseif ( (y == ymax_all) || (y == ymin_all) )
                        nodedata.type(cur_nodeid) = EXT_Y_EDGE;
                    elseif ( (xi == length(xs_self)) || (xi == 1) ) && ( (yi == length(ys_self)) || (yi == 1))
                        % "Interior corner"
                        assert(false); % This should never occur, since it should always be a shared node.
                    elseif ( (xi == length(xs_self)) || (xi == 1) )
                        nodedata.type(cur_nodeid) = INT_X_EDGE;
                    elseif ( (yi == length(ys_self)) || (yi == 1))
                        nodedata.type(cur_nodeid) = INT_Y_EDGE;
                    else
                        nodedata.type(cur_nodeid) = NORMAL;
                    end
                else
                    % if node is shared, then it has already been listed, but we
                    % still must repopulate it with the correct parent information.
                    assert(iparent ~= 0)
                    idgrid(xi, yi) = parent_nodeid;
                    nodedata.grid(parent_nodeid) = iparent;
                    nodedata.x(parent_nodeid) = x_parent;
                    nodedata.y(parent_nodeid) = y_parent;
                    nodedata.xi(parent_nodeid) = xi_parent;
                    nodedata.yi(parent_nodeid) = yi_parent;
                    nodedata.subgrid(parent_nodeid) = j;
                    nodedata.xi_sub(parent_nodeid) = xi;
                    nodedata.yi_sub(parent_nodeid) = yi;
                    nodedata.type(parent_nodeid) = SHARED;

                    if ( (x == xmax_all) || (x == xmin_all) ) && ( (y == ymax_all) || (y == ymin_all))
                        nodedata.type(parent_nodeid) = EXT_CORNER;
                    elseif ( (x == xmax_all) || (x == xmin_all) )
                        nodedata.type(parent_nodeid) = EXT_X_EDGE;
                    elseif ( (y == ymax_all) || (y == ymin_all) )
                        nodedata.type(parent_nodeid) = EXT_Y_EDGE;
                    end
                end


            end
        end

        idgrids{j} = idgrid;
    end
    nodeN = length(nodedata.x);
    assert(nodeN==cur_nodeid);


    % Set up neighbor fetching shortcuts
    get_left_node = @(node) mm2d_get_neighbor_node(node, 'left', nodedata, idgrids);
    get_right_node = @(node) mm2d_get_neighbor_node(node, 'right', nodedata, idgrids);
    get_down_node = @(node) mm2d_get_neighbor_node(node, 'down', nodedata, idgrids);
    get_up_node = @(node) mm2d_get_neighbor_node(node, 'up', nodedata, idgrids);
    get_left_parent_node = @(node) mm2d_get_neighbor_parent_node(node, 'left', nodedata, idgrids);
    get_right_parent_node = @(node) mm2d_get_neighbor_parent_node(node, 'right', nodedata, idgrids);
    get_down_parent_node = @(node) mm2d_get_neighbor_parent_node(node, 'down', nodedata, idgrids);
    get_up_parent_node = @(node) mm2d_get_neighbor_parent_node(node, 'up', nodedata, idgrids);

    if (debug_verbose);
        nodedata
        for j = 1:meshN
            idgrids{j}
        end



        % check neighbor code
        for this_node = 1:nodeN
            this_xval = nodedata.x(this_node);
            this_yval = nodedata.y(this_node);

            this_left_node = get_left_node(this_node);
            this_left_xval = NaN;
            this_left_yval = NaN;
            if (this_left_node > 0)
                this_left_xval = nodedata.x(this_left_node);
                this_left_yval = nodedata.y(this_left_node);
                assert(this_left_xval < this_xval);
                assert(abs(this_left_yval - this_yval) < alignment_tol);
            end

            this_right_node = get_right_node(this_node);
            this_right_xval = NaN;
            this_right_yval = NaN;
            if (this_right_node > 0)
                this_right_xval = nodedata.x(this_right_node);
                this_right_yval = nodedata.y(this_right_node);
                assert(this_right_xval > this_xval);
                assert(abs(this_right_yval - this_yval) < alignment_tol);
            end

            this_down_node = get_down_node(this_node);
            this_down_xval = NaN;
            this_down_yval = NaN;
            if (this_down_node > 0)
                this_down_xval = nodedata.x(this_down_node);
                this_down_yval = nodedata.y(this_down_node);
                assert(abs(this_down_xval - this_xval) < alignment_tol);
                assert(this_down_yval < this_yval);
            end

            this_up_node = get_up_node(this_node);
            this_up_xval = NaN;
            this_up_yval = NaN;
            if (this_up_node > 0)
                this_up_xval = nodedata.x(this_up_node);
                this_up_yval = nodedata.y(this_up_node);
                assert(abs(this_up_xval - this_xval) < alignment_tol);
                assert(this_up_yval > this_yval);
            end

            this_left_parent_node = mm2d_get_neighbor_parent_node(this_node, 'left', nodedata, idgrids);
            this_right_parent_node = mm2d_get_neighbor_parent_node(this_node, 'right', nodedata, idgrids);
            this_down_parent_node = mm2d_get_neighbor_parent_node(this_node, 'down', nodedata, idgrids);
            this_up_parent_node = mm2d_get_neighbor_parent_node(this_node, 'up', nodedata, idgrids);
            %disp(sprintf('Type: %d',nodedata.type(this_node)))
            %disp(sprintf('[grid,subgrid]: [%d,%d]',nodedata.grid(this_node), nodedata.subgrid(this_node)))
            %disp(sprintf('LCR: %d (%d) %d', this_left_node, this_node, this_right_node))
            %disp(sprintf('DCU: %d (%d) %d', this_down_node, this_node, this_up_node))
            %disp(sprintf('LCRp: %d (%d) %d', this_left_parent_node, this_node, this_right_parent_node))
            %disp(sprintf('DCUp: %d (%d) %d', this_down_parent_node, this_node, this_up_parent_node))
            %disp(' ');
        end
    end




    %%%%%%%%%% SET UP DIFFERENTIATION MATRICES %%%%%%%%%%

    % make sparse zero matrices to start with.
    dx = sparse(nodeN, nodeN, 0);
    d2x = sparse(nodeN, nodeN, 0);
    dy = sparse(nodeN, nodeN, 0);
    d2y = sparse(nodeN, nodeN, 0);

    % First, efficiently set up differences on all the interior "NORMAL" nodes.
    % Next, go back and deal with the edge cases one by one.

    all_nodes = [1:nodeN];
    nontrivial_nodes_dx = all_nodes( (nodedata.type ~= NORMAL) & (nodedata.type ~= INT_Y_EDGE) & ~( (nodedata.type == EXT_Y_EDGE) & (nodedata.subgrid == 0)) );
    nontrivial_nodes_dy = all_nodes( (nodedata.type ~= NORMAL) & (nodedata.type ~= INT_X_EDGE) & ~( (nodedata.type == EXT_X_EDGE) & (nodedata.subgrid == 0)) );

    % construct dx and dy for each mesh with deriv_findiff
    for j = 1:meshN
        this_xs = xss{j};
        this_ys = yss{j};
        this_idgrid = idgrids{j};

        [this_dx, this_d2x, ~, ~] = deriv_findiff(this_xs, 2);
        [this_dy, this_d2y, ~, ~] = deriv_findiff(this_ys, 2);

        % permute results to match node ids...
        for k = 1:length(this_ys)
            dx(this_idgrid(:,k), this_idgrid(:,k)) = this_dx;
            d2x(this_idgrid(:,k), this_idgrid(:,k)) = this_d2x;
        end
        for k = 1:length(this_xs)
            dy(this_idgrid(k,:), this_idgrid(k,:)) = this_dy;
            d2y(this_idgrid(k,:), this_idgrid(k,:)) = this_d2y;
        end

    end

    % Construct dx and d2x edge cases
    if (debug); disp('Constructing dx & d2x...'); end
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
        this_subgridindex = nodedata.subgrid(this_node);
        this_xi_sub       = nodedata.xi_sub(this_node);
        this_yi_sub       = nodedata.yi_sub(this_node);
        this_x            = nodedata.x(this_node);
        this_y            = nodedata.y(this_node);

        if (this_subgridindex == 0)
            this_subidgrid = [];
        else
            this_subidgrid = idgrids{this_subgridindex};
        end

        % clear this row of dx / d2x
        dx(this_node, :) = 0;
        d2x(this_node, :) = 0;

        if (this_type == NORMAL || this_type == INT_Y_EDGE || this_type == EXT_Y_EDGE || this_type == SHARED)
            % Easiest case, just do central difference

            left_node = get_left_node(this_node);
            right_node = get_right_node(this_node);
            assert((left_node ~= 0) && (right_node ~= 0));
            left_x = nodedata.x(left_node);
            right_x = nodedata.x(right_node);

            h3 = (right_x - this_x);
            h4 = (this_x - left_x);
            dx(this_node, left_node) = -(h3/h4)/(h3 + h4);
            dx(this_node, this_node) = 1/h4 - 1/h3;
            dx(this_node, right_node) = (h4/h3)/(h3 + h4);

            d2x(this_node, left_node) = 2/(h4*(h3 + h4));
            d2x(this_node, this_node) = -2/(h3*h4);
            d2x(this_node, right_node) = 2/(h3*(h3 + h4));
        elseif (this_type == EXT_X_EDGE || this_type == EXT_CORNER)
            % We're at the true edge of the whole mesh, we'll just do a first order difference.

            % use first or second order difference, depending on edge_order, x0_edge_order, and x0_cylindrical
            if (this_xi == 1)
                next_node = get_right_node(this_node);
                far_node = get_right_node(next_node);
            elseif (this_xi == size(this_idgrid,1));
                next_node = get_left_node(this_node);
                far_node = get_left_node(next_node);
            else
                assert(false)
            end
            assert((next_node ~= 0) && (far_node ~= 0));
            h3 = nodedata.x(next_node) - nodedata.x(this_node);
            h5 = nodedata.x(far_node) - nodedata.x(next_node);

            if ( (x0_cylindrical) && (abs(this_x - xmin_all) < alignment_tol) )
                dx(this_node, this_node) = 0;
                dx(this_node, next_node) = 0;

                d2x(this_node, this_node) = -2/(h3^2);
                d2x(this_node, next_node) = 2/(h3^2);
            elseif (edge_order == 1) || ( (x0_edge_order == 1) && (abs(this_x - xmin_all) < alignment_tol) )
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

        elseif (this_type == INT_X_EDGE)
            % We're at a mesh boundary within the system; we use 6 nearby points to construct
            % the correct second order difference (time for the ghost nodes).

            assert(this_subgridindex == 0);
            assert(this_gridindex ~= 0);

            % 1) get nearest neighbors above and below
            up_node = get_up_parent_node(this_node);
            down_node = get_down_parent_node(this_node);
            assert( (up_node ~= 0) && (down_node ~= 0) );
            assert(nodedata.grid(down_node) ~= this_gridindex);

            % 2) get side neighbors
            if (this_xi == 1)
                inside_node = get_right_node(this_node);
                outside_up_node = get_left_node(up_node);
                outside_down_node = get_left_node(down_node);
            elseif (this_xi == size(this_idgrid,1));
                inside_node = get_left_node(this_node);
                outside_up_node = get_right_node(up_node);
                outside_down_node = get_right_node(down_node);
            end
            assert( (outside_up_node ~= 0) && (outside_down_node ~= 0) );

            h1 = nodedata.y(up_node) - nodedata.y(this_node);
            h2 = nodedata.y(this_node) - nodedata.y(down_node);
            h3 = nodedata.x(outside_up_node) - nodedata.x(this_node);
            assert(abs(h3-(nodedata.x(outside_down_node) - nodedata.x(this_node)))<0.5*alignment_tol);
            h4 = nodedata.x(this_node) - nodedata.x(inside_node);

            dx(this_node, this_node) = h3/(h4*(h3+h4));
            dx(this_node, up_node) = -h2*h4/( (h1+h2)*h3*(h3+h4) );
            dx(this_node, down_node) = -h1*h4/( (h1+h2)*h3*(h3+h4) );
            dx(this_node, outside_up_node) = h2*h4/( (h1+h2)*h3*(h3+h4) );
            dx(this_node, outside_down_node) = h1*h4/( (h1+h2)*h3*(h3+h4) );
            dx(this_node, inside_node) = -h3/(h4*(h3+h4));

            d2x(this_node, this_node) = -2/(h4*(h3+h4));
            d2x(this_node, up_node) = -2*h2/( (h1+h2)*h3*(h3+h4) );
            d2x(this_node, down_node) = -2*h1/( (h1+h2)*h3*(h3+h4) );
            d2x(this_node, outside_up_node) = 2*h2/( (h1+h2)*h3*(h3+h4) );
            d2x(this_node, outside_down_node) = 2*h1/( (h1+h2)*h3*(h3+h4) );
            d2x(this_node, inside_node) = 2/(h4*(h3+h4));
        end

        if (show_waitbar)
            waitbar(progress/length(nontrivial_nodes_dx), f_wait, sprintf('Constructing dx & d2x operators...'));
        end
    end

    if (debug); disp('Constructing dy & d2y...'); end
    progress = 0;
    if (show_waitbar)
        waitbar(progress/length(nontrivial_nodes_dy), f_wait, sprintf('Constructing dy & d2y operators...'));
    end
    for this_node = nontrivial_nodes_dy
        progress = progress + 1;
        this_type         = nodedata.type(this_node);
        this_gridindex    = nodedata.grid(this_node);
        this_idgrid       = idgrids{this_gridindex};
        this_xi           = nodedata.xi(this_node);
        this_yi           = nodedata.yi(this_node);
        this_subgridindex = nodedata.subgrid(this_node);
        this_xi_sub       = nodedata.xi_sub(this_node);
        this_yi_sub       = nodedata.yi_sub(this_node);
        this_x            = nodedata.x(this_node);
        this_y            = nodedata.y(this_node);

        if (this_subgridindex == 0)
            this_subidgrid = [];
        else
            this_subidgrid = idgrids{this_subgridindex};
        end

        % clear this row of dy / d2y
        dy(this_node, :) = 0;
        d2y(this_node, :) = 0;

        if (this_type == NORMAL || this_type == INT_X_EDGE || this_type == EXT_X_EDGE || this_type == SHARED)
            % Easiest case, just do central difference

            down_node = get_down_node(this_node);
            up_node = get_up_node(this_node);
            assert((down_node ~= 0) && (up_node ~= 0));
            down_y = nodedata.y(down_node);
            up_y = nodedata.y(up_node);

            h3 = (up_y - this_y);
            h4 = (this_y - down_y);
            dy(this_node, down_node) = -(h3/h4)/(h3 + h4);
            dy(this_node, this_node) = 1/h4 - 1/h3;
            dy(this_node, up_node) = (h4/h3)/(h3 + h4);

            d2y(this_node, down_node) = 2/(h4*(h3 + h4));
            d2y(this_node, this_node) = -2/(h3*h4);
            d2y(this_node, up_node) = 2/(h3*(h3 + h4));
        elseif (this_type == EXT_Y_EDGE || this_type == EXT_CORNER)
            % We're at the true edge of the whole mesh, we'll just do a first order difference.

            % use first order difference
            % TODO: come up with a better scheme?
            if (this_yi == 1)
                next_node = get_up_node(this_node);
                far_node = get_up_node(next_node);
            elseif (this_yi == size(this_idgrid,2));
                next_node = get_down_node(this_node);
                far_node = get_down_node(next_node);
            else
                assert(false)
            end
            assert((next_node ~= 0) && (far_node ~= 0));
            h3 = nodedata.y(next_node) - nodedata.y(this_node);
            h5 = nodedata.y(far_node) - nodedata.y(next_node);

            if (edge_order == 1)
                dy(this_node, this_node) = -1/h3;
                dy(this_node, next_node) = 1/h3;
            else
                dy(this_node, this_node) = -(2*h3+h5)/(h3*(h3+h5));
                dy(this_node, next_node) = (1/h3 + 1/h5);
                dy(this_node, far_node) = -h3/(h5*(h3+h5));
            end

            d2y(this_node, this_node) = 2/(h3*(h3+h5));
            d2y(this_node, next_node) = -2/(h3*h5);
            d2y(this_node, far_node) = 2/(h5*(h3+h5));
        elseif (this_type == INT_Y_EDGE)
            % We're at a mesh boundary within the system; we use 6 nearby points to construct
            % the correct second order difference (time for the ghost nodes).

            assert(this_subgridindex == 0);
            assert(this_gridindex ~= 0);

            % 1) get nearest neighbors to the left and right
            right_node = get_right_parent_node(this_node);
            left_node = get_left_parent_node(this_node);
            assert( (right_node ~= 0) && (left_node ~= 0) );
            assert(nodedata.grid(left_node) ~= this_gridindex);

            % 2) get up and down neighbors
            if (this_yi == 1)
                inside_node = get_up_node(this_node);
                outside_right_node = get_down_node(right_node);
                outside_left_node = get_down_node(left_node);
            elseif (this_yi == size(this_idgrid,2));
                inside_node = get_down_node(this_node);
                outside_right_node = get_up_node(right_node);
                outside_left_node = get_up_node(left_node);
            end
            assert( (outside_right_node ~= 0) && (outside_left_node ~= 0) );

            h1 = nodedata.x(right_node) - nodedata.x(this_node);
            h2 = nodedata.x(this_node) - nodedata.x(left_node);
            h3 = nodedata.y(outside_right_node) - nodedata.y(this_node);
            assert(abs(h3-(nodedata.y(outside_left_node) - nodedata.y(this_node)))<0.5*alignment_tol);
            h4 = nodedata.y(this_node) - nodedata.y(inside_node);

            dy(this_node, this_node) = h3/(h4*(h3+h4));
            dy(this_node, right_node) = -h2*h4/( (h1+h2)*h3*(h3+h4) );
            dy(this_node, left_node) = -h1*h4/( (h1+h2)*h3*(h3+h4) );
            dy(this_node, outside_right_node) = h2*h4/( (h1+h2)*h3*(h3+h4) );
            dy(this_node, outside_left_node) = h1*h4/( (h1+h2)*h3*(h3+h4) );
            dy(this_node, inside_node) = -h3/(h4*(h3+h4));

            d2y(this_node, this_node) = -2/(h4*(h3+h4));
            d2y(this_node, right_node) = -2*h2/( (h1+h2)*h3*(h3+h4) );
            d2y(this_node, left_node) = -2*h1/( (h1+h2)*h3*(h3+h4) );
            d2y(this_node, outside_right_node) = 2*h2/( (h1+h2)*h3*(h3+h4) );
            d2y(this_node, outside_left_node) = 2*h1/( (h1+h2)*h3*(h3+h4) );
            d2y(this_node, inside_node) = 2/(h4*(h3+h4));
        end
        if (show_waitbar)
            waitbar(progress/length(nontrivial_nodes_dy), f_wait, sprintf('Constructing dy & d2y operators...'));
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
        this_idgrid = idgrids{j};

        this_x_distL = [0 diff(this_xs)];
        this_x_distR = [diff(this_xs) 0];
        this_y_height = 0.5*([0 diff(this_ys)] + [diff(this_ys) 0]);

        for xi = 1:length(this_xs)
            for yi = 1:length(this_ys)
                this_node = this_idgrid(xi, yi);

                if cylindrical_weight
                    this_x = nodedata.x(this_node);
                    distL = this_x_distL(xi);
                    distR = this_x_distR(xi);
                    height = this_y_height(yi);
                    width = 0.5*(distL + distR);

                    % Any parts with r < 0 should be omitted
                    if this_x < 0
                        if this_x + distR < 0
                            % if a full node below 0, just ignore it
                            int_weight(this_node) = 0;
                        else
                            % If the node to the right is above 0, need to be careful
                            int_weight(this_node) = height*pi*(this_x + distR)^3/(3*distR);
                            % TODO: CHECK THIS!!!!
                        end
                    elseif this_x - distL < 0
                        % If this node is above zero but the node to the left is below 0, need to be careful again
                        int_weight(this_node) = height*pi*( this_x^2 - (this_x)^3/(3*distL) + this_x*distR + distR^2/3);
                        % TODO: CHECK THIS!!!!
                    else
                        int_weight(this_node) = width*height*2*pi*(this_x + distR/3 - distL/3);
                    end
                else
                    int_weight(this_node) = this_y_height(yi)*( this_x_distL(xi) + this_x_distR(xi) )*0.5;
                end

            end
        end

    end


    progress = 0;
    if (show_waitbar)
        waitbar(progress/length(nontrivial_nodes_intwgt), f_wait, sprintf('Calculating integration weights...'));
    end

    for this_node = nontrivial_nodes_intwgt
        progress = progress + 1;

        left_node = get_left_node(this_node);
        if (left_node == 0)
            up_parent_node = get_up_parent_node(this_node);
            if (up_parent_node ~= 0)
                left_node = get_left_node(up_parent_node);
            end
        end

        right_node = get_right_node(this_node);
        if (right_node == 0)
            up_parent_node = get_up_parent_node(this_node);
            if (up_parent_node ~= 0)
                right_node = get_right_node(up_parent_node);
            end
        end

        down_node = get_down_node(this_node);
        if (down_node == 0)
            right_parent_node = get_right_parent_node(this_node);
            if (right_parent_node ~= 0)
                down_node = get_down_node(right_parent_node);
            end
        end

        up_node = get_up_node(this_node);
        if (up_node == 0)
            right_parent_node = get_right_parent_node(this_node);
            if (right_parent_node ~= 0)
                up_node = get_up_node(right_parent_node);
            end
        end


        distL = 0;
        distR = 0;
        if (left_node ~= 0)
            distL = (nodedata.x(this_node) - nodedata.x(left_node));
        end
        if (right_node ~= 0)
            distR = (nodedata.x(right_node) - nodedata.x(this_node));
        end
        width = 0.5*(distR + distL);
        assert(width > 0);


        height = 0;
        if (down_node ~= 0)
            height = height + 0.5*(nodedata.y(this_node) - nodedata.y(down_node));
        end
        if (up_node ~= 0)
            height = height + 0.5*(nodedata.y(up_node) - nodedata.y(this_node));
        end
        assert(height > 0);


        if cylindrical_weight
            this_x = nodedata.x(this_node);
            int_weight(this_node) = width*height*2*pi*(this_x + distR/3 - distL/3);

            % Any parts with r < 0 should be omitted
            if this_x < 0
                if this_x + distR < 0
                    % if a full node below 0, just ignore it
                    int_weight(this_node) = 0;
                else
                    % If the node to the right is above 0, need to be careful
                    int_weight(this_node) = height*pi*(this_x + distR)^3/(3*distR);
                    % TODO: CHECK THIS!!!!
                end
            elseif this_x - distL < 0
                % If this node is above zero but the node to the left is below 0, need to be careful again
                int_weight(this_node) = height*pi*( this_x^2 - (this_x)^3/(3*distL) + this_x*distR + distR^2/3);
                % TODO: CHECK THIS!!!!
            end
        else
            int_weight(this_node) = width*height;
        end

        if (show_waitbar)
            waitbar(progress/length(nontrivial_nodes_intwgt), f_wait, sprintf('Calculating integration weights...'));
        end
    end

    if (show_waitbar)
        delete(f_wait);
    end

end



%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%

function next_node = mm2d_get_neighbor_node(this_node, direction, nodedata, idgrids)
    % Direction: 'left', 'right', 'up', 'down'

    % Define node types
    NORMAL = 0;
    SHARED = 1;
    INT_X_EDGE = 2;
    INT_Y_EDGE = 3;
    EXT_X_EDGE = 5;
    EXT_Y_EDGE = 6;
    EXT_CORNER = 7;

    this_type         = nodedata.type(this_node);
    this_gridindex    = nodedata.grid(this_node);
    this_idgrid       = idgrids{this_gridindex};
    this_xi           = nodedata.xi(this_node);
    this_yi           = nodedata.yi(this_node);
    this_subgridindex = nodedata.subgrid(this_node);
    this_xi_sub       = nodedata.xi_sub(this_node);
    this_yi_sub       = nodedata.yi_sub(this_node);
    this_x            = nodedata.x(this_node);
    this_y            = nodedata.y(this_node);

    if strcmpi(direction,'left')
        not_on_edge = (this_xi > 1);
        sub_not_on_edge = (this_xi_sub > 1);
        xistep = -1;
        yistep = 0;
    elseif strcmpi(direction,'right')
        not_on_edge = (this_xi < size(this_idgrid,1));
        if (this_subgridindex ~= 0)
            sub_not_on_edge = (this_xi_sub < size(idgrids{this_subgridindex},1));
        end
        xistep = 1;
        yistep = 0;
    elseif strcmpi(direction,'down')
        not_on_edge = (this_yi > 1);
        sub_not_on_edge = (this_yi_sub > 1);
        xistep = 0;
        yistep = -1;
    elseif strcmpi(direction,'up')
        not_on_edge = (this_yi < size(this_idgrid,2));
        if (this_subgridindex ~= 0)
            sub_not_on_edge = (this_yi_sub < size(idgrids{this_subgridindex},2));
        end
        xistep = 0;
        yistep = 1;
    else
        assert('direction must be ''left'', ''right'', ''up'', or ''down''.');
    end

    if this_subgridindex == 0 % no subgrid
        if not_on_edge
            next_node = this_idgrid(this_xi + xistep, this_yi + yistep);
        else
            next_node = 0; % true edge
        end
    else % yes subgrid
        if sub_not_on_edge
            next_node = idgrids{this_subgridindex}(this_xi_sub + xistep, this_yi_sub + yistep);
        else
            % we're at an interior boundary; we can look at the parent mesh
            if not_on_edge
                next_node = this_idgrid(this_xi + xistep, this_yi + yistep);
            else
                next_node = 0; % edge of parent mesh too
            end
        end
    end

end

function next_parent_node = mm2d_get_neighbor_parent_node(this_node, direction, nodedata, idgrids)
    % Direction: 'left', 'right', 'up', 'down'

    % Define node types
    NORMAL = 0;
    SHARED = 1;
    INT_X_EDGE = 2;
    INT_Y_EDGE = 3;
    EXT_X_EDGE = 5;
    EXT_Y_EDGE = 6;
    EXT_CORNER = 7;

    this_type         = nodedata.type(this_node);
    this_gridindex    = nodedata.grid(this_node);
    this_idgrid       = idgrids{this_gridindex};
    this_xi           = nodedata.xi(this_node);
    this_yi           = nodedata.yi(this_node);
    this_subgridindex = nodedata.subgrid(this_node);
    this_xi_sub       = nodedata.xi_sub(this_node);
    this_yi_sub       = nodedata.yi_sub(this_node);
    this_x            = nodedata.x(this_node);
    this_y            = nodedata.y(this_node);

    next_parent_node = 0;

    if (this_gridindex == 0) && (this_subgridindex == 0)
        % Fail if this is the root mesh and not a shared node
        return
    end

    if strcmpi(direction,'left')
        not_on_edge = (this_xi > 1);
        sub_not_on_edge = (this_xi_sub > 1);
        xistep = -1;
        yistep = 0;
    elseif strcmpi(direction,'right')
        not_on_edge = (this_xi < size(this_idgrid,1));
        if (this_subgridindex ~= 0)
            sub_not_on_edge = (this_xi_sub < size(idgrids{this_subgridindex},1));
        end
        xistep = 1;
        yistep = 0;
    elseif strcmpi(direction,'down')
        not_on_edge = (this_yi > 1);
        sub_not_on_edge = (this_yi_sub > 1);
        xistep = 0;
        yistep = -1;
    elseif strcmpi(direction,'up')
        not_on_edge = (this_yi < size(this_idgrid,2));
        if (this_subgridindex ~= 0)
            sub_not_on_edge = (this_yi_sub < size(idgrids{this_subgridindex},2));
        end
        xistep = 0;
        yistep = 1;
    else
        assert('direction must be ''left'', ''right'', ''up'', or ''down''.');
    end

    %if this_type == SHARED
    if this_subgridindex ~= 0
        % our goal is to find another node with the same idgrid as the parent of this node
        goal = @(node) ((node == 0) || nodedata.grid(node) == this_gridindex || nodedata.subgrid(node) == this_gridindex);
    else
        % our goal is to find another node whose child is the same idgrid as this node
        goal = @(node) ((node == 0) || nodedata.subgrid(node) == this_gridindex);
    end

    next_node = mm2d_get_neighbor_node(this_node, direction, nodedata, idgrids);
    while ~goal(next_node)
        next_node = mm2d_get_neighbor_node(next_node, direction, nodedata, idgrids);
    end
    next_parent_node = next_node;

end
