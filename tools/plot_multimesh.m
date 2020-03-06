% PLOT_MULTIMESH Construct a plot of the nested mesh boundaries.
%
% PLOT_MULTIMESH(meshes)
% PLOT_MULTIMESH(ax, meshes)
% PLOT_MULTIMESH(ax, meshes,...
%       'PointColor',          point_color,   ...
%       'PointSize',           point_size,    ...
%       'EdgeColor',           edge_color,    ...
%       'BodyCoordinates',     body_coords,   ...
%       'BodyRadii',           body_radii,    ...
%       'BodyColor',           body_color,    ...
%       )
%
% TODO: instructions

function plot_multimesh(varargin)

hasax = isa(varargin{1},'matlab.graphics.axis.Axes');



parser = inputParser;
if (hasax)
    parser.addRequired('axes');
end
parser.addRequired('meshes', @(x) iscell(x));
parser.addParameter('PointColor', 'k', @(x) ischar(x));
parser.addParameter('PointSize', 20, @(x) isnumeric(x) && isscalar(x) && ispositive(x));
parser.addParameter('EdgeColor', 'r', @(x) ischar(x));
parser.addParameter('BodyCoordinates', [], @(x) isnumeric(x));
parser.addParameter('BodyRadii', [], @(x) isnumeric(x) && (isvector(x) || isempty(x)));
parser.addParameter('BodyColor', 'b', @(x) ischar(x) || iscell(x));
parse(parser, varargin{:});
meshes = parser.Results.meshes;
point_color = parser.Results.PointColor;
point_size = parser.Results.PointSize;
edge_color = parser.Results.EdgeColor;
body_coords = parser.Results.BodyCoordinates;
body_radii = parser.Results.BodyRadii;
body_color_temp = parser.Results.BodyColor;

Nbody = length(body_radii);

hasbodies = ~(isempty(body_coords) && isempty(body_radii));

if iscell(body_color_temp)
    body_colors = body_color_temp;
    assert(length(body_colors) == Nbody, 'If body_colors is passed as a cell of colors, it must be the same length as body_radii and body_coords.');
else
    body_colors = {};
    for j = 1:Nbody
        body_colors{end+1} = body_color_temp;
    end
end

Nmesh = length(meshes);
assert(Nmesh > 0);
dim = length(meshes{1});
assert(dim == 2 || dim == 3);

if hasbodies
    if (dim == 3)
        assert(size(body_coords, 1) == Nbody, ...
            'body_coords and body_radii must have the same length (one entry for each body)');
        assert(size(body_coords, 2) == 3);
    elseif (dim == 2)
        if (min(size(body_coords)) == 1) && (length(body_coords) == Nbody)
            % only z coordinates were passed
            if size(body_coords,1) == 1
                body_coords = body_coords';
            end
            body_coords_new = zeros(Nbody, 2);
            body_coords_new(:,2) = body_coords;
            body_coords = body_coords_new;
        elseif (size(body_coords,1) == Nbody) && (size(body_coords,2) == 2)
            % r,z coordinates passed already
        else
            assert('body_coords should be either a list of z coordinates or a table of [r,z] coordinates');
        end
    else
        assert(false);
    end
end

if (hasax)
    ax = parser.Results.axes;
else
    ax = gca;
end

holding = ishold(ax);

if (dim == 2)
    % 2D
    rmin = 0;
    rmax = 0;
    zmin = Inf;
    zmax = -Inf;
    for j = 1:Nmesh
        this_rss = meshes{j}{1};
        this_zss = meshes{j}{2};
        this_rbound = [this_rss(1) this_rss(end)];
        this_zbound = [this_zss(1) this_zss(end)];
        rmax = max(rmax,max(this_rss));
        zmin = min(zmin,min(this_zss));
        zmax = max(zmax,max(this_zss));

        for i1 = 1:2; for i2 = 1:2;
            plot(ax, [this_rbound(1) this_rbound(2)], [this_zbound(i2) this_zbound(i2)], edge_color);
            hold(ax, 'on')
            plot(ax, [this_rbound(i1) this_rbound(i1)], [this_zbound(i2) this_zbound(i2)], edge_color);
            plot(ax, [this_rbound(i1) this_rbound(i1)], [this_zbound(1) this_zbound(2)], edge_color);
        end; end
        if ~strcmpi(point_color,'none')
            [rpoints, zpoints] = ndgrid(meshes{j}{1},meshes{j}{2});
            rpoints = reshape(rpoints, 1, prod(size(rpoints)));
            zpoints = reshape(zpoints, 1, prod(size(zpoints)));
            scatter(ax, rpoints,zpoints,point_size,point_color,'.');
        end
    end

    if hasbodies
        for j = 1:Nbody
            r = body_radii(j);
            br = body_coords(j,1);
            bz = body_coords(j,2);

            theta = 0:0.01:pi;
            circ_x = r*sin(theta)+br; circ_x = [circ_x,circ_x(1)];
            circ_y = r*cos(theta)+bz; circ_y = [circ_y,circ_y(1)];
            patch(circ_x,circ_y,body_colors{j},'EdgeColor', 'none', 'Parent', ax);
        end
    end

    xlabel('r');
    ylabel('z');
    axis equal
    xlim([rmin rmax])
    ylim([zmin zmax])

elseif (dim == 3)
    % 3D
    xmin = Inf;
    xmax = -Inf;
    ymin = Inf;
    ymax = -Inf;
    zmin = Inf;
    zmax = -Inf;
    for j = 1:Nmesh
        this_xss = meshes{j}{1};
        this_yss = meshes{j}{2};
        this_zss = meshes{j}{3};
        this_xbound = [this_xss(1) this_xss(end)];
        this_ybound = [this_yss(1) this_yss(end)];
        this_zbound = [this_zss(1) this_zss(end)];
        xmin = min(xmin,min(this_xss));
        xmax = max(xmax,max(this_xss));
        ymin = min(ymin,min(this_yss));
        ymax = max(ymax,max(this_yss));
        zmin = min(zmin,min(this_zss));
        zmax = max(zmax,max(this_zss));

        for i1 = 1:2; for i2 = 1:2;
            plot3(ax, [this_xbound(1) this_xbound(2)], [this_ybound(i1) this_ybound(i1)], [this_zbound(i2) this_zbound(i2)], edge_color);
            hold(ax, 'on');
            plot3(ax, [this_xbound(i1) this_xbound(i1)], [this_ybound(1) this_ybound(2)], [this_zbound(i2) this_zbound(i2)], edge_color);
            plot3(ax, [this_xbound(i1) this_xbound(i1)], [this_ybound(i2) this_ybound(i2)], [this_zbound(1) this_zbound(2)], edge_color);
        end; end
        if ~strcmpi(point_color,'none')
            [xpoints, ypoints, zpoints] = ndgrid(meshes{j}{1},meshes{j}{2},meshes{j}{3});
            xpoints = reshape(xpoints, 1, prod(size(xpoints)));
            ypoints = reshape(ypoints, 1, prod(size(ypoints)));
            zpoints = reshape(zpoints, 1, prod(size(zpoints)));
            scatter3(ax, xpoints,ypoints,zpoints,point_size,point_color,'.');
        end
    end

    if hasbodies
        [px,py,pz] = sphere;
        for j = 1:Nbody
            r = body_radii(j);
            bx = body_coords(j,1);
            by = body_coords(j,2);
            bz = body_coords(j,3);
            surf(ax, px*r+bx,py*r+by,pz*r+bz, 'FaceColor', body_colors{j}, 'EdgeColor', 'none');
        end
    end

    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    xlim([xmin xmax])
    ylim([ymin ymax])
    zlim([zmin zmax])
end

if ~holding
    hold(ax, 'off')
end

