% CONTOURF_INTERP Construct a contour plot with interpolated data.
%
% [c, h, x, y, data] = CONTOURF_INTERP(meshes, idgrids, nodedata, data)
% [c, h, x, y, data] = CONTOURF_INTERP(ax, meshes, idgrids, nodedata, data)
% [c, h, x, y, data] = CONTOURF_INTERP(ax, meshes, idgrids, nodedata, data,...
%       'XRange',              [xmin xmax],   ...
%       'YRange',              [ymin ymax],   ...
%       'XResolution',         xres,          ...
%       'YResolution',         yres,          ...
%       'Levels',              num_levels,    ...
%       'LevelList',           [level list],  ...
%       'LineColor',           line_color,    ...
%       'Fill',                true/false,    ...
%       'AxesEqual',           true/false,    ...
%       'MeshBoundaryColor',   'none',        ...
%       'XOffset',             x_offset,      ...
%       'YOffset',             y_offset,      ...
%       'DistanceMultiplier',  dist_mult,     ...
%       )
%
%
% meshes, idgrids, and nodedata must be the corresponding variables generated by
% multimesh_2d.
% data should be a vector of real values, arranged to correspond to the nodes of
% the multimesh.
%
% CONTOURF_INTERP interpolates the data in the multimesh and displays a contour
% plot of the result. The output is the same as that of the standard CONTOURF, c
% being the contour data and h a handle to the contour object.
%
% Several optional arguments are allowed.
% Plot settings:
%  'XRange':             Min and max x values setting the plot limits
%  'YRange':             Min and max y values setting the plot limits
%  'XResolution':        Number of points in x for the interpolated grid.
%                        Default: 50
%  'YResolution':        Number of points in y for the interpolated grid.
%                        Default: 50
%  'Levels':             Number of contour levels to draw
%  'LevelList':          Specific list of values at which to draw contours
%  'LineColor':          Color of contour lines
%  'Fill':               true => plot with CONTOURF. false => plot with CONTOUR
%  'AxesEqual':          true => x and y will have the same scale
%  'MeshBoundaryColor':  Color of mesh boundaries (default: 'none')
%
%  A few additional options are provided for modifying the x and y axes, since
%  it is inconvenient to do so from within nodedata.
%  'XOffset':            Constant to add to x
%  'YOffset':            Constant to add to y
%  'DistanceMultiplier': Constant by which to multiply x and y
%

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-05-11
% Updated: 2020-01-01

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


function [c,h, xg_out, yg_out, data_interp] = contourf_interp(varargin)

hasax = isa(varargin{1},'matlab.graphics.axis.Axes');



parser = inputParser;
if (hasax)
    parser.addRequired('axes');
end
parser.addRequired('meshes', @(x) iscell(x));
parser.addRequired('idgrids', @(x) iscell(x));
parser.addRequired('nodedata', @(x) isstruct(x));
parser.addRequired('data', @(x) (isvector(x) && min(size(x))==1));
parser.addParameter('XRange', [], @(x) (isvector(x) && isnumeric(x) && length(x)==2 && min(size(x))==1));
parser.addParameter('YRange', [], @(x) (isvector(x) && isnumeric(x) && length(x)==2 && min(size(x))==1));
parser.addParameter('XResolution', 50, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('YResolution', 50, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('LevelList', [], @(x) (isvector(x) && isnumeric(x) && min(size(x))==1 && length(x) > 1));
parser.addParameter('Levels', 6, @(x) (isscalar(x) && isnumeric(x)));
parser.addParameter('LineColor', 'k', @(x) ischar(x));
parser.addParameter('Fill', true, @(x) isscalar(x));
parser.addParameter('DistanceMultiplier', 1, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('XOffset', 0, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('YOffset', 0, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('AxesEqual', false, @(x) isscalar(x));
parser.addParameter('MeshBoundaryColor', 'none', @(x) ischar(x));
parse(parser, varargin{:});
meshes = parser.Results.meshes;
idgrids = parser.Results.idgrids;
nodedata = parser.Results.nodedata;
data = parser.Results.data;
xrange = parser.Results.XRange;
yrange = parser.Results.YRange;
xres = parser.Results.XResolution;
yres = parser.Results.YResolution;
levellist = parser.Results.LevelList;
levelsN = parser.Results.Levels;
linecolor = parser.Results.LineColor;
should_fill = parser.Results.Fill;
dist_conv = parser.Results.DistanceMultiplier;
xoffset = parser.Results.XOffset;
yoffset = parser.Results.YOffset;
equal_axes = parser.Results.AxesEqual;
boundarycolor = parser.Results.MeshBoundaryColor;

assert( (xres >= 2) && (yres >= 2) && (xres <= 1000) && (yres <= 1000) , 'XResolution and YResolution must be between 2 and 1000');
xres = ceil(xres);
yres = ceil(yres);

if (hasax)
    ax = parser.Results.axes;
else
    ax = gca;
end

holding = ishold(ax);

meshN = length(meshes);
assert(length(idgrids) == meshN);

xmin = Inf;
xmax = -Inf;
ymin = Inf;
ymax = -Inf;
for j = 1:meshN
    xs = meshes{j}{1};
    ys = meshes{j}{2};
    xmin = min(xmin, min(xs));
    xmax = max(xmax, max(xs));
    ymin = min(ymin, min(ys));
    ymax = max(ymax, max(ys));
end
if length(xrange) == 0
    xrange = [xmin xmax];
end
if length(yrange) == 0
    yrange = [ymin ymax];
end
assert( (xrange(2) > xrange(1)) && (yrange(2) > yrange(1)) , 'XRange and YRange must be nonzero.');

assert(length(data) == length(nodedata.x));
assert(isvector(data));
if size(data, 1)==1
    data = data';
end

% Setup contour levels
% TODO: should this be done _after_ interpolation?
if length(levellist) == 0
    checkdata = data( (~isnan(data)) & (nodedata.x' >= xrange(1)) & (nodedata.x' <= xrange(2)) & (nodedata.y' >= yrange(1)) & (nodedata.y' <= yrange(2)));
    cmin = min(checkdata(~isinf(checkdata)));
    cmax = max(checkdata(~isinf(checkdata)));

    if ~isempty(cmin) && ~isempty(cmax)
        levellist = [0:1/(levelsN+1):1]*(cmax - cmin) + cmin;
    else
        levellist = [];
    end
else
    cmin = min(levellist);
    cmax = max(levellist);
end

% Interpolate data
xs_interp = ([0:xres-1]/(xres-1)) * (xrange(2) - xrange(1)) + xrange(1);
ys_interp = ([0:yres-1]/(yres-1)) * (yrange(2) - yrange(1)) + yrange(1);
[xg_interp, yg_interp] = ndgrid(xs_interp, ys_interp);
data_interp = griddata(nodedata.x, nodedata.y, data, xg_interp, yg_interp, 'cubic');

xg_out = (xg_interp+xoffset)*dist_conv;
yg_out = (yg_interp+yoffset)*dist_conv;

% Contour interpolated data
if should_fill
    [c, h]  = contourf(ax, xg_out, yg_out, data_interp, 'LevelList', levellist, 'LineColor', linecolor);
    if ~isempty(cmin) && ~isempty(cmax)
        caxis(ax, [cmin cmax]);
    end
else
    [c, h]  = contour(ax, xg_out, yg_out, data_interp, 'LevelList', levellist, 'LineColor', linecolor);
end
if (equal_axes)
    axis(ax,'equal');
end
xlim(ax, (xrange+xoffset)*dist_conv);
ylim(ax, (yrange+yoffset)*dist_conv);

% Draw mesh boundaries
if ~strcmpi(boundarycolor, 'none')
    for j = 1:meshN
        this_xmin = min(meshes{j}{1});
        this_xmax = max(meshes{j}{1});
        this_ymin = min(meshes{j}{2});
        this_ymax = max(meshes{j}{2});
        hold(ax, 'on')
        rectangle('Position', [this_xmin, this_ymin, (this_xmax - this_xmin), (this_ymax - this_ymin)], 'EdgeColor', boundarycolor);
    end
end

if (holding)
    hold(ax, 'on')
else
    hold(ax, 'off')
end

end
