% CONTOURF_INTERP_3D Construct contour plot of a 2D slice with interpolated data.
%
% [c,h, coord1, coord2, data] = CONTOURF_INTERP_3D(meshes, idgrids, nodedata,...
%   data, xslice, yslice, zslice)
% [c,h, coord1, coord2, data] = CONTOURF_INTERP_3D(ax, meshes, idgrids, nodedata,...
%   data, xslice, yslice, zslice)
% [c,h, coord1, coord2, data] = CONTOURF_INTERP_3D(ax, meshes, idgrids, nodedata,...
%   data, xslice, yslice, zslice,...
%  'XRange',              [xmin xmax],   ...
%  'YRange',              [ymin ymax],   ...
%  'ZRange',              [zmin zmax],   ...
%  'XResolution',         xres,          ...
%  'YResolution',         yres,          ...
%  'ZResolution',         zres,          ...
%  'Levels',              num_levels,    ...
%  'LevelList',           [level list],  ...
%  'LineColor',           line_color,    ...
%  'Fill',                true/false,    ...
%  'AxesEqual',           true/false,    ...
%  'XOffset',             x_offset,      ...
%  'YOffset',             y_offset,      ...
%  'ZOffset',             z_offset,      ...
%  'DistanceMultiplier',  dist_mult,     ...
%  )
%
%
% meshes, idgrids, and nodedata must be the corresponding variables generated by
% multimesh_3d.
% data should be a vector of real values, arranged to correspond to the nodes of
% the multimesh.
%
% CONTOURF_INTERP_3D interpolates the data in the multimesh and displays a
% contour plot of the result. The output is the same as that of the standard
% CONTOURF, c being the contour data and h a handle to the contour object.
%
% Currently, two of xslice, yslice, and zslice must be [], and the other must
% be a single number. For example, if xslice = [], yslice = [], and zslice = 1,
% then the output will be a slice in the x-y plane at z=1.
%
% Several optional arguments are allowed.
% Plot settings:
%  'XRange':             Min and max x values setting the plot limits
%  'YRange':             Min and max y values setting the plot limits
%  'ZRange':             Min and max z values setting the plot limits
%  'XResolution':        Number of points in x for the interpolated grid.
%                        Default: 50
%  'YResolution':        Number of points in y for the interpolated grid.
%                        Default: 50
%  'ZResolution':        Number of points in z for the interpolated grid.
%                        Default: 50
%  'Levels':             Number of contour levels to draw
%  'LevelList':          Specific list of values at which to draw contours
%  'LineColor':          Color of contour lines
%  'Fill':               true => plot with CONTOURF. false => plot with CONTOUR
%  'AxesEqual':          true => x and y will have the same scale
%
%  A few additional options are provided for modifying the x and y axes, since
%  it is inconvenient to do so from within nodedata.
%  'XOffset':            Constant to add to x
%  'YOffset':            Constant to add to y
%  'ZOffset':            Constant to add to z
%  'DistanceMultiplier': Constant by which to multiply x, y, and z
%
% Note that only the arguments relevant to the plotted plane are used. E.g.,
% if the plot will be in the x-y plen, then 'ZRange', 'ZResolution', etc. are
% ignored.

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-08-30
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


function [c, h, coord1_out, coord2_out, data_interp] = contourf_slice_3d(varargin)

hasax = isa(varargin{1},'matlab.graphics.axis.Axes');



parser = inputParser;
if (hasax)
    parser.addRequired('axes');
end
parser.addRequired('meshes', @(x) iscell(x));
parser.addRequired('idgrids', @(x) iscell(x));
parser.addRequired('nodedata', @(x) isstruct(x));
parser.addRequired('data', @(x) isvector(x));
parser.addRequired('xslice', @(x) (isvector(x) || isempty(x)));
parser.addRequired('yslice', @(x) (isvector(x) || isempty(x)));
parser.addRequired('zslice', @(x) (isvector(x) || isempty(x)));
parser.addParameter('XRange', [], @(x) (isvector(x) && isnumeric(x) && length(x)==2));
parser.addParameter('YRange', [], @(x) (isvector(x) && isnumeric(x) && length(x)==2));
parser.addParameter('ZRange', [], @(x) (isvector(x) && isnumeric(x) && length(x)==2));
parser.addParameter('XResolution', 50, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('YResolution', 50, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('ZResolution', 50, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('LevelList', [], @(x) (isvector(x) && isnumeric(x) && length(x) > 1));
parser.addParameter('Levels', 6, @(x) (isscalar(x) && isnumeric(x)));
parser.addParameter('LineColor', 'k', @(x) ischar(x));
parser.addParameter('Fill', true, @(x) isscalar(x));
parser.addParameter('DistanceMultiplier', 1, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('XOffset', 0, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('YOffset', 0, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('ZOffset', 0, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('AxesEqual', false, @(x) isscalar(x));
parse(parser, varargin{:});
meshes = parser.Results.meshes;
idgrids = parser.Results.idgrids;
nodedata = parser.Results.nodedata;
data = parser.Results.data;
xslice = parser.Results.xslice;
yslice = parser.Results.yslice;
zslice = parser.Results.zslice;
xrange = parser.Results.XRange;
yrange = parser.Results.YRange;
zrange = parser.Results.ZRange;
xres = parser.Results.XResolution;
yres = parser.Results.YResolution;
zres = parser.Results.ZResolution;
levellist = parser.Results.LevelList;
levelsN = parser.Results.Levels;
linecolor = parser.Results.LineColor;
should_fill = parser.Results.Fill;
dist_conv = parser.Results.DistanceMultiplier;
xoffset = parser.Results.XOffset;
yoffset = parser.Results.YOffset;
zoffset = parser.Results.ZOffset;
equal_axes = parser.Results.AxesEqual;

if (length(xslice) == 1) && (length(yslice) == 0) && (length(zslice) == 0)
    dir = 1;
    assert( (yres >= 2) && (zres >= 2) && ...
            (yres <= 1000) && (zres <= 1000),...
            'YResolution and ZResolution must be between 2 and 1000');
    yres = ceil(yres);
    zres = ceil(zres);

    i1 = 2;
    i2 = 3;
    res1 = yres;
    res2 = zres;
    offset1 = yoffset;
    offset2 = zoffset;
    range1 = yrange;
    range2 = zrange;
    nd1 = nodedata.y;
    nd2 = nodedata.z;
    nd3 = nodedata.x;
    slice = xslice;
elseif (length(xslice) == 0) && (length(yslice) == 1) && (length(zslice) == 0)
    dir = 2;
    assert( (xres >= 2) && (zres >= 2) && ...
            (xres <= 1000) && (zres <= 1000),...
            'XResolution and ZResolution must be between 2 and 1000');
    zres = ceil(zres);
    xres = ceil(yres);

    i1 = 3;
    i2 = 1;
    res1 = zres;
    res2 = xres;
    offset1 = zoffset;
    offset2 = xoffset;
    range1 = zrange;
    range2 = xrange;
    nd1 = nodedata.z;
    nd2 = nodedata.x;
    nd3 = nodedata.y;
    slice = yslice;
elseif (length(xslice) == 0) && (length(yslice) == 0) && (length(zslice) == 1)
    dir = 3;
    assert( (yres >= 2) && (xres >= 2) && ...
            (yres <= 1000) && (xres <= 1000),...
            'XResolution and YResolution must be between 2 and 1000');
    xres = ceil(xres);
    yres = ceil(yres);

    i1 = 1;
    i2 = 2;
    res1 = xres;
    res2 = yres;
    offset1 = xoffset;
    offset2 = yoffset;
    range1 = xrange;
    range2 = yrange;
    nd1 = nodedata.x;
    nd2 = nodedata.y;
    nd3 = nodedata.z;
    slice = zslice;
else
    assert(false,'Two of xslice, yslice, and zslice must be [], and the other must be a scalar.');
end

if (hasax)
    ax = parser.Results.axes;
else
    ax = gca;
end

holding = ishold(ax);

meshN = length(meshes);
assert(length(idgrids) == meshN);

min1 = Inf;
max1 = -Inf;
min2 = Inf;
max2 = -Inf;
min3 = Inf;
max3 = -Inf;
for j = 1:meshN
    cs1 = meshes{j}{i1};
    cs2 = meshes{j}{i2};
    cs3 = meshes{j}{dir};
    min1 = min(min1, min(cs1));
    max1 = max(max1, max(cs1));
    min2 = min(min2, min(cs2));
    max2 = max(max2, max(cs2));
    min3 = min(min3, min(cs3));
    max3 = max(max3, max(cs3));
end
if length(range1) == 0
    range1 = [min1 max1];
end
if length(range2) == 0
    range2 = [min2 max2];
end

assert( (range1(2) > range1(1)) && (range2(2) > range2(1)),...
    'Relevant input ranges (of XRange, YRange, and ZRange) must be nonzero.');
[~, min3] = nearest_index(nd3, slice, 'below');
[~, max3] = nearest_index(nd3, slice, 'above');
range3 = [min3 max3];
assert((min3 <= slice) && (max3 >= slice), 'slice is outside the data range.');


% Setup contour levels
% TODO: should this be done _after_ interpolation?
if length(levellist) == 0
    checkdata = data( (~isnan(data)) & (nd1' >= range1(1)) & (nd1' <= range1(2)) & (nd2' >= range2(1)) & (nd2' <= range2(2)) & (nd3' >= range3(1)) & (nd3' <= range3(2)) );
    cmin = min(checkdata);
    cmax = max(checkdata);

    if ((length(cmin) ==1) && (length(cmax)==1) && ~isnan(cmin) && ~isnan(cmax))
        levellist = [0:1/(levelsN+1):1]*(cmax - cmin) + cmin;
    else
        levellist = [];
    end
else
    cmin = min(levellist);
    cmax = max(levellist);
end

% Interpolate data
rs_interp1 = ([0:res1-1]/(res1-1)) * (range1(2) - range1(1)) + range1(1);
rs_interp2 = ([0:res2-1]/(res2-1)) * (range2(2) - range2(1)) + range2(1);
[rg_interp1, rg_interp2, rg_interp3] = meshgrid(rs_interp1, rs_interp2, slice);
data_interp = griddata(nd1, nd2, nd3, data, rg_interp1, rg_interp2, rg_interp3, 'natural');

coord1_out = (rg_interp1+offset1)*dist_conv;
coord2_out = (rg_interp2+offset2)*dist_conv;

% Contour interpolated data
if should_fill
    [c, h]  = contourf(ax, (rg_interp1+offset1)*dist_conv, (rg_interp2+offset2)*dist_conv, data_interp, 'LevelList', levellist, 'LineColor', linecolor);
    if ((length(cmin) ==1) && (length(cmax)==1) && ~isnan(cmin) && ~isnan(cmax))
        caxis(ax, [cmin cmax]);
    end
else
    [c, h]  = contour(ax, (rg_interp1+offset1)*dist_conv, (rg_interp2+offset2)*dist_conv, data_interp, 'LevelList', levellist, 'LineColor', linecolor);
end
if (equal_axes)
    axis(ax,'equal');
end
xlim(ax, (range1+offset1)*dist_conv);
ylim(ax, (range2+offset2)*dist_conv);

if (holding)
    hold(ax, 'on')
else
    hold(ax, 'off')
end

end


