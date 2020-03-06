% DERIV_FINDIFF Return differentiation matrices using finite difference
%
% [pd,pd2,pd3,pd4] = DERIV_FINDIFF(x,order)
%
% Returns first derivative matrix (pd), second derivative matrix (pd2), third
% (pd3), and fourth (pd4).
% Note: If there are insufficient points or the order is too low for higher
% derivatives, they will be returned as zero matrices.
%
% x is the domain vector, not necessarily equispaced.
% order is the differentiation order, which must be an even integer.
%
% Example:
%    >> x = linspace(0,4,5)
%    [ 0 1 2 3 4 ]
%    >> DERIV_FINDIFF(x,2)
%    -1.5   2    -0.5   0     0
%    -0.5   0     0.5   0     0
%     0    -0.5   0     0.5   0
%     0     0    -0.5   0     0.5
%     0     0     0.5  -2     1.5
%

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2016-05-03
% Updated: 2019-05-31

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


function [pd,pd2,pd3,pd4] = deriv_findiff(x, order);

% parse arguments

assert(min(size(x))==1, 'x must be a vector, not a matrix.');
if (size(x,2)==1); x = transpose(x); end; % make x a row vector
N = length(x);

% require even order
assert(round(order/2)==order/2, 'Derivative order must be even.');

% 2nd order needs at least 3 points, etc.
assert(N > order, 'The number of input points must be greater than the derivative order.');

pd = zeros(N);
pd2 = zeros(N);
pd3 = zeros(N);
pd4 = zeros(N);

% Construct derivative matrices by:
% 1) Constructing Taylor expansion at each point in the mesh
% 2) Inverting it
% 3) Combining the results from every point

% ix = x index
% Our convention will be regular letters (e.g. i or j) for _relative_ indices in loops,
% and letters with x appended (e.g. ix, jx) for the corresponding mesh indices.
for (ix=1:N);

    % A will be the Taylor expansion matrix at point ix
    % Let D be the column vector [f(ix), f'(ix), f''(ix), ...], the derivatives of some f at ix.
    % Then, A.D=F should give a column vector F, which is the Taylor expansion approximation of
    % [f(1), f(2), ..., f(N)] around the point f(ix).
    %
    % For example, say x = [1, 3, 11], and we're considering the first point (ix=1).
    % Then, f(1) = f(1),
    %       f(3) ~ f(1) + 2*f'(1) + 0.5*(2)^2*f''(1)
    %       f(11) ~ f(1) + 10*f'(1) + 0.5*(10)^2*f''(1)
    % So,
    % A = [[ 1, 0, 0],
    %      [ 1, 2, 2],
    %      [ 1, 10, 50]]
    %
    % Then the inverse of this matrix, Am=inv(A), can be multiplied by [f(1),f(3),f(11)]
    % and will return [f(1), f'(1), f''(1)].
    %
    % Combining the second row of the Am matrices at each point ix will give pd, etc.

    A = zeros(order+1);

    % starti is the offset index for the leftmost point used to compute the derivative.
    % For example, say order = 2 and our mesh has 10 points.
    % Then we'll compute the derivative at ix=1 from indices 1,2,3; hence starti = 0.
    % At some point inside the mesh, say ix=5, we'll compute the derivative from indices 4,5,6;
    % hence starti = -1.
    % Finally, at the right end, ix=10 and we compute the derivative from 8,9,10, so starti = -2.
    starti = -order/2;
    starti = max((1 - ix),starti); % adjust in case ix is too far left
    starti = min(N-(order + ix),starti); % adjust for ix too far right

    % Construct the Taylor expansion A matrix, as described above.
    facs = factorial([0:order]);
    for j=1:(order+1);
        jx = j + starti - 1;
        h_here = x(ix+jx) - x(ix);
        jxhs = ones(1,order+1)*h_here;
        A(j,:) = ((jxhs).^[0:order])./facs;
    end;

    % As mentioned above, the second row of Am is used to get f', the third row will give f'', etc.
    Am = inv(A);
    pdline = Am(2,:);
    pd2line = Am(3,:);
    if (order > 2);
        pd3line = Am(4,:);
    else;
        pd3line = 0*pdline;
    end;
    if (order > 3);
        pd4line = Am(5,:);
    else;
        pd4line = 0*pdline;
    end;
    pdline = [zeros(1,ix+starti-1), pdline, zeros(1,N-(ix+starti-1)-(order+1))];
    pd2line = [zeros(1,ix+starti-1), pd2line, zeros(1,N-(ix+starti-1)-(order+1))];
    pd3line = [zeros(1,ix+starti-1), pd3line, zeros(1,N-(ix+starti-1)-(order+1))];
    pd4line = [zeros(1,ix+starti-1), pd4line, zeros(1,N-(ix+starti-1)-(order+1))];

    pd(ix,:) = pdline;
    pd2(ix,:) = pd2line;
    pd3(ix,:) = pd3line;
    pd4(ix,:) = pd4line;
end;
