% SORTTREE Sort a tree graph so that parents always come before their children.
%
% i = SORTTREE(t)
%
% i is the list of indices in sorted order
% t must be a vector with each element indicating the index of its parent, and
% 0 being the root.
%
%
% For example, t = [2 0 1 7 2 1 0] represents the graph:
%
%     2      7
%    / \     |
%   1   5    4
%  / \
% 3   6
%
% In this case, the usage and result would be:
%    >> i = SORTTREE([2 0 1 7 2 1 0])
%    i = [ 2 7 1 5 4 3 6 ]

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-03-18
% Updated: 2019-03-18

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


function i = sorttree(t)

validateattributes(t, {'numeric'}, {'2d'});
assert(min(size(t)) == 1)
assert(sum(t==0)>=1);
N = length(t);
assert(all(t == round(t)));
assert(max(t) <= N);
assert(min(t) == 0);

i = [];
level1_i = [0];
for j=1:N;
    level2_i = [];
    for k = level1_i
        level2_i = [level2_i find(t == k)];
    end
    if length(level2_i) == 0
        break
    else
        i = [i level2_i];
        level1_i = level2_i;
    end
end
assert(length(i) == length(t), 'Graph must be a tree (no cycles)');
