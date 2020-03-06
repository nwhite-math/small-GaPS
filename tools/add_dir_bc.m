% ADD_DIR_BC Apply Dirichlet boundary condition to a discretized differential operator matrix.
%
% [Areduced, K, ninds, addbfactor, addxfactor] = ADD_DIR_BC(inds, A, bcval)
% [Areduced, K, ninds] = ADD_DIR_BC(inds, A)
%
% inds: the indices where the BC should be implemented
% A: the operator to be reduced (its order will go down by |inds|)
% bcval: the BC value (such that f = bcval on the boundary)
%
% Areduced is the reduced operator.
% ninds is the set of indices not in inds.
%
% K is the restoration matrix. Suppose we reduced A to Ared.
% Then we could find some eigenvalue Vred of Ared, but Vred would have size |ninds|.
% The full restored V is: V = K*Vred.
%
% By default, bcval is assumed to be 0, i.e. the boundary conditions are homogeneous.
% If bcval is nonzero, then addbfactor is returned. Then, the problem (A f = b) with
% boundary conditions becomes (Areduced * f(ninds) = b(ninds) + addbfactor).
%
% Note that to solve eigenvalue problems, it is necessary to use a homogeneous BC.
%
% Example 1: Heat equation with homogeneous Dirichlet conditions.
%    >> N = 11;
%    >> x = linspace(0,1,N);
%    >> [pd,pd2,pd3,pd4] = deriv_findiff(x,2);
%    >>
%    >> A = pd2; % the heat equation differential operator
%    >>
%    >> inds = [1 N];         % Indices where we apply BC
%    >> [Ared, K, ninds] = ADD_DIR_BC(inds, A);
%    >>
%    >> [Vred,ered] = eigs(Ared,1,'lr'); % evec of largest real eval
%    >> V = K*Vred; % restored full vector
%    >> plot(x,V);
%
%
% Example 2: Poisson equation with inhomogeneous Dirichlet conditions.
%            Solve d^2 f/dx^2 = -1, f(0) = 1, f(1) = 0.
%
%    >> N = 11;
%    >> x = linspace(0,1,N);
%    >> [pd,pd2,pd3,pd4] = deriv_findiff(x,2);
%    >>
%    >> A = pd2; % the heat equation differential operator
%    >>
%    >> m = [1 0];              % Left BC = 1, right = 0
%    >> b = repmat(-1.0, N, 1); % Poisson equation value
%    >>
%    >>
%    >> inds = [1 N];         % Indices where we apply BC
%    >> [Ared, K, ninds, addbfac, addxfac] = ADD_DIR_BC(inds, A, m);
%    >>
%    >> fred = Ared\(b(ninds) + addbfac);
%    >> f = K*fred + addxfac;
%    >>
%    >> plot(x,f);

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-02-15
% Updated: 2019-02-15

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


function [Ared,K,ninds, varargout] = add_dir_bc(inds, A, bcval);

hasbcval = true;
if (nargin <= 2);
    hasbcval = false;
end;

N = length(A);

% Make sure the input is good
assert(min(inds) >= 1);
assert(max(inds) <= N);
assert(N > length(inds));
assert(length(size(A)) == 2);
assert(size(A,1) == size(A,2));
if (hasbcval);
    assert(min(size(bcval)) == 1);
    if(length(bcval) == N)
        bcval = bcval(inds);
    else
        assert(length(bcval) == length(inds));
    end
    if size(bcval,1) < size(bcval,2)
        bcval = bcval';
    end
end;

% set up "not inds" array
ninds = [1:N];
ninds(inds) = 0;
ninds = find(ninds);

% K acting only on inds
K = spalloc(N, length(ninds), length(ninds));
K(ninds,:) = speye(length(ninds));

Ared = A(ninds,ninds);

varargout = {};
if (hasbcval);
    addbfactor = -A(ninds,inds)*bcval;
    varargout{1} = addbfactor;
    addxfactor = repmat(0, N, 1);
    addxfactor(inds) = bcval;
    varargout{2} = addxfactor;
end;
