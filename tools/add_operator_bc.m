% ADD_OPERATOR_BC Apply boundary condition(s) to a discretized differential operator matrix.
%
% [Areduced, K, ninds, addbfactor, addxfactor] = ADD_OPERATOR_BC(Q, inds, A, bcval)
% [Areduced, K, ninds] = ADD_OPERATOR_BC(Q, inds, A)
% [NaN, K, ninds] = ADD_OPERATOR_BC(Q, inds)
%
% Q: a BC operator (e.g. flux operator, or just [1 0 0 0; ...] for Dirichlet)
% inds: the indices where the BC should be implemented
% A: the operator to be reduced (its order will go down by |inds|)
% bcval: the BC value (such that Q*f = bcval on the boundary)
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
% Example 1: Heat equation with homogeneous Dirichlet condition on the left,
%            and Neumann on the right.
%    >> N = 11;
%    >> x = linspace(0,1,N);
%    >> [pd,pd2,pd3,pd4] = deriv_findiff(x,2);
%    >>
%    >> A = pd2; % the heat equation differential operator
%    >>
%    >> Q = zeros(N);
%    >> Q(1) = 1;             % Dirichlet BC on the first row.
%    >> Q(end,:) = pd(end,:); % Neumann BC on the last row.
%    >>
%    >> inds = [1 N];         % Indices where we apply BC
%    >> [Ared, K, ninds] = ADD_OPERATOR_BC(Q, inds, A);
%    >>
%    >> [Vred,ered] = eigs(Ared,1,'lr'); % evec of largest real eval
%    >> V = K*Vred; % restored full vector
%    >> plot(x,V);
%
%
% Example 2: Poisson equation with inhomogeneous Dirichlet condition on the left,
%            and Neumann on the right.
%            Solve d^2 f/dx^2 = -1, f(0) = 1, f'(1) = 0.
%
%    >> N = 11;
%    >> x = linspace(0,1,N);
%    >> [pd,pd2,pd3,pd4] = deriv_findiff(x,2);
%    >>
%    >> A = pd2; % the heat equation differential operator
%    >>
%    >> Q = zeros(N);
%    >> Q(1) = 1;             % Dirichlet BC on the first row.
%    >> Q(end,:) = pd(end,:); % Neumann BC on the last row.
%    >> m = [1 0];            % Dirichlet BC = 1, Neumann = 0
%    >> b = repmat(-1.0, N, 1); % Poisson equation value
%    >>
%    >>
%    >> inds = [1 N];         % Indices where we apply BC
%    >> [Ared, K, ninds, addbfac, addxfac] = ADD_OPERATOR_BC(Q, inds, A, m);
%    >>
%    >> fred = Ared\(b(ninds) + addbfac);
%    >> f = K*fred + addxfac;
%    >>
%    >> plot(x,f);

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2016-06-23
% Updated: 2019-01-25

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


function [Ared,K,ninds, varargout] = add_operator_bc(Q, inds, A, bcval);

hasA = true;
hasbcval = true;
if (nargin <= 3);
    hasbcval = false;
end;
if (nargin <= 2);
    hasA = false;
end;

N = size(Q,2);

% Make sure the input is good
assert(min(inds) >= 1);
assert(max(inds) <= N);
assert(N > length(inds));
assert(size(Q,1) >= length(inds));
assert(size(Q,1) <= N);
if (hasA);
    assert(size(A,1) == size(A,2));
    assert(size(A,1) == N);
end;
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

% you can give Q as a full square matrix;
% the irrelevant rows will be thrown out here
if (size(Q,1)==N);
    Q = Q(inds,:);
end;

% set up "not inds" array
ninds = [1:N];
ninds(inds) = 0;
ninds = find(ninds);

% K acting only on inds
Qinv = inv(Q(:,inds));
Kb = -Qinv*Q(:,ninds);

K = spalloc(N, length(ninds), length(ninds)+length(inds)*N);
K(ninds,:) = speye(length(ninds));
K(inds,:) = Kb;

if (hasA);
    Ared = A(ninds,ninds)+A(ninds,inds)*Kb;
else;
    Ared = NaN;
end;

varargout = {};
if (hasbcval);
    addbfactor = -A(ninds,inds)*Qinv*bcval;
    varargout{1} = addbfactor;
    addxfactor = repmat(0, N, 1);
    addxfactor(inds) = Qinv*bcval;
    varargout{2} = addxfactor;
end;
