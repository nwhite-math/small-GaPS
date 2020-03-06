% VAINSHTEINONELAP Return the Laplacian of the Vainshtein field of a single
% spherical body.
%
% phi_nd = VAINSHTEINONELAP(k, rho_nd, rs_nd, r_nd)
%
% VAINSHTEINONELAP returns the nondimensionalized Laplacian of the Vainshtein
% field of a spherical body with nondimensional density rho_nd and non-
% dimensional radius rs_nd, where the equation has constant k.
%
% The original equation is:
%  ,--------------------------------------------------------------------------,
%  |  0 = 3 lap phi + r_c^2[  (lap phi)^2 - d_ij phi d^ij phi ] - 8 pi G rho  |
%  '--------------------------------------------------------------------------'
%
% which is nondimensionalized by a given length scale d and reference density rho_c,
% yielding:
%  ,-------------------------------------------,
%  |  phi_c  = d^2 sqrt(8 pi G rho_c ) / r_c   |
%  |  k      = 3 / [ sqrt(8 pi G rho_c) r_c ]  |
%  |  r_nd   = r / d                           |
%  |  rs_nd  = rs / d                          |
%  |  rho_nd = rho / rho_c                     |
%  '-------------------------------------------'
%                       ||
%                       \/
%  ,-----------------------------------------------------------------,
%  |  0 = k lap phi + [  (lap phi)^2 - d_ij phi d^ij phi ] - rho_nd  |
%  '-----------------------------------------------------------------'
%
% If rho = 0, then 0 is returned.
%
% If k = 0, then the nonlinear-only analytic solution is returned.

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-04-11
% Updated: 2019-04-11

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


function lapphi = vainshteineonelap(k, rho, rs, r)

    usage_str = 'vainshteineonelap(k, rho_nd, rs_nd, r_nd)';

    % Note: COMSOL makes calls with giant vectors for no reason, so we can't do these assertions....
    %assert(max(size(k)) == 1, sprintf('%s\n%s',usage_str, 'k must be a scalar.'));
    %assert(max(size(rho)) == 1, sprintf('%s\n%s',usage_str, 'rho_nd must be a scalar.'));
    %assert(max(size(rs)) == 1, sprintf('%s\n%s',usage_str, 'rs_nd must be a scalar.'));

    k = k(1);
    rho = rho(1);
    rs = rs(1);

    % If rho == 0, return 0
    if (rho == 0)
        lapphi = 0*r;
        return
    end

    % If k == 0, return the nonlinear-only analytic solution
    if (k == 0)
        %lapphi = sqrt(rho*3/2)*( (1).*(r<=rs) +(2*(rs./r).^(3/2)).*(r>rs) );
        lapphi = sqrt(rho*3/2)*( (1) + 0*r );
        lapphi(r>rs) = sqrt(rho*3/2)*( (2*(rs./r(r>rs)).^(3/2)) );
        return
    end

    % rho ~= 0 and k ~= 0, so return the full solution
    %lapphi = 0.75*k*(     -1 + ( (sqrt(9*k^2+24*rho)/(3*k)).*(r<=rs) + ( (3*k^2*r.^3 + 4*rs^3*rho)./(k*r.^(3/2).*sqrt(9*k^2*r.^3 + 24*rs^3*rho)) ).*(r>rs) )    );
    lapphi = 0.75*k*(     -1 + ( (sqrt(9*k^2+24*rho)/(3*k)) + 0*r )    );
    lapphi(r>rs) = 0.75*k*(     -1 + (  + ( (3*k^2*r(r>rs).^3 + 4*rs^3*rho)./(k*r(r>rs).^(3/2).*sqrt(9*k^2*r(r>rs).^3 + 24*rs^3*rho)) ) )    );

end
