% VAINSHTEINONE Return the Vainshtein field of a single spherical body.
%
% phi_nd = VAINSHTEINONE(k, rho_nd, rs_nd, r_nd, [use_hypergeom], [reverse])
%
% VAINSHTEINONE returns the nondimensionalized Vainshtein field of a
% spherical body with nondimensional density rho_nd and nondimensional radius
% rs_nd, where the equation has constant k.
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
% Computing hypergeometric functions in MATLAB is very slow, so by default we use
% a patched polynomial that is good to within 0.1% everywhere.
% If use_hypergeom is passed as true, then hypergeometric functions will be used.
% NOTE: if symbolic toolbox is not installed, hypergeom is not available and will
% not be used.
%
% If reverse=false, the main branch will be returned (default). Else, the inverted
% (branch 2) solution will be returned. Note that the branch 2 solution is unique
% only up to a constant.
%
% If rho = 0, then 0 is returned.
%
% If k = 0, then the nonlinear-only analytic solution is returned. Because this
% solution goes to infinity as r goes to infinity, phi is instead normalized
% to phi(0) = 0.

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2018-06-20
% Updated: 2018-10-10

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


function phi = vainshteineone(k, rho, rs, r, use_hypergeom, reverse)

    if (nargin < 6)
        reverse = false;
        if (nargin < 5)
            use_hypergeom = false;
        end
    end
    assert(nargin >= 4);

    usage_str = 'vainshteineone(k, rho_nd, rs_nd, r_nd, [use_hypergeom], [reverse])';

    % Note: COMSOL makes calls with giant vectors for no reason, so we can't do these assertions....
    %assert(max(size(k)) == 1, sprintf('%s\n%s',usage_str, 'k must be a scalar.'));
    %assert(max(size(rho)) == 1, sprintf('%s\n%s',usage_str, 'rho_nd must be a scalar.'));
    %assert(max(size(rs)) == 1, sprintf('%s\n%s',usage_str, 'rs_nd must be a scalar.'));
    %assert(max(size(use_hypergeom)) == 1, sprintf('%s\n%s',usage_str, 'use_hypergeom must be a boolean.'));

    k = k(1);
    rho = rho(1);
    rs = rs(1);
    use_hypergeom = use_hypergeom(1);
    assert(islogical(use_hypergeom) || use_hypergeom == 0 || use_hypergeom == 1, sprintf('%s\n%s',usage_str, 'use_hypergeom must be a boolean.'));
    reverse = reverse(1);
    assert(islogical(reverse) || reverse == 0 || reverse == 1, sprintf('%s\n%s',usage_str, 'reverse must be a boolean.'));


    toolbox_ver = ver;
    has_symtools = license('test', 'symbolic_toolbox') && any(strcmp(cellstr(char(toolbox_ver.Name)), 'Symbolic Math Toolbox'));

    if (use_hypergeom && ~has_symtools)
        warning('Symbolic toolbox is not installed, so use_hypergeom will be ignored.');
        use_hypergeom = false;
    end

    % If rho == 0, return 0
    if (rho == 0)
        phi = 0*r;
        return
    end

    % If k == 0, return the nonlinear-only analytic solution
    if (k == 0)
        phi = sqrt(rho/24)*( (r.^2).*(r<=rs) + rs^(3/2)*(4*sqrt(r) - 3*sqrt(rs)).*(r>rs) );
        return
    end

    % rho ~= 0 and k ~= 0, so return the full solution

    % Inside the body it's a polynomial

    if (use_hypergeom)
        phi_in = (1/24)*( -3*k*r.^2 + 3*k*rs^2*hypergeom([-0.5,-2/3],1/3,-(8*rho)/(3*k^2) ) + (r.^2 - rs^2)*sqrt(9*k^2 + 24*rho) ) .* (r <= rs);
    else
        % without symbolic toolbox, not "hypergeom" function
        phi_in = (1/24)*( 3*k*rs^2*hypergeom_term(k^(2/3)/rho^(1/3)) + (r.^2 - rs^2)*sqrt(9*k^2 + 24*rho) ) .* (r <= rs);
    end

    if (use_hypergeom)
        phi = ( (k/8) * (max(r,rs).^2) .* (-1 + hypergeom([-0.5,-2/3],1/3,-(8*rs^3*rho)/(3*k^2) * max(r,rs).^(-3) ) ) ) .* (r > rs) + phi_in;
    else
        phi = ( (k/8) * (max(r,rs).^2) .* hypergeom_term((k^(2/3))/(rs*rho^(1/3)) * max(r,rs)) ) .* (r > rs) + phi_in;
    end

    if (reverse)
        phi = - (k/4)*r.^2 - phi;
    end

end


% Approximation of hypergeom([-0.5, -2/3], 1/3, -(8/3) * x.^(-3) )
function hcombo = hypergeom_term(x)

    % <0.1% error from 0 to 1
    h2 = -2*gamma(1/6)*gamma(1/3)/(3^(2/3)*sqrt(pi)) * x.^(-2) + 8 *sqrt(2/3)*x.^(-3/2) - 1 + sqrt(3/2)/7 * x.^(3/2) - (3/416)*sqrt(3/2).*x.^(9/2);

    % <0.1% error from 1 to 1.6
    h3 = h2 + 0.00591542 + 0.0305915*(x - 1.3) + 0.0645799*(x - 1.3).^2 + 0.0686049*(x - 1.3).^3;

    % <0.1% error from 1.55 to infinity
    h4 = -(8/3)*x.^(-3) + (4/9)*x.^(-6) - (64/189)*x.^(-9) + (32/81)*x.^(-12);

    hcombo = h2.*(x<1) + h3.*(x >= 1).*(x<1.6) + h4.*(x>=1.6);
end
