% NEAREST_INDEX Return the index and value of the point in the array nearest each number in input_vals.
% In case of ties, returns FIRST matching value in the array.
%
% [nindex, nvals] = NEAREST_INDEX(the_array, input_vals, (method));
%
% the_array must be a 1D array of values.
% input_vals can be a number or an array.
% method may be 'nearest' (default), 'below', or 'above'
% 'nearest' returns the nearest value, as one would expect.
% 'below' and 'above' return the nearest value below or above the provided
% value. Note: If the value is beyond the endpoints, then the nearest value is
% returned instead.
%
%
% Example:
%    >> x = [11 12 13 14 15];
%    >> [nindex, nvals] = NEAREST_INDEX(x, [13.14 12.99 0 20 14.5])
%    nindex = [ 3    3    1    5   4]
%    nvals =  [13   13   11   15  14]
%
%    >> [nindex, nvals] = NEAREST_INDEX(x, [13.14 12.99 0 20 14.5], 'below')
%    nindex = [ 3    2    1    5   4]
%    nvals =  [13   12   11   15  14]
%
%    >> [nindex, nvals] = NEAREST_INDEX(x, [13.14 12.99 0 20 14.5], 'above')
%    nindex = [ 4    3    1    5   5]
%    nvals =  [14   13   11   15  15]

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2016-08-04
% Updated: 2019-06-17

function [nindex, nvals] = nearest_index(the_array, input_vals, method);

usage_str = 'Usage: nearest_index(the_array, input_vals, [method])';
assert((nargin == 2) || (nargin == 3), usage_str);

if (nargin == 2)
    method = 'nearest';
end

assert(strcmpi(method, 'nearest') || strcmpi(method, 'above') || strcmpi(method, 'below'), 'method (third argument) must be ''nearest'', ''above'', or ''below''.')

nindex = input_vals*0;
nvals = input_vals*NaN;

assert(min(size(the_array))==1, 'First argument (array to search) must be a vector.');


[rin cin] = size(input_vals);

arr_min = min(the_array);
arr_max = max(the_array);

for r = [1:rin];
    for c = [1:cin];
        input_val = input_vals(r,c);

        if strcmpi(method, 'nearest') || (input_val < arr_min) || (input_val > arr_max)
            [~, besti] = min(abs(input_val - the_array));
        elseif strcmpi(method, 'above')
            [~, besti] = min(abs(input_val - the_array) + (arr_max - arr_min + 1)*(the_array < input_val));
        elseif strcmpi(method, 'below')
            [~, besti] = min(abs(input_val - the_array) + (arr_max - arr_min + 1)*(the_array > input_val));
        else
            assert(false, 'method (third argument) must be ''nearest'', ''above'', or ''below''.')
        end
        nindex(r,c) = besti;
        nvals(r,c) = the_array(besti);

    end;
end;
