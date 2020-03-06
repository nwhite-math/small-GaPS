% CHECK_AND_DISP_PARAMS Check that variables are set, validate them, and display
%  their values.
%
% This is a script, NOT a function (due to scoping), so the input arguments must
% be set as variables before running.
% The input arguments are:
%  cadp_names        (required)
%  cadp_descriptions (optional)
%  cadp_attributes   (optional)
%
% cadp_names should be a cell of variable names, e.g. {'a', 'b', 'c'}.
% cadp_descriptions should be a cell of strings describing these variables (the
%  descriptions will be printed alongside the variables).
% cadp_attributes should be a Nx2 cell of attribute cells, which can be passed
% directly to validateattributes.
% If the first argument is 'validatestring', then the second will be passed
% to validatestring instead of validateattributes.
%
% Example cadp_attributes:
% { {'double'}, {'positive','2d'}; 'validatestring', {'option1', 'option2', 'option3'}}.
%
% WARNING: This will call 'eval' on the strings in cadp_names!

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2018-10-01
% Updated: 2018-10-01

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


% Check that variables exist
for cadp_i = 1:length(cadp_names);
    if (~exist(cadp_names{cadp_i},'var'));
        if exist('cadp_descriptions', 'var')
            error(['Variable ' cadp_names{cadp_i} ' (' cadp_descriptions{cadp_i} ') must be set.']);
        else
            error(['Variable ' cadp_names{cadp_i} ' must be set.']);
        end
    end;
end;

% Check that variables satisfy the correct attributes
% This will call eval on cadp_names, but it's after ensuring that they are all
% variables so it should be fine.
if exist('cadp_attributes', 'var')
    for cadp_i = 1:length(cadp_names);
        if ischar(cadp_attributes{cadp_i, 1}) && strcmp(cadp_attributes{cadp_i, 1}, 'validatestring')
            validatestring(eval(cadp_names{cadp_i}), cadp_attributes{cadp_i, 2}, '', cadp_names{cadp_i});
        else
            validateattributes(eval(cadp_names{cadp_i}), cadp_attributes{cadp_i, 1}, cadp_attributes{cadp_i, 2}, '', cadp_names{cadp_i});
        end
    end;
end

% Finally, print out the inputs

cadp_this_str = '...';
for cadp_i = 1:length(cadp_names);
    cadp_this_var = eval(cadp_names{cadp_i});
    if ischar(cadp_this_var)
        cadp_this_str = cadp_this_var;
    elseif iscell(cadp_this_var)
        cadp_this_str = ['{' num2str(size(cadp_this_var,1)) 'x' num2str(size(cadp_this_var,2)) ' cell}'];
        if (length(cadp_this_var) <= 5)
            try
                cadp_this_str = ['{' strjoin(cadp_this_var, ', ') '}'];
            end
        end
    else
        cadp_this_str =  ['[' num2str(size(cadp_this_var,1)) 'x' num2str(size(cadp_this_var,2)) ' array]'];
        try
            if (min(size(cadp_this_var)) == 1)
                if (max(size(cadp_this_var)) == 1)
                    cadp_this_str = num2str(cadp_this_var);
                elseif (length(cadp_this_var) <= 10)
                    cadp_this_str = ['[' num2str(cadp_this_var) ']'];
                end
            end
        end
    end

    if exist('cadp_descriptions', 'var')
        fprintf('%-20s %-20s (%s)\n', cadp_names{cadp_i}, cadp_this_str, cadp_descriptions{cadp_i});
    else
        fprintf('%-20s %-20s\n', cadp_names{cadp_i}, cadp_this_str);
    end
end

clear('cadp_i',...
      'cadp_j',...
      'cadp_this_var',...
      'cadp_this_str')
