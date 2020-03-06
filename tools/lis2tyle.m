% LIS2TYLE Give a plot the default LIS2T style.
%
% LIS2TYLE(fig_h)
% LIS2TYLE(fig_h, 'Preset', 'print')
% LIS2TYLE(fig_h, 'Preset', 'print', 'AutoArrange', false)
% LIS2TYLE(fig_h, 'Preset', 'print', 'PreserveTicks', false)
%
% LIS2TYLE can take many parameters as arguments. Any parameters not passed
% explicitly will be set by the choice of preset (default: print). If you want
% to avoid setting a variable, but also don't want to pass a specific value,
% then pass the parameter as 'ignore'.
% For example, if you set various markers with different MarkerSize and don't
% want them all to be overridden by a single value, you could call LIS2TYLE as
% LIS2TYLE(fig_h, 'MarkerSize', 'ignore')
%
% To ignore just a single element in a plot, give it the tag 'lis2tyle_ignore':
% set(obj, 'Tag', 'lis2tyle_ignore');
%
% The 'none' preset doesn't change any default values; this is useful for only
% calling AutoArrange and nothing else.
%
% A second feature, outside of the preset parameters, is AutoArrange.
% If AutoArrange is set to true, then LIS2TYLE will adjust the spacing of the
% plot slightly. For example, it will give the colorbar and some of the labels
% a bit more room.
% AutoArrange is false by default.
%
% If PreserveTicks is set to true, then the existing tickmarks and labels will
% be preserved, even if increasing the font would otherwise make MATLAB reduce
% their number.
%
% Available presets: 'print', 'jcap', 'none'
%
% Parameters: Parameter           | 'print' preset value | description
%             Units               | 'inches'             | units of measure for width and height
%             Width               | 5.6                  | width
%             Height              | 4.2                  | height
%             AxesWidth           | 'ignore'             | width of axes (takes precedence over figure)
%             AxesHeight          | 'ignore'             | height of axes (takes precedence over figure)
%             AxesLineWidth       | 0.75                 |
%             LineWidth           | 1.5                  | width for plot liens
%             MarkerSize          | 8                    | size of markers in the plot
%             Box                 | 'on'                 | 'on' or 'off' (should there be axes on all sides)
%             AxesNumFont         | 'Myriad Pro'         | font for numbers along axes
%             AxesNumFontSize     | 10                   |
%             AxesNumFontWeight   | 'normal'             | 'light', 'normal', 'demi', or 'bold'
%             AxesNumFontAngle    | 'normal'             | 'normal', 'italic', or 'oblique'
%             AxesLabelFont       | 'Myriad Pro'         | font for axes labels (e.g. 'x' or 'y' beside the axis)
%             AxesLabelFontSize   | 12                   |
%             AxesLabelFontWeight | 'normal'             | 'light', 'normal', 'demi', or 'bold'
%             AxesLabelFontAngle  | 'normal'             | 'normal', 'italic', or 'oblique'
%             AxesTextFont        | 'Myriad Pro'         | font for other text objects appearing in axes (e.g. annotations)
%             AxesTextFontSize    | 12                   |
%             AxesTextFontWeight  | 'normal'             | 'light', 'normal', 'demi', or 'bold'
%             AxesTextFontAngle   | 'normal'             | 'normal', 'italic', or 'oblique'
%             LegendFont          | 'Myriad Pro'         | font for legends
%             LegendFontSize      | 11                   |
%             LegendFontWeight    | 'normal'             | 'light', 'normal', 'demi', or 'bold'
%             LegendFontAngle     | 'normal'             | 'normal', 'italic', or 'oblique'
%             TitleFont           | 'Myriad Pro'         | font for plot title
%             TitleFontSize       | 12                   |
%             TitleFontWeight     | 'normal'             | 'light', 'normal', 'demi', or 'bold'
%             TitleFontAngle      | 'normal'             | 'normal', 'italic', or 'oblique'
%             TickDir             | 'in'                 | 'in' or 'out' (ticks overlap the plot or not)
%             TickLength          | 'ignore'             | [2D length, 3D length]
%             ColorbarTickLength  | 'ignore'             | length as % of height
%             ColorbarTickLengthW | 'ignore'             | length as % of width
%             ColorbarOffset      | 'ignore'             | distance from axis to colorbar
%             XMinorTick          | 'on'                 | 'on' or 'off'
%             YMinorTick          | 'on'                 | 'on' or 'off'
%             ZMinorTick          | 'on'                 | 'on' or 'off'
%

% Inspired in part by: https://dgleich.wordpress.com/2013/06/04/creating-high-quality-graphics-in-matlab-for-papers-and-presentations/

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2017-08-08
% Updated: 2018-08-22

function lis2tyle(varargin)

%%%%%%%%%% Parse input %%%%%%%%%%

parser = inputParser;
parser.addRequired( 'fig_h', @(x) isa(x, 'matlab.ui.Figure'));
parser.addParameter('Preset', 'print', @(x) (ischar(x)));
parser.addParameter('AutoArrange', false, @(x) (islogical(x) && length(x) == 1));
parser.addParameter('PreserveTicks', false, @(x) (islogical(x) && length(x) == 1));
parser.addParameter('XPaperMargins', [0 0], @(x) (isnumeric(x) && length(x) == 2));
parser.addParameter('YPaperMargins', [0 0], @(x) (isnumeric(x) && length(x) == 2));

% NOTE when adding parameters: Currently there's an assertion to make sure all numeric parameters are > 0
numeric_parameters = {...
                        'Width',...
                        'Height',...
                        'AxesWidth',...
                        'AxesHeight',...
                        'AxesLineWidth',...
                        'LineWidth',...
                        'MarkerSize',...
                        'AxesNumFontSize',...
                        'AxesLabelFontSize',...
                        'AxesTextFontSize',...
                        'LegendFontSize',...
                        'TitleFontSize',...
                        'ColorbarTickLength',...
                        'ColorbarTickLengthW',...
                        'ColorbarOffset',...
                      };

char_parameters = {...
                        'Units',...
                        'Box',...
                        'AxesNumFont',...
                        'AxesNumFontWeight',...
                        'AxesNumFontAngle',...
                        'AxesLabelFont',...
                        'AxesLabelFontWeight',...
                        'AxesLabelFontAngle',...
                        'AxesTextFont',...
                        'AxesTextFontWeight',...
                        'AxesTextFontAngle',...
                        'LegendFont',...
                        'LegendFontWeight',...
                        'LegendFontAngle',...
                        'TitleFont',...
                        'TitleFontWeight',...
                        'TitleFontAngle',...
                        'TickDir',...
                        'XMinorTick',...
                        'YMinorTick',...
                        'ZMinorTick',...
                      };

numeric_list_parameters = {...
                        'TickLength',...
                      };


funcchar_parameters = {...
                        'XTickFormat',...
                        'YTickFormat',...
                        'ZTickFormat',...
                      };

all_parameter_names = [numeric_parameters, char_parameters, numeric_list_parameters, funcchar_parameters];

for j = 1:length(numeric_parameters);
    param = numeric_parameters{j};
    parser.addParameter(param, [], @(x) ((isnumeric(x) && (length(x)==1) && (x>0)) || (ischar(param_val) && strcmpi(param_val, 'ignore'))));
end
for j = 1:length(char_parameters);
    param = char_parameters{j};
    parser.addParameter(param, [], @(x) (ischar(x)));
end
for j = 1:length(numeric_list_parameters);
    param = numeric_list_parameters{j};
    parser.addParameter(param, [], @(x) ((isnumeric(x)) || (ischar(param_val) && strcmpi(param_val, 'ignore'))));
end
for j = 1:length(funcchar_parameters);
    param = funcchar_parameters{j};
    parser.addParameter(param, [], @(x) (ischar(x) || isa(x, 'function_handle')));
end

parse(parser, varargin{:});
fig_h = parser.Results.fig_h;
preset = parser.Results.Preset;
auto_arrange = parser.Results.AutoArrange;
preserve_ticks = parser.Results.PreserveTicks;
xpaper_margins = parser.Results.XPaperMargins;
ypaper_margins = parser.Results.YPaperMargins;

parameters = containers.Map;
for j = 1:length(all_parameter_names)
    param = all_parameter_names{j};
    parameters(param) = parser.Results.(param);
end

assert_strval(preset, {'print', 'none', 'jcap'}, true, 'Preset', 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Set up default values based on presets %%%%%%%%%%

default_parameters = containers.Map;
if strcmpi(preset, 'print')
    default_parameters('Units')               = 'inches';
    default_parameters('Width')               = 5.6;
    default_parameters('Height')              = 4.2;
    default_parameters('AxesWidth')           = 'ignore';
    default_parameters('AxesHeight')          = 'ignore';
    default_parameters('Box')                 = 'on';
    default_parameters('AxesLineWidth')       = 0.75;
    default_parameters('LineWidth')           = 1.5;
    default_parameters('MarkerSize')          = 8;
    default_parameters('AxesNumFont')         = 'Myriad Pro';
    default_parameters('AxesNumFontSize')     = 10;
    default_parameters('AxesNumFontWeight')   = 'normal';
    default_parameters('AxesNumFontAngle')    = 'normal';
    default_parameters('AxesLabelFont')       = 'Myriad Pro';
    default_parameters('AxesLabelFontSize')   = 12;
    default_parameters('AxesLabelFontWeight') = 'normal';
    default_parameters('AxesLabelFontAngle')  = 'normal';
    default_parameters('AxesTextFont')        = 'Myriad Pro';
    default_parameters('AxesTextFontSize')    = 12;
    default_parameters('AxesTextFontWeight')  = 'normal';
    default_parameters('AxesTextFontAngle')   = 'normal';
    default_parameters('LegendFont')          = 'Myriad Pro';
    default_parameters('LegendFontSize')      = 11;
    default_parameters('LegendFontWeight')    = 'normal';
    default_parameters('LegendFontAngle')     = 'normal';
    default_parameters('TitleFont')           = 'Myriad Pro';
    default_parameters('TitleFontSize')       = 12;
    default_parameters('TitleFontWeight')     = 'normal';
    default_parameters('TitleFontAngle')      = 'normal';
    default_parameters('TickDir')             = 'in';
    default_parameters('TickLength')          = 'ignore';
    default_parameters('ColorbarTickLength')  = 'ignore';
    default_parameters('ColorbarTickLengthW') = 'ignore';
    default_parameters('ColorbarOffset')      = 'ignore';
    default_parameters('XMinorTick')          = 'on';
    default_parameters('YMinorTick')          = 'on';
    default_parameters('ZMinorTick')          = 'on';
    default_parameters('XTickFormat')         = 'ignore';
    default_parameters('YTickFormat')         = 'ignore';
    default_parameters('ZTickFormat')         = 'ignore';
elseif strcmpi(preset, 'jcap')
    default_parameters('Units')               = 'cm';
    default_parameters('Width')               = 18.125;
    default_parameters('Height')              = 14.5;
    default_parameters('AxesWidth')           = 'ignore';
    default_parameters('AxesHeight')          = 'ignore';
    default_parameters('Box')                 = 'on';
    default_parameters('AxesLineWidth')       = 1.5;
    default_parameters('LineWidth')           = 3;
    default_parameters('MarkerSize')          = 16;
    default_parameters('AxesNumFont')         = 'Times New Roman';
    default_parameters('AxesNumFontSize')     = 25;
    default_parameters('AxesNumFontWeight')   = 'normal';
    default_parameters('AxesNumFontAngle')    = 'normal';
    default_parameters('AxesLabelFont')       = 'Times New Roman';
    default_parameters('AxesLabelFontSize')   = 27;
    default_parameters('AxesLabelFontWeight') = 'normal';
    default_parameters('AxesLabelFontAngle')  = 'normal';
    default_parameters('AxesTextFont')        = 'Times New Roman';
    default_parameters('AxesTextFontSize')    = 25;
    default_parameters('AxesTextFontWeight')  = 'normal';
    default_parameters('AxesTextFontAngle')   = 'normal';
    default_parameters('LegendFont')          = 'Times New Roman';
    default_parameters('LegendFontSize')      = 25;
    default_parameters('LegendFontWeight')    = 'normal';
    default_parameters('LegendFontAngle')     = 'normal';
    default_parameters('TitleFont')           = 'Times New Roman';
    default_parameters('TitleFontSize')       = 25;
    default_parameters('TitleFontWeight')     = 'normal';
    default_parameters('TitleFontAngle')      = 'normal';
    default_parameters('TickDir')             = 'in';
    default_parameters('TickLength')          = [0.018, 0.018];
    default_parameters('ColorbarTickLength')  = 'ignore';
    default_parameters('ColorbarTickLengthW') = 0.5;
    default_parameters('ColorbarOffset')      = 0.6;
    default_parameters('XMinorTick')          = 'on';
    default_parameters('YMinorTick')          = 'on';
    default_parameters('ZMinorTick')          = 'on';
    default_parameters('XTickFormat')         = 'ignore';
    default_parameters('YTickFormat')         = 'ignore';
    default_parameters('ZTickFormat')         = 'ignore';
elseif strcmpi(preset, 'none')
    default_parameters('Units')               = 'ignore';
    default_parameters('Width')               = 'ignore';
    default_parameters('Height')              = 'ignore';
    default_parameters('AxesWidth')           = 'ignore';
    default_parameters('AxesHeight')          = 'ignore';
    default_parameters('Box')                 = 'ignore';
    default_parameters('AxesLineWidth')       = 'ignore';
    default_parameters('LineWidth')           = 'ignore';
    default_parameters('MarkerSize')          = 'ignore';
    default_parameters('AxesNumFont')         = 'ignore';
    default_parameters('AxesNumFontSize')     = 'ignore';
    default_parameters('AxesNumFontWeight')   = 'ignore';
    default_parameters('AxesNumFontAngle')    = 'ignore';
    default_parameters('AxesLabelFont')       = 'ignore';
    default_parameters('AxesLabelFontSize')   = 'ignore';
    default_parameters('AxesLabelFontWeight') = 'ignore';
    default_parameters('AxesLabelFontAngle')  = 'ignore';
    default_parameters('AxesTextFont')        = 'ignore';
    default_parameters('AxesTextFontSize')    = 'ignore';
    default_parameters('AxesTextFontWeight')  = 'ignore';
    default_parameters('AxesTextFontAngle')   = 'ignore';
    default_parameters('LegendFont')          = 'ignore';
    default_parameters('LegendFontSize')      = 'ignore';
    default_parameters('LegendFontWeight')    = 'ignore';
    default_parameters('LegendFontAngle')     = 'ignore';
    default_parameters('TitleFont')           = 'ignore';
    default_parameters('TitleFontSize')       = 'ignore';
    default_parameters('TitleFontWeight')     = 'ignore';
    default_parameters('TitleFontAngle')      = 'ignore';
    default_parameters('TickDir')             = 'ignore';
    default_parameters('TickLength')          = 'ignore';
    default_parameters('ColorbarTickLength')  = 'ignore';
    default_parameters('ColorbarTickLengthW') = 'ignore';
    default_parameters('ColorbarOffset')      = 'ignore';
    default_parameters('XMinorTick')          = 'ignore';
    default_parameters('YMinorTick')          = 'ignore';
    default_parameters('ZMinorTick')          = 'ignore';
    default_parameters('XTickFormat')         = 'ignore';
    default_parameters('YTickFormat')         = 'ignore';
    default_parameters('ZTickFormat')         = 'ignore';
% TODO: more presets
end

for j = 1:length(all_parameter_names)
    param = all_parameter_names{j};
    if length(parameters(param)) == 0
        parameters(param) = default_parameters(param);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%% Validate input %%%%%%%%%%

units = parameters('Units');
assert_strval(units, {'inches', 'in', 'centimeters', 'cm', 'points', 'pt', 'pixels', 'px', 'normalized','ignore'}, true, 'Units', 2);
if strcmpi(units,'in'); units = 'inches'; end;
if strcmpi(units,'cm'); units = 'centimeters'; end;
if strcmpi(units,'pt'); units = 'points'; end;
if strcmpi(units,'px'); units = 'pixels'; end;
parameters('Units') = units;

assert_strval(parameters('TickDir'), {'out','in','ignore'}, true, 'TickDir', 2);
assert_strval(parameters('XMinorTick'), {'on','off','ignore'}, true, 'XMinorTick', 2);
assert_strval(parameters('YMinorTick'), {'on','off','ignore'}, true, 'YMinorTick', 2);
assert_strval(parameters('ZMinorTick'), {'on','off','ignore'}, true, 'ZMinorTick', 2);

if isa(parameters('XTickFormat'), 'function_handle')
    xtickformat_f = parameters('XTickFormat');
elseif strcmpi(parameters('XTickFormat'), 'ignore')
    xtickformat_f = @(num) num2str(num);
else
    xtickformat_f = @(num) sprintf(parameters('XTickFormat'), num);
end
if isa(parameters('YTickFormat'), 'function_handle')
    ytickformat_f = parameters('YTickFormat');
elseif strcmpi(parameters('YTickFormat'), 'ignore')
    ytickformat_f = @(num) num2str(num);
else
    ytickformat_f = @(num) sprintf(parameters('YTickFormat'), num);
end
if isa(parameters('ZTickFormat'), 'function_handle')
    ztickformat_f = parameters('ZTickFormat');
elseif strcmpi(parameters('ZTickFormat'), 'ignore')
    ztickformat_f = @(num) num2str(num);
else
    ztickformat_f = @(num) sprintf(parameters('ZTickFormat'), num);
end

%TODO FINISH VALIDATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Get figure info and set dimensions %%%%%%%%%%

% get all axes and legends of the figure
axes_hs = [];
legend_hs = [];
fig_children = fig_h.Children();
for j = 1:length(fig_children)
    fig_child = fig_children(j);
    if hasIgnoreTag(fig_child); continue; end;

    if isa(fig_child, 'matlab.graphics.axis.Axes')
        if (length(axes_hs) == 0)
            axes_hs = [fig_child];
        else
            axes_hs(end+1) = fig_child;
        end
    elseif isa(fig_child, 'matlab.graphics.illustration.Legend')
        if (length(legend_hs) == 0)
            legend_hs = [fig_child];
        else
            legend_hs(end+1) = fig_child;
        end
    end
end
axesN = length(axes_hs);
legendN = length(legend_hs);

% Get all annotation axes of the figure
% (hidden method: https://www.mathworks.com/matlabcentral/answers/167708-accessing-the-annotation-handle )
annotation_axes_hs = [];
annotation_axes_hs = findall(fig_h, 'Tag', 'scribeOverlay');
annotation_axesN = length(annotation_axes_hs);

% Save ticks of every axis.
if preserve_ticks
    axes_Xticks = {};
    axes_XtickLabels = {};
    axes_Yticks = {};
    axes_YtickLabels = {};
    axes_Zticks = {};
    axes_ZtickLabels = {};
    for ax = axes_hs
        axes_Xticks{end+1} = ax.XTick;
        axes_XtickLabels{end+1} = ax.XTickLabel;
        axes_Yticks{end+1} = ax.YTick;
        axes_YtickLabels{end+1} = ax.YTickLabel;
        axes_Zticks{end+1} = ax.ZTick;
        axes_ZtickLabels{end+1} = ax.ZTickLabel;
    end
end

% set figure window size to match the size we'll print later
%NOTE: width and height still have to be set again when we print the figure
setOrIgnore(fig_h, 'Units', units);
pos = get(fig_h, 'Position');
width = replaceIgnore(parameters('Width'), pos(3));
height = replaceIgnore(parameters('Height'), pos(4));
dheight = height - pos(4);
set(fig_h, 'Position', [(pos(1)) (pos(2)-dheight) width height]);

% Set printing properties to have the same width and height
setOrIgnore(fig_h, 'PaperUnits', units);

set(fig_h, 'PaperSize', [width*(1+sum(xpaper_margins)) height*(1+sum(ypaper_margins))]);
fig_h.PaperPosition = [0 0 width*(1+sum(xpaper_margins)) height*(1+sum(ypaper_margins))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% AutoArrange plot %%%%%%%%%%

if (auto_arrange)

    cb = findall(fig_h, 'Tag', 'Colorbar');

    axpos_orig = {};
    for j = 1:axesN
        axes_h.Units = 'normalized';
        axes_h = axes_hs(j);
        axpos_orig{j} = axes_h.Position;
    end

    % Move the colorbar label right, to avoid overlapping the tick labels.
    % Then, move the colorbar left to compensate so that the label is within
    % the figure bounds.
    label_dx = 0;
    if (length(cb) == 1) && ~hasIgnoreTag(cb)
        cbtick_pts = cb.FontSize;
        cb_units_orig = cb.Units;
        cblabel_pts = cb.Label.FontSize;
        cblabel_units_orig = cb.Label.Units;
       
        % Move the label to the right by cbtickpts/2 points
        cb.Label.Units = 'points';
        cblabel_pos = cb.Label.Position;
        cb.Label.Position = cblabel_pos + [(0.5*cbtick_pts) 0 0];
        cblabel_pos = cb.Label.Position;

        % Now we have to figure out how far off the edge we've gone.
        fig_units_orig = fig_h.Units;
        fig_h.Units = 'points';
        fig_width = fig_h.Position(3);

        cb.Units = 'points';
        cb_pos = cb.Position;

        buffer_pts = 6;
        label_dx = cb_pos(1) + cblabel_pos(1) + cblabel_pts + buffer_pts - fig_width;

        % Move the colorbar left by label_dx
        cb.Position = cb_pos + [-label_dx, 0, 0, 0];

        % Move the colorbar an extra 1% left, bringing it closer to the axes.
        cb.Units = 'normalized';
        cb_pos = cb.Position;
        cb.Position = cb_pos + [-0.01, 0, 0, 0];

        % Move the colorbar an extra 1% up, giving x label more room
        cb.Units = 'normalized';
        cb_pos = cb.Position;
        cb.Position = cb_pos + [0, 0.01, 0, 0];


        cb.Units = cb_units_orig;
        cb.Label.Units = cblabel_units_orig;
        fig_h.Units = fig_units_orig;
    end

    % Changing cb screws up the axes position, so reset it here.
    for j = 1:axesN
        axes_h = axes_hs(j);
        axes_h.Position = axpos_orig{j};
    end

    % Now adjust axes positions.
    for j = 1:axesN
        axes_h = axes_hs(j);

        axes_units_orig = axes_h.Units;
        axes_h.Units = 'points';

        axpos_pts = axes_h.Position;
        % adjust to accommodate colorbar
        axes_h.Position = axpos_pts + [0, 0, -label_dx, 0];


        % Move all axes up slightly (giving x label more room)
        %axes_h.Units = 'normalized';
        label_dy = replaceIgnore(parameters('AxesLabelFontSize'), 10);
        axpos = axes_h.Position;
        axes_h.Position = axpos + [0 1 0 0]*label_dy*0.5;
        axpos = axes_h.Position;

        % If there's no colorbar and only one axes, then move axes and xlabel
        % up together some more.
        if ((length(cb) == 0) && (length(axes_h)==1))
            % check that axes_h has XLabel property, and that it is nonempty
            if (length(axes_h.findprop('XLabel')) > 0)
                xlabel_h = axes_h.XLabel;
                if ( length(xlabel_h.String) > 0)
                    axes_h.Position = axpos + [0 0 0 -0.02];
                    axpos = axes_h.Position;

                    xlabel_units_orig = xlabel_h.Units;
                    xlabel_h.Units = 'normalized';
                    xlpos = xlabel_h.Position;
                    xlabel_h.Position = xlpos + [0 0.01 0];
                    xlpos = xlabel_h.Position;
                    xlabel_h.Units = xlabel_units_orig;
                end
            end
        end

        axes_h.Units = axes_units_orig;

    end

    % TODO: more arrangement optimizations

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Set plot style %%%%%%%%%%

% Loop through all the axes and set parameters
for j = 1:axesN
    axes_h = axes_hs(j);
    setOrIgnore(axes_h...
           ,'Box',        parameters('Box')...
           ,'LineWidth',  parameters('AxesLineWidth')...
           ,'FontName',   parameters('AxesNumFont')...
           ,'FontSize',   parameters('AxesNumFontSize')...
           ,'FontWeight', parameters('AxesNumFontWeight')...
           ,'FontAngle',  parameters('AxesNumFontAngle')...
           ,'TickDir',    parameters('TickDir')...
           ,'TickLength', parameters('TickLength')...
           ,'XMinorTick', parameters('XMinorTick')...
           ,'YMinorTick', parameters('YMinorTick')...
           ,'ZMinorTick', parameters('ZMinorTick')...
            );

    title_h = axes_h.Title;
    setOrIgnore(title_h...
           ,'FontName',   parameters('TitleFont')...
           ,'FontSize',   parameters('TitleFontSize')...
           ,'FontWeight', parameters('TitleFontWeight')...
           ,'FontAngle',  parameters('TitleFontAngle')...
            );

    xlabel_h = axes_h.XLabel;
    setOrIgnore(xlabel_h...
           ,'FontName',   parameters('AxesLabelFont')...
           ,'FontSize',   parameters('AxesLabelFontSize')...
           ,'FontWeight', parameters('AxesLabelFontWeight')...
           ,'FontAngle',  parameters('AxesLabelFontAngle')...
            );

    ylabel_h = axes_h.YLabel;
    setOrIgnore(ylabel_h...
           ,'FontName',   parameters('AxesLabelFont')...
           ,'FontSize',   parameters('AxesLabelFontSize')...
           ,'FontWeight', parameters('AxesLabelFontWeight')...
           ,'FontAngle',  parameters('AxesLabelFontAngle')...
            );

        % TODO Z axis (check if exists though??)

    axes_children = axes_h.Children();
    for k = 1:length(axes_children)
        axes_child = axes_children(k);
        if hasIgnoreTag(axes_child); continue; end;

        % lines
        if isa(axes_child, 'matlab.graphics.chart.primitive.Line')
            line_h = axes_child;
            setOrIgnore(line_h...
                   ,'LineWidth', parameters('LineWidth')...
                   ,'MarkerSize', parameters('MarkerSize')...
                    );
        % extra text
        elseif isa(axes_child, 'matlab.graphics.primitive.Text')
            text_h = axes_child;
            setOrIgnore(text_h...
                   ,'FontName',   parameters('AxesTextFont')...
                   ,'FontSize',   parameters('AxesTextFontSize')...
                   ,'FontWeight', parameters('AxesTextFontWeight')...
                   ,'FontAngle',  parameters('AxesTextFontAngle')...
                    );
        % TODO: deal with more objects in axes
        end
    end
end

% Loop through all the legends and set parameters
for j = 1:legendN
    legend_h = legend_hs(j);
    if hasIgnoreTag(legend_h); continue; end;

    setOrIgnore(legend_h...
           ,'FontName',   parameters('LegendFont')...
           ,'FontSize',   parameters('LegendFontSize')...
           ,'FontWeight', parameters('LegendFontWeight')...
           ,'FontAngle',  parameters('LegendFontAngle')...
            );
end

% Loop through all the colorbars and set parameters
cb_hs = findall(fig_h, 'Tag', 'Colorbar');
for j = 1:length(cb_hs)
    cb_h = cb_hs(j);
    if hasIgnoreTag(cb_h); continue; end;

    setOrIgnore(cb_h...
           ,'Units',   parameters('Units')...
           ,'FontName',   parameters('LegendFont')...
           ,'FontSize',   parameters('LegendFontSize')...
           ,'FontWeight', parameters('LegendFontWeight')...
           ,'FontAngle',  parameters('LegendFontAngle')...
           ,'LineWidth',  parameters('AxesLineWidth')...
           ,'TickLength',  parameters('ColorbarTickLength')...
            );

    % set TickLength as width percentage
    if ~ischar(parameters('ColorbarTickLengthW'))
        cur_units = cb_h.Units;

        cb_h.Units = 'centimeters';
        cb_h.TickLength = parameters('ColorbarTickLengthW')*cb_h.Position(3)/cb_h.Position(4);

        cb_h.Units = cur_units;
    end
end

% Loop through all the annotations and set parameters
% We loop through both annotation_axes_hs and axes_hs, in case some annotations
% have been set as axes children.
for j = 1:annotation_axesN + axesN
    if (j <= annotation_axesN)
        axes_h = annotation_axes_hs(j);
    else
        axes_h = axes_hs(j - annotation_axesN);
    end
    if hasIgnoreTag(axes_h); continue; end;

    axes_children = axes_h.Children();
    for k = 1:length(axes_children)
        axes_child = axes_children(k);
        if hasIgnoreTag(axes_child); continue; end;

        % lines
        if (isa(axes_child, 'matlab.graphics.shape.Line') ...
           || isa(axes_child, 'matlab.graphics.shape.Arrow') ...
           || isa(axes_child, 'matlab.graphics.shape.DoubleArrow') ...
           || isa(axes_child, 'matlab.graphics.shape.TextArrow') ...
           || isa(axes_child, 'matlab.graphics.shape.Rectangle') ...
           || isa(axes_child, 'matlab.graphics.primitive.Rectangle') ...
           || isa(axes_child, 'matlab.graphics.shape.Ellipse') ...
           || isa(axes_child, 'matlab.graphics.shape.TextBox') )
            line_h = axes_child;
            setOrIgnore(line_h...
                   ,'LineWidth', parameters('LineWidth')...
                    );
        end
        % extra text
        if ( isa(axes_child, 'matlab.graphics.shape.TextArrow') ...
           || isa(axes_child, 'matlab.graphics.shape.TextBox') )
            text_h = axes_child;
            setOrIgnore(text_h...
                   ,'FontName',   parameters('AxesTextFont')...
                   ,'FontSize',   parameters('AxesTextFontSize')...
                   ,'FontWeight', parameters('AxesTextFontWeight')...
                   ,'FontAngle',  parameters('AxesTextFontAngle')...
                    );
        end
        % TODO: add more annotation variables
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% SET TICKS %%%%%%%%%%
% Set ticks of every axis.
for axi = 1:axesN
    ax = axes_hs(axi);
    if preserve_ticks
        ax.XTick = axes_Xticks{axi};
        ax.YTick = axes_Yticks{axi};
        ax.ZTick = axes_Zticks{axi};
        %ax.XTickLabel = axes_XtickLabels{axi};
        %ax.YTickLabel = axes_YtickLabels{axi};
        %ax.ZTickLabel = axes_ZtickLabels{axi};
    end
    ax.XTickLabel = cellfun(xtickformat_f, num2cell(ax.XTick), 'UniformOutput', false)';
    ax.YTickLabel = cellfun(ytickformat_f, num2cell(ax.YTick), 'UniformOutput', false)';
    ax.ZTickLabel = cellfun(ztickformat_f, num2cell(ax.ZTick), 'UniformOutput', false)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% SET AXES SIZE %%%%%%%%%%

if (axesN == 1)
    axes_h = axes_hs(1);
    setOrIgnore(axes_h, 'Units', units);

    apos = axes_h.Position;
    pos = fig_h.Position;
    ppos = fig_h.PaperPosition;
    psize = fig_h.PaperSize;

    ax_width = replaceIgnore(parameters('AxesWidth'), apos(3));
    ax_height = replaceIgnore(parameters('AxesHeight'), apos(4));
    ax_width_diff = ax_width - apos(3);
    ax_height_diff = ax_height - apos(4);

    set(axes_h, 'Position', [apos(1), apos(2), apos(3) + ax_width_diff, apos(4) + ax_height_diff]);
    set(fig_h, 'Position', [pos(1), pos(2), pos(3) + ax_width_diff, pos(4) + ax_height_diff]);
    set(fig_h, 'PaperSize', [psize(1) + ax_width_diff, psize(2) + ax_height_diff]);
    fig_h.PaperPosition = [ppos(1), ppos(2), ppos(3) + ax_width_diff, ppos(4) + ax_height_diff];

    % set colorbar offset
    cb_hs = findall(fig_h, 'Tag', 'Colorbar');
    if length(cb_hs) == 1
        cb_h = cb_hs(1);

        if ~ischar(parameters('ColorbarOffset'))
            cb_offset = parameters('ColorbarOffset');

            cb_h.Position(1) = axes_h.Position(1) + axes_h.Position(3) + cb_offset;

        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% FIX MARGINS %%%%%%%%%%

if (axesN == 1)
    axes_h = axes_hs(1);
    %setOrIgnore(axes_h, 'Units', units);
    apos = axes_h.Position;
    set(axes_h, 'Position', [(apos(1) + xpaper_margins(1)*apos(3)), (apos(2) + ypaper_margins(1)*apos(4)), apos(3), apos(4)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % end of main function

%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%

% Test if an object has the tag 'lis2tyle_ignore'
function has_ignore_tag = hasIgnoreTag(obj)
    has_ignore_tag = false;

    try % in case null or non-object variable is passed

        if obj.isprop('Tag')
            this_tag = obj.Tag;
            if ischar(this_tag) && strcmpi(this_tag, 'lis2tyle_ignore')
                has_ignore_tag = true;
            end
        end

    end

end

% Return val, unless it is 'ignore' in which case return default
function result = replaceIgnore(val, default)
    result = val;
    if ischar(val) && strcmpi(val, 'ignore')
        result = default;
    end
end

% Set object properties, unless the values given are 'ignore'
function setOrIgnore(obj, varargin);
    numargs = nargin;
    assert(numargs > 1);
    assert(mod(numargs,2) == 1);
    assert(mod(length(varargin),2) == 0); % should be equivalent to previous line
    
    for j = 1:length(varargin)/2
        param_name = varargin{2*j-1};
        param_val = varargin{2*j};
        if ~(ischar(param_val) && strcmpi(param_val, 'ignore'))
            set(obj, param_name, param_val);
        end
    end

end

% Check that str is in valid_strs, and give a warning or error if it's not.
% This is really long because it was copied from NCW's function library.
% It's just used here for parsing the input.
function valid = assert_strval(str, valid_strs, ignore_case, varname, warning_level);

    valid = false;
    assert(nargin == 2 || nargin == 3 || nargin == 4 || nargin == 5);
    assert(strcmp(class(valid_strs), 'cell'));
    assert(strcmp(class(str), 'char'));

    if (nargin <= 2);
        ignore_case = false;
    end;

    if (nargin >= 4);
        assert(strcmp(class(varname), 'char'));
    else;
        varname = inputname(1);
    end;

    if (nargin <= 4)
        warning_level = 1;
    end

    N = length(valid_strs);

    % Set the warning message based on the valid_strs cell.
    if (N == 0);
        message = [varname ' is invalid.'];
    elseif (N == 1);
        message = [varname ' must be ' valid_strs{1} '.'];
    elseif (N == 2);
        message = [varname ' must be ' valid_strs{1} ' or ' valid_strs{2} '.'];
    else;
        message = [varname ' must be ' strjoin({valid_strs{1:end-1}},', ') ', or ' valid_strs{end} '.'];
    end;

    % Set the string comparison function to ignore or not ignore case.
    if (ignore_case);
        strcmp_func = @(s1,s2) strcmpi(s1,s2);
    else;
        strcmp_func = @(s1,s2) strcmp(s1,s2);
    end;

    % Check the input string against each valid string.
    for j=1:N;
        this_str = valid_strs{j};
        if (strcmp_func(str,this_str));
            valid = true;
            break;
        end;
    end;

    if (~valid);
        if (warning_level == 0)
        elseif (warning_level == 1)
            warnStruct.message = message;
            warnStruct.identifier = 'assert_strval:failed';
            warnStruct.stack = dbstack(1);
            warning(warnStruct);
        else
            errorStruct.message = message;
            errorStruct.identifier = 'assert_strval:failed';
            errorStruct.stack = dbstack(1);
            error(errorStruct);
        end
    end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
