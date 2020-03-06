% VISUALIZE_2D A variety of plotting functions for visualizing vain_findiff.m
% output. These functions may be called individually, or the wrapper function
% visualize_2d may be called with a list of plots to create.
%
% [axs, fs] = VISUALIZE_2D(meshes, idgrids, nodedata, plots_to_show, ...
%       'axs', axs, ...
%       'phiv', phiv, 'phi_initv', phi_initv, 'newton_psiv', newton_psiv, ...
%       'lap_phi_sumv', lap_phi_sumv, 'lap_phi_onevs', lap_phi_onevs, ...
%       'Resv', Resv, 'endensv', endensv, ...
%       'intres_hist', intres_hist, 'energy_hist', energy_hist, ...
%       'phiv_refs', phiv_refs, 'phi_ref_titles', phi_ref_titles, ...
%       'drv', drv, 'dzv', dzv, 'lapv', lapv, 'this_iter', this_iter, ...
%       'nBodies', nBodies, 'Rss', Rss, 'zbs', zbs, ...
%       'rrange_close', rrange_close, 'zrange_close', zrange_close ...
%       )
%
% VISUALIZE_2D creates all the plots listed in plots_to_show. All variables
% needed for each of those plots must be passed to VISUALIZE_2D as shown above.
% For example, if 'residual_history' is in the cell plots_to_show, then
% intres_hist must be passed as well. The variable axs should be a
% containers.Map with keys being plot types (e.g., 'residual_history') and
% values being axes. This way, plots can be updated in the same axes while
% vain_findiff.m is running.
%
% The valid plot types are:
%  'phi_difference', 'phi', 'initial_phi', 'residual', 'residual_history',
%  'energy_history', 'energy_density', 'log_energy_density',
%  'hiramatsu_screening_stat', 'laplacian_phi', 'newton_fraction',
%  'phi_difference_close', 'phi_close', 'initial_phi_close', 'residual_close',
%  'energy_density_close', 'log_energy_density_close',
%  'hiramatsu_screening_stat_close', 'laplacian_phi_close', 'newton_fraction_close'
%

% Author: Nicholas White <nwhite@caltech.edu>
%         Laboratory of Interfacial and Small Scale Transport (LIS2T)
%         Dept. of Applied Physics, Caltech
%         https://www.troian.caltech.edu
% Created: 2019-05-20
% Updated: 2019-05-21

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


function [axs, fs] = visualize_2d(varargin)

parser = inputParser;

parser.addRequired('meshes',          @ (x) iscell(x));
parser.addRequired('idgrids',         @ (x) iscell(x));
parser.addRequired('nodedata',        @ (x) isstruct(x));
parser.addRequired('plots_to_show',   @ (x) iscell(x));

parser.addParameter('axs',             [], @ (x) isa(x, 'containers.Map'));
parser.addParameter('phiv',            [], @ (x) isvector(x));
parser.addParameter('phi_initv',       [], @ (x) isvector(x));
parser.addParameter('newton_psiv',     [], @ (x) isvector(x));
parser.addParameter('lap_phi_sumv',    [], @ (x) isvector(x));
parser.addParameter('lap_phi_onevs',   {}, @ (x) iscell(x));
parser.addParameter('Resv',            [], @ (x) isvector(x));
parser.addParameter('endensv',         [], @ (x) isvector(x));
parser.addParameter('intres_hist',     [], @ (x) isvector(x));
parser.addParameter('energy_hist',     [], @ (x) isvector(x));
parser.addParameter('phiv_refs',       {}, @ (x) iscell(x));
parser.addParameter('phi_ref_titles',  {}, @ (x) iscell(x));
parser.addParameter('drv',             [], @ (x) ismatrix(x));
parser.addParameter('dzv',             [], @ (x) ismatrix(x));
parser.addParameter('lapv',            [], @ (x) ismatrix(x));
parser.addParameter('this_iter',       [], @ (x) isscalar(x));
parser.addParameter('nBodies',          0, @ (x) isscalar(x));
parser.addParameter('Rss',             [], @ (x) (isvector(x) || isempty(x)));
parser.addParameter('zbs',             [], @ (x) (isvector(x) || isempty(x)));
parser.addParameter('rrange_close',    [], @ (x) (isvector(x) && (length(x) == 2)));
parser.addParameter('zrange_close',    [], @ (x) (isvector(x) && (length(x) == 2)));

parse(parser, varargin{:});
meshes         = parser.Results.meshes;
idgrids        = parser.Results.idgrids;
nodedata       = parser.Results.nodedata;
plots_to_show  = parser.Results.plots_to_show;
axs            = parser.Results.axs;
phiv           = parser.Results.phiv;
phi_initv      = parser.Results.phi_initv;
newton_psiv    = parser.Results.newton_psiv;
lap_phi_sumv   = parser.Results.lap_phi_sumv;
lap_phi_onevs  = parser.Results.lap_phi_onevs;
Resv           = parser.Results.Resv;
endensv        = parser.Results.endensv;
intres_hist    = parser.Results.intres_hist;
energy_hist    = parser.Results.energy_hist;
phiv_refs      = parser.Results.phiv_refs;
phi_ref_titles = parser.Results.phi_ref_titles;
drv            = parser.Results.drv;
dzv            = parser.Results.dzv;
lapv           = parser.Results.lapv;
this_iter      = parser.Results.this_iter;
nBodies        = parser.Results.nBodies;
Rss            = parser.Results.Rss;
zbs            = parser.Results.zbs;
rrange_close   = parser.Results.rrange_close;
zrange_close   = parser.Results.zrange_close;


if isempty(this_iter)
    iter_str = '';
else
    iter_str = sprintf(' (iteration %d)', this_iter);
end

% Set colorscheme for Hiramatsu screening statistic plot (reverse jet)
persistent rjet;
if isempty(rjet)
    % set colormap to match Hiramatsu et al 2013
    f_jet = figure; rjet = jet; rjet=rjet(end:-1:1,:); close(f_jet);
end

% Filter out any invalid plot types
valid_plot_types = {...
'phi_difference', 'phi', 'initial_phi', 'residual', 'residual_history',...
'energy_history', 'energy_density', 'log_energy_density',...
'hiramatsu_screening_stat', 'laplacian_phi', 'newton_fraction',...
'phi_difference_close', 'phi_close', 'initial_phi_close', 'residual_close',...
'energy_density_close', 'log_energy_density_close',...
'hiramatsu_screening_stat_close', 'laplacian_phi_close',...
'newton_fraction_close'};

plots_to_show = intersect(valid_plot_types, plots_to_show);

% Only show plots with the correct variables defined
if isempty(phiv)
    plots_to_show = setdiff(plots_to_show, ...
    {'phi', 'phi_difference', 'laplacian_phi', ...
     'phi_close', 'phi_difference_close', 'laplacian_phi_close'});
end
if isempty(phi_initv)
    plots_to_show = setdiff(plots_to_show, {'initial_phi', 'initial_phi_close'});
end
if isempty(newton_psiv) || isempty(drv) || isempty(dzv)
    plots_to_show = setdiff(plots_to_show, {'newton_fraction', 'newton_fraction_close'});
end
if isempty(lap_phi_sumv) || isempty(lap_phi_onevs)
    plots_to_show = setdiff(plots_to_show, ...
    {'hiramatsu_screening_stat', 'hiramatsu_screening_stat_close'});
end
if isempty(Resv)
    plots_to_show = setdiff(plots_to_show, {'residual', 'residual_close'});
end
if isempty(endensv)
    plots_to_show = setdiff(plots_to_show, {'energy_density', 'energy_density_close'});
end
if isempty(intres_hist)
    plots_to_show = setdiff(plots_to_show, {'residual_history'});
end
if isempty(energy_hist)
    plots_to_show = setdiff(plots_to_show, {'energy_history'});
end
if isempty(phiv_refs) || isempty(phi_ref_titles)
    plots_to_show = setdiff(plots_to_show, {'phi_difference', 'phi_difference_close'});
end
if isempty(lapv)
    plots_to_show = setdiff(plots_to_show, {'laplacian_phi', 'laplacian_phi_close'});
end
if (nBodies <= 1) || isempty(Rss) || isempty(zbs)
    plots_to_show = setdiff(plots_to_show, ...
    {'hiramatsu_screening_stat', 'hiramatsu_screening_stat_close'});
end
if isempty(rrange_close) || isempty(zrange_close)
    plots_to_show = setdiff(plots_to_show, ...
    {'phi_difference_close', 'phi_close', 'initial_phi_close', 'residual_close',...
    'energy_density_close', 'log_energy_density_close',...
    'hiramatsu_screening_stat_close', 'laplacian_phi_close',...
    'newton_fraction_close'});
end


% Set up axs and fs (axes and figures)
if isempty(axs);
    axs = containers.Map;
    fs = containers.Map;
    for j = 1:length(plots_to_show)
        plot_str = plots_to_show{j};
        fs(plot_str) = figure;
        axs(plot_str) = axes('Parent', fs(plot_str));
    end
else
    plots_to_show = intersect(plots_to_show, axs.keys);
    closed_plots = {};
    fs = containers.Map;
    for j = 1:length(plots_to_show)
        plot_str = plots_to_show{j};
        if ~ishandle(axs(plot_str))
            closed_plots{end+1} = plot_str;
            fs(plot_str) = [];
        else
            this_object = axs(plot_str);
            while( ~isa(this_object, 'matlab.ui.Figure') && isprop(this_object, 'Parent') && (length(this_object.Parent) > 0) )
                this_object = this_object.Parent;
            end
            if isa(this_object, 'matlab.ui.Figure')
                fs(plot_str) = this_object;
            else
                fs(plot_str) = [];
            end
        end
    end
    plots_to_show = setdiff(plots_to_show, closed_plots);
end

% Set up some coordinate data
rv = nodedata.x';
zv = nodedata.y';

% For MATLAB to make a nice contour plots, we block out the body edges
eblockv = 0*rv + 1;
if (length(Rss) == nBodies) && (length(zbs) == nBodies)
    for bj = 1:nBodies
        eblockv( abs((rv.^2 + (zv-zbs(bj)).^2) - Rss(bj)^2) <= Rss(bj)^2*0.25) = NaN;
    end
end


% Construct all the plots
should_show = @(plot_name) ( length(intersect(plots_to_show, {plot_name}))>0);

if should_show('phi')
    try
    this_ax = axs('phi');
    contourf_interp(this_ax, meshes, idgrids, nodedata, phiv, 'XResolution', 100, 'YResolution', 100, 'AxesEqual', true);
    %if ~isempty(phiv_refs)
        %phiv_ref = phig_refs{1};
        %phi_ref_title = phi_ref_titles{1};
        %hold(this_ax, 'on')
        %contourf(this_ax, meshes, idgrids, nodedata, phiv_ref, 'XResolution', 100, 'YResolution', 100);
        %hold(this_ax, 'off')
        %legend(this_ax, {'phi', phi_ref_title});
    %end
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, ['phi' iter_str]);
    colorbar('peer', this_ax);
    drawnow
    catch ME
        warning('Error plotting phi');
    end
end

if should_show('phi_close')
    try
    this_ax = axs('phi_close');
    contourf_interp(this_ax, meshes, idgrids, nodedata, phiv, 'XResolution', 100, 'YResolution', 100, 'XRange', rrange_close, 'YRange', zrange_close, 'AxesEqual', true);
    %if ~isempty(phiv_refs)
        %phiv_ref = phig_refs{1};
        %phi_ref_title = phi_ref_titles{1};
        %hold(this_ax, 'on')
        %contourf(this_ax, meshes, idgrids, nodedata, phiv_ref, 'XResolution', 100, 'YResolution', 100, 'XRange', rrange_close, 'YRange', zrange_close);
        %hold(this_ax, 'off')
        %legend(this_ax, {'phi_close', phi_ref_title});
    %end
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, ['phi' iter_str]);
    colorbar('peer', this_ax);
    drawnow
    catch ME
        warning('Error plotting phi_close');
    end
end

if should_show('phi_difference')
    try
    if ~isempty(phiv_refs)
        % Plot (phi - phi_ref)/phi_ref
        phiv_ref = phiv_refs{1};
        phi_ref_title = phi_ref_titles{1};
        this_ax = axs('phi_difference');
        phidiffv = (phiv - phiv_ref)./phiv_ref;
        contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*phidiffv, 'AxesEqual', true);
        colorbar('peer',this_ax)
        xlabel(this_ax, 'r');
        ylabel(this_ax, 'z');
        title(this_ax, sprintf('(phi - %s)/%s%s', phi_ref_title, phi_ref_title, iter_str));
        drawnow
    end
    catch ME
        warning('Error plotting phi_difference');
    end
end

if should_show('phi_difference_close')
    try
    if ~isempty(phiv_refs)
        % Plot (phi - phi_ref)/phi_ref
        phiv_ref = phiv_refs{1};
        phi_ref_title = phi_ref_titles{1};
        this_ax = axs('phi_difference_close');
        phidiffv = (phiv - phiv_ref)./phiv_ref;
        contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*phidiffv, 'AxesEqual', true, 'XRange', rrange_close, 'YRange', zrange_close);
        colorbar('peer',this_ax)
        xlabel(this_ax, 'r');
        ylabel(this_ax, 'z');
        title(this_ax, sprintf('(phi - %s)/%s%s', phi_ref_title, phi_ref_title, iter_str));
        drawnow
    end
    catch ME
        warning('Error plotting phi_difference_close');
    end
end

if should_show('initial_phi')
    try
    this_ax = axs('initial_phi');
    contourf_interp(this_ax, meshes, idgrids, nodedata, phi_initv, 'XResolution', 100, 'YResolution', 100, 'AxesEqual', true);
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, 'phi (initial condition)');
    colorbar('peer', this_ax);
    drawnow
    catch ME
        warning('Error plotting initial_phi');
    end
end

if should_show('initial_phi_close')
    try
    this_ax = axs('initial_phi_close');
    contourf_interp(this_ax, meshes, idgrids, nodedata, phi_initv, 'XResolution', 100, 'YResolution', 100, 'XRange', rrange_close, 'YRange', zrange_close);
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, 'phi (initial condition)');
    colorbar('peer', this_ax);
    drawnow
    catch ME
        warning('Error plotting initial_phi_close');
    end
end

if should_show('laplacian_phi')
    try
    this_ax = axs('laplacian_phi');
    contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*(lapv*phiv), 'AxesEqual', true);
    colorbar('peer',this_ax)
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, sprintf('Laplacian phi%s', iter_str));
    drawnow
    catch ME
        warning('Error plotting laplacian_phi');
    end
end

if should_show('laplacian_phi_close')
    try
    this_ax = axs('laplacian_phi_close');
    contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*(lapv*phiv), 'AxesEqual', true, 'XRange', rrange_close, 'YRange', zrange_close);
    colorbar('peer',this_ax)
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, sprintf('Laplacian phi%s', iter_str));
    drawnow
    catch ME
        warning('Error plotting laplacian_phi_close');
    end
end

if should_show('residual')
    try
    this_ax = axs('residual');
    contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*Resv, 'AxesEqual', true);
    colorbar('peer',this_ax)
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, sprintf('Residual%s', iter_str));
    drawnow
    catch ME
        warning('Error plotting residual');
    end
end

if should_show('residual_close')
    try
    this_ax = axs('residual_close');
    contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*Resv, 'AxesEqual', true, 'XRange', rrange_close, 'YRange', zrange_close);
    colorbar('peer',this_ax)
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, sprintf('Residual%s', iter_str));
    drawnow
    catch ME
        warning('Error plotting residual_close');
    end
end

if should_show('energy_density')
    try
    this_ax = axs('energy_density');
    contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*endensv, 'AxesEqual', true);
    colorbar('peer',this_ax)
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, sprintf('Energy Density%s', iter_str));
    drawnow
    catch ME
        warning('Error plotting energy_density');
    end
end

if should_show('energy_density_close')
    try
    this_ax = axs('energy_density_close');
    contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*endensv, 'AxesEqual', true, 'XRange', rrange_close, 'YRange', zrange_close);
    colorbar('peer',this_ax)
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, sprintf('Energy Density%s', iter_str));
    drawnow
    catch ME
        warning('Error plotting energy_density_close');
    end
end

if should_show('log_energy_density')
    try
    this_ax = axs('log_energy_density');
    contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*log(abs(endensv)), 'AxesEqual', true);
    colorbar('peer',this_ax)
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, sprintf('log(abs(energy density))%s', iter_str));
    drawnow
    catch ME
        warning('Error plotting log_energy_density');
    end
end

if should_show('log_energy_density_close')
    try
    this_ax = axs('log_energy_density_close');
    contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*log(abs(endensv)), 'AxesEqual', true, 'XRange', rrange_close, 'YRange', zrange_close);
    colorbar('peer',this_ax)
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, sprintf('log(abs(energy density))%s', iter_str));
    drawnow
    catch ME
        warning('Error plotting log_energy_density_close');
    end
end

if should_show('hiramatsu_screening_stat')
    try
    if (nBodies > 0)
        this_ax = axs('hiramatsu_screening_stat');
        phi_hm_screen_stat = (lapv*phiv - lap_phi_sumv)./(lap_phi_onevs{1});
        contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*phi_hm_screen_stat, 'Levels', 100, 'LineColor','none');
        hold(this_ax, 'on');
        contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*phi_hm_screen_stat, 'LevelList', [-1.6:0.2:0], 'LineColor','k', 'Fill', false);
        hold(this_ax, 'off');
        caxis(this_ax, [-1.6 0]);
        colormap(this_ax, rjet);
        colorbar('peer',this_ax)
        xlabel(this_ax, 'r');
        ylabel(this_ax, 'z');
        title(this_ax, sprintf('Hiramatsu et al. screening statistic%s', iter_str));
        drawnow
    end
    catch ME
        warning('Error plotting hiramatsu_screening_stat');
    end
end

if should_show('hiramatsu_screening_stat_close')
    try
    if (nBodies > 0)
        this_ax = axs('hiramatsu_screening_stat_close');
        phi_hm_screen_stat = (lapv*phiv - lap_phi_sumv)./(lap_phi_onevs{1});
        contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*phi_hm_screen_stat, 'Levels', 100, 'LineColor','none', 'XRange', rrange_close, 'YRange', zrange_close);
        hold(this_ax, 'on');
        contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*phi_hm_screen_stat, 'LevelList', [-1.6:0.2:0], 'LineColor','k', 'Fill', false, 'XRange', rrange_close, 'YRange', zrange_close);
        hold(this_ax, 'off');
        caxis(this_ax, [-1.6 0]);
        colormap(this_ax, rjet);
        colorbar('peer',this_ax)
        xlabel(this_ax, 'r');
        ylabel(this_ax, 'z');
        title(this_ax, sprintf('Hiramatsu et al. screening statistic%s', iter_str));
        drawnow
    end
    catch ME
        warning('Error plotting hiramatsu_screening_stat_close');
    end
end

if should_show('newton_fraction')
    try
    this_ax = axs('newton_fraction');
    gradphiv = sqrt((drv*phiv).^2 + (dzv*phiv).^2);
    gradpsiv = sqrt((drv*newton_psiv).^2 + (dzv*newton_psiv).^2);
    contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*(gradphiv./gradpsiv), 'AxesEqual', true);
    colorbar('peer',this_ax)
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, sprintf('|grad phi|/|grad Newt. pot.| %s', iter_str));
    drawnow
    catch ME
        warning('Error plotting newton_fraction');
    end
end

if should_show('newton_fraction_close')
    try
    this_ax = axs('newton_fraction_close');
    gradphiv = sqrt((drv*phiv).^2 + (dzv*phiv).^2);
    gradpsiv = sqrt((drv*newton_psiv).^2 + (dzv*newton_psiv).^2);
    contourf_interp(this_ax, meshes, idgrids, nodedata, eblockv.*(gradphiv./gradpsiv), 'AxesEqual', true, 'XRange', rrange_close, 'YRange', zrange_close);
    colorbar('peer',this_ax)
    xlabel(this_ax, 'r');
    ylabel(this_ax, 'z');
    title(this_ax, sprintf('|grad phi|/|grad Newt. pot.| %s', iter_str));
    drawnow
    catch ME
        warning('Error plotting newton_fraction_close');
    end
end

if should_show('residual_history')
    try
    this_ax = axs('residual_history');
    if length(intres_hist) > 0
        intres = intres_hist(end);
    else
        intres = NaN;
    end
    semilogy(this_ax, [0:this_iter], intres_hist);
    xlabel(this_ax, 'Iteration');
    ylabel(this_ax, '||Res||^2');
    title(this_ax, sprintf('Integrated objective: %g%s', intres, iter_str));
    xlim(this_ax, [0 this_iter+1]);
    drawnow
    catch ME
        warning('Error plotting residual_history');
    end
end

if should_show('energy_history')
    try
    this_ax = axs('energy_history');
    if length(energy_hist) > 0
        energy = energy_hist(end);
    else
        energy = NaN;
    end
    semilogy(this_ax, [0:this_iter], energy_hist);
    xlabel(this_ax, 'Iteration');
    ylabel(this_ax, 'E');
    title(this_ax, sprintf('Total energy: %g%s', energy, iter_str));
    xlim(this_ax, [0 this_iter+1]);
    drawnow
    catch ME
        warning('Error plotting energy_history');
    end
end




end
