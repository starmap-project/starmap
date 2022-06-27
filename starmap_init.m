function par = starmap_init(prob)
%STARMAP_INIT
%   Converts user-defined variables in prob to format compatible with
%   STARMAP_SOLVER, a second order staggered grid finite
%   difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%   Sets default variables and plotting routines.
%   For 1D and 2D examples, variables are assumed to be defined to
%   match dimensionality.
%   A copy of user-defined prob is attached to output structure par.
%
%   Version 2.0
%   Copyright (c) 06/28/2022 Benjamin Seibold, Martin Frank, and
%                            Rujeko Chinomona
%   http://www.math.temple.edu/~seibold
%   https://www.scc.kit.edu/personen/martin.frank.php
%   https://rujekoc.github.io/
%
%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0).
%
%   StaRMAP project website:
%   https://github.com/starmap-project

%   For license, see files LICENSE.txt or starmap_solver.m, as published on
%   https://github.com/starmap-project/starmap

%========================================================================
% Set struct
%========================================================================
par = prob;                                           % Copy prob to par.
par.prob = prob;                                 % Store original struct.

%========================================================================
% Set dimension of problem
%========================================================================
% Detect dimension of problem
if isfield(prob,'ax'), par.dim = length(prob.ax)/2;
elseif isfield(prob,'n'), par.dim = length(prob.n);
elseif isfield(prob,'bc'), par.dim = length(prob.bc);
else, error('At least one of ax, n, or bc should be defined');
end
% Check if variables are compatible
if isfield(prob,'ax')
    if isfield(prob,'n')
        if abs(length(prob.ax)/2 - length(prob.n))>0
            error('Incompatible problem variables ax and n.')
        end
    end
    if isfield(prob,'bc')
        if abs(length(prob.ax)/2 - length(prob.bc))>0
            error('Incompatible problem variables ax and bc.')
        end
    end
elseif isfield(prob,'n')
    if isfield(prob,'bc')
        if abs(length(prob.n) - length(prob.bc))>0
            error('Incompatible problem variables n and bc.')
        end
    end
end
%========================================================================
% Set Moment Matrices and system components
%========================================================================
switch par.closure                      % Define closure matrix function.
    case  'P', closurefun = 'starmap_closure_pn';
    case 'SP', closurefun = 'starmap_closure_spn';
end
% Compute moment matrices.
[par.mom_order,par.Mx,par.My,par.Mz] = feval(closurefun,par.n_mom,par.dim);
n_sys = length(par.Mz);                    % Number of system components.

%========================================================================
% Set data structures for starmap_solver
%========================================================================
switch par.dim
    case 1
        if isfield(prob,'ax'), par.ax = [prob.ax prob.ax prob.ax]; end %Axes
        if isfield(prob,'n'), par.n = [1 1 prob.n]; end % No. of grid cells
        if isfield(prob,'bc'), par.bc = prob.bc*[1 1 1]; end  % Boundary.
        if isfield(prob,'ic')                       % Initial conditions.
           if nargin(prob.ic)==1
               par.ic = @(~,~,x) prob.ic(x);
           elseif nargin(prob.ic)==2
               par.ic = @(~,~,x,k) prob.ic(x,k);
           else
               error('ic has incompatible number of arguments.')
           end
        end
        if isfield(prob,'sigma_a')              % Absorption coefficient.
           if nargin(prob.sigma_a)==1
               par.sigma_a = @(~,~,x) prob.sigma_a(x);
           elseif nargin(prob.sigma_a)==2               % Time dependent.
               par.sigma_a = @(~,~,x,t) prob.sigma_a(x,t);
           else
               error('sigma_a has incompatible number of arguments.')
           end
        end
        if isfield(prob,'sigma_s0')   % Isotropic scattering coefficient.
           if nargin(prob.sigma_s0)==1
               par.sigma_s0 = @(~,~,x) prob.sigma_s0(x);
           elseif nargin(prob.sigma_s0)==2              % Time dependent.
               par.sigma_s0 = @(~,~,x,t) prob.sigma_s0(x,t);
           else
               error('sigma_s0 has incompatible number of arguments.')
           end
        end
        if isfield(prob,'sigma_sm') % Anisotropic scattering coefficient.
           if nargin(prob.sigma_sm)==2
               par.sigma_sm = @(~,~,m,x) prob.sigma_sm(x,m);
           elseif nargin(prob.sigma_sm) == 3            % Time dependent.
               par.sigma_sm = @(~,~,x,m,t) prob.sigma_sm(x,m,t);
           else
               error('sigma_sm has incompatible number of arguments.')
           end
        end
        if isfield(prob,'source')                          % Source term.
           if nargin(prob.source)==2
               par.source = @(~,~,x,t) prob.source(x,t);
           elseif nargin(prob.source)==3             % Moments specified.
               par.source = @(~,~,x,t,k) prob.source(x,t,k);
           elseif nargin(prob.source) == 1
               par.source = @(~,~,x) prob.source(x);
           else
               error('source has incompatible number of arguments.')
           end
        end
        if isfield(prob,'filterfunction')              % Filter function.
           if nargin(prob.filterfunction)==4
               par.filterfunction = @(par,dt,k,~,~,x) ...
                   prob.filterfunction(par,dt,k,x);
           elseif nargin(prob.filterfunction) == 5      % Time dependent.
               par.filterfunction = @(par,dt,k,~,~,x,t) ...
                   prob.filterfunction(par,dt,k,x,t);
           elseif nargin(prob.filterfunction) == 3
               par.filterfunction = prob.filterfunction;
           else
               error('filterfunction has incompatible number of arguments.')
           end
        end
        if ~isfield(par,'output'),   par.output = @default_output_1d; end
    case 2
        if isfield(prob,'ax'), par.ax = [prob.ax(1:2) prob.ax]; end %Axes
        if isfield(prob,'n'), par.n = [1 prob.n]; end % Number of grid cells
        if isfield(prob,'bc'), par.bc =[prob.bc(1) prob.bc]; end % Boundary.
        if isfield(prob,'ic')                        % Initial conditions.
           if nargin(prob.ic)==2
               par.ic = @(~,x,y) prob.ic(x,y);
           elseif nargin(prob.ic)==3                 % Moments specified.
               par.ic = @(~,x,y,k) prob.ic(x,y,k);
           else
               error('ic has incompatible number of arguments.')
           end
        end
        if isfield(prob,'sigma_a')              % Absorption coefficient.
           if nargin(prob.sigma_a)==2
               par.sigma_a = @(~,x,y) prob.sigma_a(x,y);
           elseif nargin(prob.sigma_a)==3               % Time dependent.
               par.sigma_a = @(~,x,y,t) prob.sigma_a(x,y,t);
           else
               error('sigma_a has incompatible number of arguments.')
           end
        end
        if isfield(prob,'sigma_s0')   % Isotropic scattering coefficient.
           if nargin(prob.sigma_s0)==2
               par.sigma_s0 = @(~,x,y) prob.sigma_s0(x,y);
           elseif nargin(prob.sigma_s0)==3              % Time dependent.
               par.sigma_s0 = @(~,x,y,t) prob.sigma_s0(x,y,t);
           else
               error('sigma_s0 has incompatible number of arguments.')
           end
        end
        if isfield(prob,'sigma_sm') % Anisotropic scattering coefficient.
           if nargin(prob.sigma_sm)==3
               par.sigma_sm = @(~,x,y,m) prob.sigma_sm(x,y,m);
           elseif nargin(prob.sigma_sm) == 4            % Time dependent.
               par.sigma_sm = @(~,x,y,m,t) prob.sigma_sm(x,y,m,t);
           else
               error('sigma_sm has incompatible number of arguments.')
           end
        end
        if isfield(prob,'source')                          % Source term.
           if nargin(prob.source)==3
               par.source = @(~,x,y,t) prob.source(x,y,t);
           elseif nargin(prob.source)==4             % Moments specified.
               par.source = @(~,x,y,t,k) prob.source(x,y,t,k);
           elseif nargin(prob.source) == 2
               par.source = @(~,x,y) prob.source(x,y);
           else
               error('source has incompatible number of arguments.')
           end
        end
        if isfield(prob,'filterfunction')              % Filter function.
           if nargin(prob.filterfunction)==5
               par.filterfunction = @(par,dt,k,~,x,y) ...
                   prob.filterfunction(par,dt,k,x,y);
           elseif nargin(prob.filterfunction) == 6      % Time dependent.
               par.filterfunction = @(par,dt,k,~,x,y,t) ...
                   prob.filterfunction(par,dt,k,x,y,t);
           elseif nargin(prob.filterfunction) == 3
               par.filterfunction = prob.filterfunction;
           else
               error('filterfunction has incompatible number of arguments.')
           end
        end
        if ~isfield(par,'output'),   par.output = @default_output_2d; end
    case 3
        if ~isfield(par,'output'),   par.output = @default_output_3d; end
end

%========================================================================
% Set other defaults
%========================================================================
if ~isfield(par,'name'),     par.name = 'Output'; end
if ~isfield(par,'ic'),       par.ic = @zero; end
if ~isfield(par,'sigma_a'),  par.sigma_a = @zero; end
if ~isfield(par,'sigma_s0'), par.sigma_s0 = @zero; end
if ~isfield(par,'source'),  par.source = @zero; end
if ~isfield(par,'source_ind') % If no source moment vector specified, use
    par.source_ind = 1:((nargin(par.source)>=0)+...   % zeroth moment for
        (nargin(par.source)>4)*(n_sys-1));   % isotropic, and all moments
end                                             % for anisotropic source.
if ~isfield(par,'ax'),       par.ax = [0 1 0 1 0 1]; end % Problem domain.
if ~isfield(par,'n'),        par.n = [100 100 100]; end %No. of grid cells.
if ~isfield(par,'timestep'), par.timestep = 0; end %0=hyp; 1=par; 2=par0.
if ~isfield(par,'cn'),       par.cn = 0.99; end     % Default CFL number.
if ~isfield(par,'bc'),       par.bc = [0 0 0]; end  % 0=periodic; 1=extrap.
if ~isfield(par,'t_plot'),   par.t_plot = 0:.1:1; end % Plotting times
if ~isfield(par,'tfinal'),   par.tfinal = par.t_plot(end); end
if ~isfield(par,'mom_output'), par.mom_output = 1; end % Moments plotted.
if isfield(par,'filterfunction'), par.closure_mod = 'F'; % Closure modifier
else, par.closure_mod = '';end

%========================================================================
% Utility functions
%========================================================================
function f = zero(varargin)
% Zero function.
f = zeros(size(varargin{1}));

function default_output_1d(par,x,U,step)
% Default plotting routine.
clf, plot(x,U), xlim(par.ax(1:2)), xlabel('x');
title(sprintf('%s with %s%d at t = %0.2f',par.name,par.closure,...
    par.n_mom,par.t_plot(step)))
drawnow


function default_output_2d(par,x,y,U,step)
% Default plotting routine.
clf, imagesc(x,y,U'), axis xy equal tight
title(sprintf('%s with %s%d at t = %0.2f',par.name,par.closure,...
    par.n_mom,par.t_plot(step)))
drawnow

function default_output_3d(par,x,y,z,U,step)
% Default plotting routine.
clf, imagesc(x,y,U(:,:,ceil(length(z)/2))'), axis xy equal tight
title(sprintf('%s with %s%d at t = %0.2f',par.name,par.closure,...
    par.n_mom,par.t_plot(step)))
drawnow





