function starmap_2d_ex_lattice_filter
%STARMAP_2DEX_LATTICE_FILTER
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: Lattice/Checkerboard geometry, as defined in
%   [Brunner, Holloway, JCP 210 (2005) 386-399].
%   Computation up to time t=3.2. This file is an adaptation
%   of the test case STARMAP_2D_EX_LATTICE for the use of the
%   FPN approximation.
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

%   For license, see file starmap_solver.m, as published on
%   https://github.com/starmap-project/starmap

%========================================================================
% Problem Parameters
%========================================================================
prob = struct(...
'name','Lattice Test',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',5,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'source',@source,... % source term (defined below)
'ax',[0 7 0 7],... % coordinates of computational domain
'n',[250 250],... % numbers of grid cells in each coordinate direction
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,3.2,51),... % output times
'output',@output,... % problem-specific output routine (defined below)
'filter','exp',... % type of filter (can be 'exp' or 'spherical')
'f_position','substep'... % position of filter ('substep','full', or '')
);
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
switch prob.filter                     % Filter function (defined below).
    case 'exp', prob.filterfunction = @exp_filter;
    case 'spherical', prob.filterfunction = @spherical_spline_filter;
end
par = starmap_init(prob); % Configure data structures for starmap solver.
starmap_solver(par)                                         % Run solver.
%========================================================================
% Problem Specific Functions
%========================================================================
function f = sigma_a(x,y)
% Absorption coefficient.
cx = ceil(x); cy = ceil(y);
g = (ceil((cx+cy)/2)*2==(cx+cy)).*(1<cx&cx<7&1<cy&cy-2*abs(cx-4)<4);
f = (1-g)*0+g*10;

function f = sigma_s0(x,y)
% Isotropic scattering coefficient.
cx = ceil(x); cy = ceil(y);
g = (ceil((cx+cy)/2)*2==(cx+cy)).*(1<cx&cx<7&1<cy&cy-2*abs(cx-4)<4);
f = (1-g)*1+g*0;

function f = source(x,y)
% Radiation source (only for zeroth moment).
f = 3<x&x<4&3<y&y<4;

function f = exp_filter(par,dt,k)
% Filter function, multiplies each decay coefficient by an exponential.
p = 2;                                             % Order of the filter.
f_eff = 5;                                    % Filter effective opacity.
f_b = -f_eff/log(exp(log(eps)*(par.n_mom/(par.n_mom+1))^p));
f = exp(log(eps)*(k/(par.n_mom+1)).^p).^(dt*f_b);

function f = spherical_spline_filter(par,dt,k,x,y)
% Filter function, multiplies each decay coefficient by a spherical
% spline filter term.
L = 1;                                     % Characteristic length scale.
fw = dt./max((par.ax([2 4])-par.ax([1 3]))./par.n)/2;
f_alpha = fw./par.n_mom.^2./...
    ((sigma_a(x,y)+sigma_s0(x,y)).*L+par.n_mom).^2;
f = 1./(1+f_alpha*k^2.*(k+1).^2);

function output(par,x,y,U,step)
% Plotting routine.
cax = [-7 0];                    % Colormap range used for log10 scaling.
vneg = cax(1)-diff(cax)/254;        % Value assigned where U is negative.
V = log10(max(U,1e-20));% Cap U s.t. U>0 and use logarithmic color scale.
Vcm = max(V,cax(1));                      % Cap colormap plot from below.
Vcm(U<0) = vneg;              % Assign special value where U is negative.
clf, subplot(1,3,1:2)
imagesc(x,y,Vcm'), axis xy equal tight, caxis([vneg cax(2)])
title(sprintf('%s with FP%d at t = %0.2f',par.name,par.n_mom,...
    par.t_plot(step)))
cm = jet(256); cm(1,:) = [1 1 1]*.5; % Change lowest color entry to gray.
colormap(cm), colorbar('ylim',cax)
xlabel('x'), ylabel('y')
subplot(1,3,3)
plot(y,interp2(x,y,V',y*0+3.5,y))   % Evaluate solution along line x=3.5.
axis([par.ax(1:2) cax+[-1 1]*.1])
title('Cut at x=3.5'), xlabel('y')
drawnow
