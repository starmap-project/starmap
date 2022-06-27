function starmap_2d_ex_controlrod
%STARMAP_2D_EX_CONROLROD
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: Control rod, as defined in Section 4.2 in
%   [Olbrant, Larsen, Frank, Seibold, JCP 238 (2013) 315-336].
%   This test case is a model for the temporal evolution of
%   the neutron density in a nuclear reactor, when an
%   absorbing rod is moved into and out of the computational
%   domain. The particularity of this test is its
%   time-dependent material parameters. This example
%   demonstrates that the StaRMAP code can deal with
%   time-dependent material parameters in an accurate fashion.
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
'name','Control Rod Test',... % name of example
'closure','SP',... % type of closure (can be 'P' or 'SP')
'n_mom',3,... % order of moment approximation
'sigma_a',@effective_absorb,... % effective absorption coeff. (below)
'sigma_s0',@effective_scatter,... % effective scattering coeff. (below)
'ic',@initial,... % initial conditions (defined below)
'ax',[-1 1 -1 1],... % coordinates of computational domain
'n',[251 251],... % numbers of grid cells in each coordinate direction
'bc',[0 0],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,0.9,46),... % output times
'output',@output... % problem-specific output routine (defined below)
);
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
par = starmap_init(prob); % Configure data structures for starmap solver.
starmap_solver(par)                                         % Run solver.
%========================================================================
% Problem Specific Functions
%========================================================================
function f = initial(x,y)
% Initial conditions (no k, thus only for zeroth moment).
f = 1e6*ones(size(x));

function f = effective_absorb(x,y,t)
% Effective absorption coefficient.
f = rod_geometry(x,y)*time_dependence(t);

function f = effective_scatter(x,y)
% Effective scattering coefficient.
sigma_s = 1;                           % physical scattering coefficient.
sigma_f = 2;                                       % fission coefficient.
nu = 0.9;                     % mean number of fission neutrons released.
f = (sigma_s+nu*sigma_f)*ones(size(x));

function chi = rod_geometry(x,y)
% Characteristic function of control rod geometry.
chi = -0.5<=x&x<=0.2&-0.5<=y&y<=0 | -0.3<=x&x<=0&0<y&y<=0.6 | ...
    0<x&x<=0.5&0<y&y<0.3;

function s = time_dependence(t)
T = 0.6;                                      % duration of rod movement.
dT = 0.2;                  % time of rod in-motion and of rod out-motion.
smax = 10;
s = smax/dT*t.*(0<t&t<=dT)+smax*(dT<t&t<=T-dT)+...
    smax/dT*(T-t).*(T-dT<t&t<=T);

function output(par,x,y,U,step)
% Plotting routine.
clf, subplot(1,3,1:2)
imagesc(x,y,U'), axis xy equal tight, caxis([0 1e6])
title(sprintf('%s with %s%d at t = %0.2f',par.name,par.closure,...
    par.n_mom,par.t_plot(step)))
colormap jet(255); colorbar
xlabel('x'), ylabel('y')
subplot(1,3,3)
plot(x,interp2(x,y,U',x,y*0))     % Plot solution evaluated along x-axis.
axis([-1 1 0 1.01e6])
title('Cut at y=0'), xlabel('x')
drawnow
