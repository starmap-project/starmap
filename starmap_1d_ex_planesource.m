function starmap_1d_ex_planesource
%STARMAP_1D_EX_PLANESOURCE
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: plane source, represented by a narrow Gaussian as
%   initial condition.
%
%   Remarks:
%   1) Since the numerical scheme does not possess slope
%   limiters, the initial Gaussian must not be chosen too
%   narrow (or else oscillations will arise).
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
'name','Plane Source Test',... % name of example
'n_mom',7,... % order of moment approximation
'closure','SP',... % type of closure
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'ic',@initial,... % initial conditions (defined below)
'ax',[-1 1]*0.6,... % coordinates of computational domain
'n',1000,... % numbers of grid cells in each coordinate direction
'bc',0,... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,0.4,21),... % output times
'output',@output... % problem-specific output routine (defined below)
);
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
par = starmap_init(prob); % Configure data structures for starmap solver.
starmap_solver(par);                                        % Run solver.

%========================================================================
% Problem Specific Functions
%========================================================================
function f = initial(x,k)
% Initial conditions (for (k-1)-st moment).
t = 3.2e-3;              % Pseudo-time for smoothing initial Dirac delta.
f = 1/(4*pi*t)*exp(-(x.^2)/(4*t))*(k==1);

function f = sigma_a(x)
% Absorption coefficient.
f = 0;

function f = sigma_s0(x)
% Isotropic scattering coefficient.
f = 1;

function output(par,x,U,step)
% Plotting routine.
plot(x,U)
xlim(par.ax(1:2))
xlabel('x')
title(sprintf('%s with %s%d at t = %0.2f',par.name,par.closure,...
    par.n_mom,par.t_plot(step)))
ylim([0,25])
drawnow
