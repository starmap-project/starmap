function starmap_1d_ex_gauss
%STARMAP_1D_EX_GAUSS
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: plane source, represented by a narrow Gaussian as
%   initial condition from the paper "Optimal prediction for moment models:
%   crescendo diffusion and reordered equations"
%   by B. Seibold & M. Frank https://doi.org/10.1007/s00161-009-0111-7.
%   Uses default plotting routine.
%
%   Version 2.0
%   Copyright (c) 06/28/2022 Benjamin Seibold, Martin Frank, and
%                            Rujeko Chinomona
%   http://www.math.temple.edu/~seibold
%   https://www.scc.kit.edu/personen/martin.frank.php
%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0)
%
%   StaRMAP project website:
%   https://github.com/starmap-project

%   For license, see file starmap_solver.m, as published on
%   https://github.com/starmap-project/starmap

%========================================================================
% Problem Parameters
%========================================================================
prob = struct(...
'name','1D Gaussian',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',20,... % order of moment approximation
'n',1000,... % number of grid cells
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'ic',@initial,... % initial conditions (defined below)
'ax',[0 1],... % coordinates of computational domain
'bc',0,... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,0.4,21)... % use default plotting routine
);
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
par = starmap_init(prob);     % Configure data structures for starmap solver
starmap_solver(par);          % Run solver

%-----------------------------------------------------------------------%


%========================================================================
% Problem Specific Functions
%========================================================================
function f = initial(x,k)
% first moment initial condition nonzero, all others zero
f = exp(-500*(x-.5).^2)*(k==1);

function f = sigma_a(x)
% Absorption coefficient.
f = 1.5;

function f = sigma_s0(x)
% Isotropic scattering coefficient.
f = 1.5;
