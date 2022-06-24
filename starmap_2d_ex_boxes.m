function starmap_2d_ex_boxes
%STARMAP_2D_EX_BOXES
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: Numerical demonstration of SPN-PN equivalence,
%   as defined in
%   [McClarren, Trans. Theory Stat. Phys. 39 (2011) 73-109].
%   For the given moment order (specified below), the example
%   is computed first with PN, and then with SPN, and the
%   difference between the solutions is evaluated. This test
%   demonstrates (a) that under the given assumptions, StaRMAP
%   yields identical (up to machine accuracy) solutions for PN
%   and SPN; and (b) that the SPN systems can be solved much
%   faster than the PN systems.
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
'name','Boxes Test',... % name of example
'n_mom',9,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'source',@source,... % source term (defined below)
'ax',[0 5 0 5],... % coordinates of computational domain
'n',[200 200],... % numbers of grid cells in each coordinate direction
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,1,21),... % output times
'output',@output... % problem-specific output routine (defined below)
);
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
clf
prob.closure = 'P';  % Set to PN.
par = starmap_init(prob);     % Configure data structures for starmap solver
solP = starmap_solver(par);                      % Run solver with PN.
prob.closure = 'SP';    % set to SPN.
par = starmap_init(prob);     % Configure data structures for starmap solver
solSP = starmap_solver(par);                    % Run solver with SPN.

% Plot differences between PN and SPN
subplot(1,2,2)
imagesc(solP(1).x,solP(1).y,abs(solP(1).U-solSP(1).U)')
axis xy equal tight
colormap jet(255); colorbar
title(sprintf('Difference between P%d and SP%d',par.n_mom,par.n_mom))
%========================================================================
% Problem Specific Functions
%========================================================================
function f = sigma_a(x,y)
% Absorption coefficient.
f = 0.9;

function f = sigma_s0(x,y)
% Isotropic scattering coefficient.
f = 0.1;

function f = source(x,y,t)
% Radiation source (only for zeroth moment).
f = (2+sin(4*pi*t)*exp(-t/3)) * ...
    (1.75<x&x<=2.25&1.75<y&y<=2.25 | 2.75<x&x<=3.25&1.5<y&y<=2.50 | ...
     1.75<x&x<=2.25&2.75<y&y<=3.25 | 3.50<x&x<=4.25&3.5<y&y<=3.75);

function output(par,x,y,U,step)
% Plotting routine, specific to this test case.
subplot(1,2,2-strcmp(par.closure(1),'P'))
imagesc(x,y,U'), axis xy equal tight, caxis([0 0.8])
title(sprintf('%s with %s%d at t = %0.2f',par.name,par.closure,...
    par.n_mom,par.t_plot(step)))
colormap jet(255); colorbar, drawnow
