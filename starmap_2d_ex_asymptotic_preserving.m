function starmap_2d_ex_asymptotic_preserving
%STARMAP_2D_EX_ASYMPTOTIC_PRESERVING
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometries.
%
%   Example: Demonstrates asymptotic preserving properties of
%   STARMAP_SOLVER. This test case has a narrow Gaussian as
%   initial condition. Comparisons are made between the use
%   of a hyperbolic timestep and a parabolic timestep in regimes
%   of large scattering coefficients. The final time and isotropic
%   scattering depend on a parameter epsilon. For each value of
%   epsilon, the solver is called with a hyperbolic timestep then a
%   parabolic timestep. Four figures are generated.
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
%   http://www.math.temple.edu/~seibold/research/starmap

%========================================================================
% Problem Parameters
%========================================================================
prob = struct(...
'name','Asymptotic Preserving Test',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',1,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'ic',@initial,... % initial conditions (defined below)
'ax',[-1 1 -1 1]*0.6,... % coordinates of computational domain
'n',[1 1]*150,... % numbers of grid cells in each coordinate direction
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
'output',@output... % problem-specific output routine (defined below)
);
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
eps = [1e-3,1e-4]; % scaling parameters
for ep=eps
    prob.t_plot = linspace(0,0.005,21)/ep; % output times
    prob.sigma_s0 = @(x,y) 1/ep; % Isotropic scattering coefficient.

    prob.timestep = 0; % use hyperbolic timestep
    figure
    par = starmap_init(prob);  % Configure data structures for starmap solver
    starmap_solver(par)                                     % Run solver.

    prob.timestep = 1; % use parabolic timestep
    figure
    par = starmap_init(prob); % Configure data structures for starmap solver
    starmap_solver(par)                                     % Run solver.
end

%========================================================================
% Problem Specific Functions
%========================================================================
function f = initial(x,y,k)
% Initial conditions (for (k-1)-st moment).
t = 3.2e-4;              % Pseudo-time for smoothing initial Dirac delta.
f = 1/(4*pi*t)*exp(-(x.^2+y.^2)/(4*t))*(k==1);

function f = sigma_a(x,y)
% Absorption coefficient.
f = 0;

function output(par,x,y,U,step)
% Plotting routine, showing the 2D results and a 1D cross-section.
t = par.t_plot(step);                                   % Time of output.
phi_max = max(max(max(U))/2);           % Axis/color scaling.
% Timestep used - hyperbolic or parabolic
if par.timestep >0, type = 'Parabolic'; else, type = 'Hyperbolic'; end
txt = strcat('Isotropic scattering = ',num2str(par.prob.sigma_s0(0,0)),...
    ',  ',type, ' timestepping');
clf, subplot(1,3,1:2), imagesc(x,y,U')        % 2D plot of approximation.
axis xy equal tight, caxis([-1 1]*phi_max*1.6)
title(sprintf('%s with %s%d at t = %0.2f',par.name,par.closure,par.n_mom,t))
colormap jet(255); colorbar
subplot(1,3,3)
r = linspace(0,max(x),100); phi = interp2(x,y,U,r,0); % 1D cross section.
hold on, plot(r,phi,'r-','linewidth',1), hold off   % Plot approximation.
axis([0 max(x) [-0.75 2.5]*phi_max])
legend(sprintf('%s%d model',par.closure,par.n_mom),'Location','SouthEast')
title('Cut at x>0, y=0')
annotation('textbox',[0,0,1, 1],'String',txt);
drawnow
