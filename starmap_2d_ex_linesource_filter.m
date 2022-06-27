function starmap_2d_ex_linesource_filter
%STARMAP_2D_EX_LINESOURCE_FILTER
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: Line source benchmark test case, as defined in
%   [Ganapol, JCP 210 (2005) 386-399].
%   This test case is adapted to using a narrow Gaussian as
%   initial condition. This file is an adaptation of the
%   test case STARMAP_EX_LINESOURCE for the use of the FPN
%   approximation. Feel free to modify the code to PN by
%   removing the filtering, to see the increase in quality
%   due to the filtering.
%
%   Remarks:
%   Since the numerical scheme does not possess slope
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

%   For license, see files LICENSE.txt or starmap_solver.m, as published on
%   https://github.com/starmap-project/starmap


%========================================================================
% Problem Parameters
%========================================================================
prob = struct(...
'name','Line Source Test',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',13,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'ic',@initial,... % initial conditions (defined below)
'ax',[-1 1 -1 1]*0.6,... % coordinates of computational domain
'n',[1 1]*150,... % numbers of grid cells in each coordinate direction
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,0.5,21),... % output times
'output',@output,... % problem-specific output routine (defined below)
'filter','exp',... % type of filter (can be 'exp' or 'spherical')
'f_position','substep'... % position of filter ('substep','full', or '')
);
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
switch prob.filter                             % Define filter function.
    case 'exp', prob.filterfunction = @exp_filter;
    case 'spherical', prob.filterfunction = @spherical_spline_filter;
end
par = starmap_init(prob); % Configure data structures for starmap solver.
starmap_solver(par)                                         % Run solver.
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

function f = sigma_s0(x,y)
% Isotropic scattering coefficient.
f = 1;

function f = exp_filter(par,dt,k)
% Filter function, multiplies each decay coefficient by an exponential.
p = 2;                                             % Order of the filter.
f_eff = 20;                                   % Filter effective opacity.
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
% Plotting routine, showing the 2D results and a 1D cross-section.
t = par.t_plot(step);                                   % Time of output.
[r0,phi0] = line_source_solution(t);             % Compute true solution.
phi_max = max(phi0(1)+eps,max(max(U))/4);           % Axis/color scaling.
clf, subplot(1,3,1:2), imagesc(x,y,U')        % 2D plot of approximation.
axis xy equal tight, caxis([-1 1]*phi_max*1.6)
title(sprintf('%s with FP%d at t = %0.2f',par.name,par.n_mom,t))
colormap jet(255); colorbar
subplot(1,3,3)
r = linspace(0,max(x),100); phi = interp2(x,y,U,r,0); % 1D cross section.
plot(r0,phi0,'-','linewidth',2,'color',[.6 .6 1])   % Plot true solution.
hold on, plot(r,phi,'r-','linewidth',1), hold off   % Plot approximation.
axis([0 max(x) [-0.75 2.5]*phi_max])
legend('true solution',sprintf('FP%d model',par.n_mom),'Location','SouthEast')
title('Cut at x>0, y=0')
drawnow

%========================================================================
% Reference solution of line source problem
%========================================================================
function [rho,f] = line_source_solution(t)
% Evaluates the true solution of the line source test case,
% using Ganapol's quadrature formulas.
rho = [(1-linspace(1,0,25).^2)*t 1];
f = phi_l(rho,t);

function f = phi_l(rho,t)
eta = rho/t;
ind = find(eta<1);
f = rho*0; phi_l0 = f;
for k = ind                                           % Compute integral.
    f(k) = quad(@(w) phi_pt(t*sqrt(eta(k).^2+w.^2),t),...
    0,sqrt(1-eta(k).^2),1e-5);
end
phi_l0(ind) = exp(-t)/(2*pi*t^2)./sqrt(1-eta(ind).^2);
f = phi_l0+(2*t)*f;

function f = phi_pt(r,t)
r = max(r,1e-10);           % Move zero value into small positive regime.
eta = r/t;
ind = find(eta<1);
g = r*0;
for k = ind                                           % Compute integral.
    g(k) = quad(@(u) sec(u/2).^2.*real(...
        (eta(k)+1i*tan(u/2)).*...
        xi(eta(k),u).^3.*exp(t/2*(1-eta(k).^2).*xi(eta(k),u)) ...
        ),0,pi,1e-2);
end
f = 1/(2*pi)*exp(-t)./(4*pi*r.*t.^2)*(t/2)^2.*(1-eta.^2).*g;   % Assemble
f(ind) = f(ind)+phi_pt1(r(ind),t);                        % final result.

function f = xi(eta,u)
q = (1+eta)./(1-eta);
q = min(max(q,-1e80),1e80);              % Remove singularity at eta = 1.
f = (log(q)+1i*u)./(eta+1i*tan(u/2));

function y = phi_pt1(r,t)
eta = r/t;
q = (1+eta)./(1-eta);
q = min(max(q,-1e80),1e80);              % Remove singularity at eta = 1.
y = exp(-t)./(4*pi*r*t^2)*t.*log(q);
