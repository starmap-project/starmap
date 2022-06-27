function starmap_3d_ex_pointsource
%STARMAP_3D_EX_POINTSOURCE
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: Point Source benchmark test. 
%   Via displaying the numerical solution along different rays, 
%   this test case showcases the spatial approximation error to 
%   rotational symmetry of the numerical method. One can see how 
%   choosing a finer spatial grid will reduce that approximation error.
%
%   Remarks:
%   1) Since the numerical scheme does not possess slope
%   limiters, the initial Gaussian must not be chosen too
%   narrow (or else oscillations will arise).
%   2) For this problem, the SPN approximation is equivalent
%   to the PN approximation, see
%   [McClarren, Trans. Theory Stat. Phys. 39 (2011) 73-109],
%   [Olbrant, Larsen, Frank, Seibold, JCP 238 (2013) 315-336].
%   Since SPN contains much fewer moments than PN, it is
%   natural to use here.
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
'name','Point Source Test',... % name of example
'closure','SP',... % type of closure (can be 'P' or 'SP')
'n_mom',39,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'ic',@initial,... % initial conditions (defined below)
'ax',[-1 1 -1 1 -1 1]*0.6,... % coordinates of computational domain
'n',[1 1 1]*50,... % numbers of grid cells in each coordinate direction
'bc',[0 0 0],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,0.5,21),... % output times
'output',@output... % problem-specific output routine (defined below)
);
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
par = starmap_init(prob);  % Configure data structures for starmap solver
starmap_solver(par)                                         % Run solver.
%========================================================================
% Problem Specific Functions
%========================================================================
function f = initial(x,y,z,k)
% Initial conditions (for (k-1)-st moment).
t = 3.2e-4;              % Pseudo-time for smoothing initial Dirac delta.
f = 1/(4*pi*t)*exp(-(x.^2+y.^2+z.^2)/(4*t))*(k==1);

function f = sigma_a(x,y,z)
% Absorption coefficient.
f = 0;

function f = sigma_s0(x,y,z)
% Isotropic scattering coefficient.
f = 1;

function output(par,x,y,z,U,step)
% Plotting routine, showing the 3D results and a 1D cross-section.
t = par.t_plot(step);                                   % Time of output.
phi_max = max(max(max(U)));                         % Axis/color scaling.
xslice = 0; yslice = 0; zslice = 0;                    % Slice locations.
clf, subplot(1,3,1:2) 
slice(x,y,z,U,xslice,yslice,zslice)         % Volume slice approximation.
xlabel('x'), ylabel('y'), zlabel('z')
axis vis3d equal, caxis([-1 1]*phi_max)
title(sprintf('%s with %s%d at t = %0.2f',par.name,par.closure,par.n_mom,t))
colormap jet(255); colorbar
subplot(1,3,3)
r = linspace(0,max(x),100); 
% 1D cross sections
phix = interp3(x,y,z,U,r,r*0,r*0);                         % x>0,y=0,z=0.
phixy = interp3(x,y,z,U,r,r,r*0);                          % x>0,x=y,z=0.
phixyz = interp3(x,y,z,U,r,r,r);                           % x>0,x=y=z>0.
% Plot 1D cross sections.
plot(r,phix,'r-','DisplayName','x>0,y=0,z=0'), hold on 
plot(r,phixy,'b-','DisplayName','x>0,x=y,z=0'), hold on 
plot(r,phixyz,'k-','DisplayName','x>0,x=y=z'), hold off
axis([0 max(x) [-1 1]*phi_max])
legend('Location','SouthEast')
title('Cuts along different rays')
drawnow
