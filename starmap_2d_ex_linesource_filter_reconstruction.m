 function starmap_2d_ex_linesource_filter_reconstruction
%STARMAP_2D_EX_LINESOURCE_RECONSTRUCTION is a modification of
%STARMAP_2D_EX_LINESOURCE.
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: Line source benchmark test case, as defined in
%   [Ganapol, JCP 210 (2005) 386-399].
%   Test case is adapted to using a narrow Gaussian as initial
%   condition. It includes the angular reconstruction in a
%   self-chosen point. The positive and negative part of the
%   solution is shown to investigate the quality of the
%   solution in this point. Feel free to modify the code to
%   PN by removing the filtering.
%   This routine is specific to the PN equations.
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
'name','Line Source Test',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',13,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'ic',@initial,... % initial conditions (defined below)
'ax',[-1 1 -1 1]*0.6,... % coordinates of computational domain
'n',[150 150],... % numbers of grid cells in each coordinate direction
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,0.5,51),... % output times
'output',@output,... % problem-specific output routine (defined below)
'mom_output','all',... % moments plotted
'filterfunction',@exp_filter,... % filter (defined below)
'f_position','substep',...% position of filter ('substep', 'full', or '')
'point',[0.1 0.2],... % point for angular reconstruction
'n_recon',100 ... % resolution of the reconstruction
);
%========================================================================
% Initialize variables for angular reconstruction
%========================================================================
% Assemble transformation matrix
n_sys = (prob.n_mom+1)*(prob.n_mom+2)/2;
prob.S = sparse(zeros(n_sys)); s = size(prob.S);
for m = 2:prob.n_mom+1
    i = 1:floor(m/2); r = m*(m-1)/2+2*i;
    prob.S(sub2ind(s,r-1,m*(m-1)/2+i)) = 1;
    prob.S(sub2ind(s,r,m*(m-1)/2+i)) = -1i;
    prob.S(sub2ind(s,r-1,m*(m+1)/2+1-i)) = (-1)^(m+1);
    prob.S(sub2ind(s,r,m*(m+1)/2+1-i)) = (-1)^(m+1)*1i;
end
prob.S = prob.S/sqrt(2);
m = 1:2:prob.n_mom+1;
prob.S(sub2ind(s,m.*(m+1)/2,(m.^2+1)/2)) = 1;
% Create grid and compute spherical harmonics
moment_order = ceil(sqrt(2*(1:n_sys) + 1/4)-3/2);
moment_vec = 2*((1:n_sys)-1)-moment_order.*(moment_order+2);
delta = pi/prob.n_recon;
theta = 0:delta:pi; phi = 0:2*delta:2*pi;
prob.Ynm = cell(1,n_sys);
for j = 1:n_sys
    prob.Ynm{j} = sph(moment_order(j),moment_vec(j),theta,phi);
end
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
par = starmap_init(prob);  % Configure data structures for starmap solver
starmap_solver(par)                 % moment matrices and run solver.

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

function output(par,x,y,U,step)
% Plotting routine, showing the 2D results and the 3D angular
% reconstruction in the self-chosen point. The positive part of the
% solution is shown in red and the negative part in green (the black
% lines show the resolution of the reconstruction).
t = par.t_plot(step);                                   % Time of output.
n_sys = length(U);
U_point = zeros(n_sys,1);        % Compute solution in given space-point.
for j = 1:n_sys
    U_point(j) = interp2(y{j},x{j},U{j},par.point(2),par.point(1));
end
U_point = par.S\U_point;            % Transform to complex sph expansion.
U_sol = real(sumcell(par.Ynm, U_point));
clf, subplot(1,5,1:3), imagesc(x{1},y{1},U{1}'),             % 2D plot of
axis xy equal tight,                                     % approximation.
title(sprintf('%s with FP%d at t = %0.2f',par.name,par.n_mom,t))
colormap jet(255); colorbar
hold on, plot(par.point(1),par.point(2),'.r'), hold off
text(par.point(1),par.point(2),['\color{red} \leftarrow (',...
    num2str(par.point(1)),',',num2str(par.point(2)),')'])
subplot(1,5,4:5),  hold on           % 3D plot of angular reconstruction.
[Surf_pos,Surf_neg] = Surfaces(U_sol,par.n_recon);
set(Surf_pos,'FaceColor',[1 0 0]), set(Surf_neg,'FaceColor',[0 1 0])
title(sprintf('Reconstruction in (%0.1f,%0.1f)',...
    par.point(1),par.point(2)))
axis tight, view(2), hold off % Set the view to 2D to show the x-y-plane.
drawnow

function y = sph(l,m,theta,phi)
% Spherical harmonics.
z = legendre(l,-cos(theta));
if m>=0, y = sqrt((2*l+1)./(4*pi).*factorial(l-m)./factorial(l+m)).*...
        (-1)^m.*z(m+1,:)'*exp(1i*m*phi);
else, y=(-1)^m.*conj(sph(l,-m,theta,phi));
end

function S = sumcell(A,w)
% Add vector of cells A, weighted with vector w.
S = A{1}*w(1); for j = 2:length(w), S = S+A{j}*w(j); end

function [Surf_pos,Surf_neg] = Surfaces(U,n)
% Compute positive and negative part of U and generate surfaces.
U_pos = 0.5*(abs(U)+U);                % Split solution into positive and
U_neg = 0.5*(abs(U)-U);                                  % negative part.
delta = pi/n; theta = 0:delta:pi; phi = 0:2*delta:2*pi;    % Create grid.
[phi,theta] = meshgrid(phi,-theta);
r = U_pos.*sin(theta);                 % Change to cartesian coordinates.
x = r.*cos(phi); y = r.*sin(phi); z = U_pos.*cos(theta);
Surf_pos = surf(x,y,z);                               % Generate surface.
r = U_neg.*sin(theta);
x = r.*cos(phi); y = r.*sin(phi); z = U_neg.*cos(theta);
Surf_neg = surf(x,y,z);
