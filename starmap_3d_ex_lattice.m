function starmap_3d_ex_lattice
%STARMAP_3D_EX_LATTICE
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: 3D Lattice/Checkerboard geometry
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
close all
prob = struct(...
'name','Lattice Test 3D',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',7,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'source',@source,... % initial conditions (defined below)
'ax',[0 7 0 7 0 7],... % coordinates of computational domain
'n',[1 1 1]*40,... % numbers of grid cells in each coordinate direction
'bc',[1 1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,2.8,51),... % output times
'output',@output ... % problem-specific output routine (defined below)
);

%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
par = starmap_init(prob); % Configure data structures for starmap solver.
starmap_solver(par)                                         % Run solver.
%========================================================================
% Problem Specific Functions
%========================================================================
function f = source(x,y,z)
% Initial conditions (for (k-1)-st moment).
f = (3<x&x<4).*(3<y&y<4).*(3<z&z<4);

function f = sigma_a(x,y,z)
% Absorption coefficient.
cx = ceil(x); cy = ceil(y); cz = ceil(z);
g = (ceil((cx+cy)/2)*2==(cx+cy)).*(1<cx&cx<7&1<cy&cy-2*abs(cx-4)<4).*...
    (ceil((cx+cz)/2)*2==(cx+cz)).*(1<cx&cx<7&1<cz&cz<7);
f = (1-g)*0+g*10;

function f = sigma_s0(x,y,z)
% Isotropic scattering coefficient.
cx = ceil(x); cy = ceil(y); cz = ceil(z);
g = (ceil((cx+cy)/2)*2==(cx+cy)).*(1<cx&cx<7&1<cy&cy-2*abs(cx-4)<4).*...
    (ceil((cx+cz)/2)*2==(cx+cz)).*(1<cx&cx<7&1<cz&cz<7);
f =  (1-g)*1+g*0;

function output(par,x,y,z,U,step)
clf
subplot(1,6,1:2)
isovalue = 2e-3*max(max(max(U)));                 % Value for isosurface.
isosurface(x,y,z,U,isovalue)
grid on
axis equal, axis(par.ax)
view(3)
camlight
% light('Position',[3.5 3.5 -7],'Style','infinite')    % Add light sources.
% light('Position',[-7 3.5 3.5],'Style','infinite')
xlabel('x'), ylabel('y'), zlabel('z')
title(sprintf('Isosurface at %0.2d',isovalue))

subplot(1,6,3:5)
cax = [-7 0];                    % Colormap range used for log10 scaling.
vneg = cax(1)-diff(cax)/254;        % Value assigned where U is negative.
V = log10(max(U,1e-20));% Cap U s.t. U>0 and use logarithmic color scale.
Vcm = max(V,cax(1));                      % Cap colormap plot from below.
Vcm(U<0) = vneg;              % Assign special value where U is negative.
xslice = 3.5; yslice = 3.5; zslice = 3.5;              % Slice locations.
slice(x,y,z,Vcm,xslice,yslice,zslice)       % Volume slice approximation.
cm = jet(256); cm(1,:) = [1 1 1]*.5; % Change lowest color entry to gray.
axis vis3d equal
colormap(cm), colorbar('ylim',cax)
caxis([vneg cax(2)])
xlabel('x'), ylabel('y'), zlabel('z')
title(sprintf('Volume slice planes at \n x=3.5, y=3.5, z=3.5'))

subplot(1,6,6)
plot(x,interp3(x,y,z,V,x,y*0+3.5,z*0+3.5))
axis([par.ax(1:2) cax+[-1 1]*.1])
title('Cut at y=z=3.5'), xlabel('x')
sgtitle(sprintf('%s with %s%d at t = %0.2f',par.name,par.closure,...
    par.n_mom,par.t_plot(step)))
drawnow


