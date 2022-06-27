function starmap_2d_create_beam
%STARMAP_2D_CREATE_BEAM
%   Creates the file STARMAP_2D_EX_BEAM_AUTO.M, which is an
%   example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   An initial beam in angle and a spatial distribution can be
%   specified. The routine generates an initial condition and
%   writes default functions for all other problem parameters
%   that can later be changed in the generated example file.
%   This routine is specific to the PN equations.
%
%   Note: By default, the domain contains a void. In that
%   case, the moment approximation cannot be used to solve for
%   a steady state solution.
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
'name','Beam Test',... % name of example
'n_mom',9,... % order of moment approximation
'ax',[-1 1 -1 1]*0.6,... % coordinates of computational domain
'n',[1 1]*150,... % numbers of grid cells in each coordinate direction
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,0.6,61)... % output times
);

% The beam is a Dirac in angle, times a Gaussian in space.
phi_beam = pi/6;                           % Angle of beam w.r.t. x-axis.
space_beam = @(x,y) 1/(4*pi*3.2e-4)*...            % Spatial distribution
    exp(-(x.^2+y.^2)/(4*3.2e-4));                              % of beam.

%========================================================================
% Compute Initial Condition
%========================================================================
n_mom = prob.n_mom;
n_sys = (n_mom+1)*(n_mom+2)/2;
mu_beam = 0;

% Assemble transformation matrix
M = sparse(zeros(n_sys)); s = size(M);
for m = 2:n_mom+1
    i = 1:floor(m/2); r = m*(m-1)/2+2*i;
    M(sub2ind(s,r-1,m*(m-1)/2+i)) = 1;
    M(sub2ind(s,r,m*(m-1)/2+i)) = -1i;
    M(sub2ind(s,r-1,m*(m+1)/2+1-i)) = (-1)^(m+1);
    M(sub2ind(s,r,m*(m+1)/2+1-i)) = (-1)^(m+1)*1i;
end
M = M/sqrt(2);
m = 1:2:n_mom+1;
M(sub2ind(s,m.*(m+1)/2,(m.^2+1)/2)) = 1;

% Assemble untransformed moment vector
u = zeros(n_sys,1);
for l = 0:n_mom
    u(l*(l+1)/2+1:(l+1)*(l+2)/2) = sph_cc(mu_beam,phi_beam+pi,l,-l:2:l)';
end

% Transform to StarMAP variables
StarMAPmoments = M*u;

%========================================================================
% Write Example File
%========================================================================
function_name = 'starmap_2d_ex_beam_auto';
fprintf('Writing example file %s.',[function_name,'.m']);
fid = fopen([function_name,'.m'],'w');
fprintf(fid,'function %s\n',function_name);
fprintf(fid,'%s%s\n','%',upper(function_name));
fprintf(fid,'%s\n','%   Example case for STARMAP_SOLVER, a second order staggered');
fprintf(fid,'%s\n','%   grid finite difference solver for linear hyperbolic moment');
fprintf(fid,'%s\n','%   approximations to radiative transfer in 2D slab geometry.');
fprintf(fid,'%s\n','%');
fprintf(fid,'%s%s%s\n','%   Created by the file ',mfilename,'.m');
fprintf(fid,'%s\n','%');
fprintf(fid,'%s\n','%   Version 2.0');
fprintf(fid,'%s\n','%   Copyright (c) 06/28/2022 Benjamin Seibold, Martin Frank, and');
fprintf(fid,'%s\n','%                            Rujeko Chinomona');
fprintf(fid,'%s\n','%   http://www.math.temple.edu/~seibold');
fprintf(fid,'%s\n','%   https://www.scc.kit.edu/personen/martin.frank.php');
fprintf(fid,'%s\n','%   https://rujekoc.github.io/');
fprintf(fid,'%s\n','%   ');
fprintf(fid,'%s\n','%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0).');
fprintf(fid,'%s\n','%   ');
fprintf(fid,'%s\n','%   StaRMAP project website:');
fprintf(fid,'%s\n','%   https://github.com/starmap-project');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%   For license, see file starmap_solver.m, as published on');
fprintf(fid,'%s\n','%   https://github.com/starmap-project/starmap');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Problem Parameters');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','prob = struct(...');
fprintf(fid,'%s\n','''name'',''Beam Test'',... % name of example');
fprintf(fid,'%s\n','''closure'',''P'',... % type of closure (can be ''P'' or ''SP'')');
fprintf(fid,'%s\n',['''n_mom'',',num2str(prob.n_mom),',... % order of moment approximation']);
fprintf(fid,'%s\n','''sigma_a'',@sigma_a,... % absorption coefficient (defined below)');
fprintf(fid,'%s\n','''sigma_s0'',@sigma_s0,... % isotropic scattering coefficient (def. below)');
fprintf(fid,'%s\n','''sigma_sm'',@sigma_sm,... % aniso. scattering coefficient (defined below)');
fprintf(fid,'%s\n','''source'',@source,... % source term (defined below)');
fprintf(fid,'%s\n','''ic'',@initial,... % initial condition');
fprintf(fid,'%s\n',['''ax'',[',num2str(prob.ax,'%g '),'],... % coordinates of computational domain']);
fprintf(fid,'%s\n',['''n'',[',num2str(prob.n,'%g '),'],... % numbers of grid cells in each coordinate direction']);
fprintf(fid,'%s\n',['''bc'',[',num2str(prob.bc,'%g '), '],... % type of boundary cond. (0 = periodic, 1 = extrapolation)']);
fprintf(fid,'%s\n',['''t_plot'',[',num2str(prob.t_plot,'%g '),'],... % output times']);
fprintf(fid,'%s\n','''output'',@output... % output routine (defined below)');
fprintf(fid,'%s\n',');');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Moment System Setup and Solver Execution');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','par = starmap_init(prob);     % Configure data structures for starmap solver');
fprintf(fid,'%s\n','starmap_solver(par)                                          % run solver');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Problem Specific Functions');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','function f = sigma_a(x,y)');
fprintf(fid,'%s\n','% Absorption coefficient.');
fprintf(fid,'%s\n','f = 0;');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = sigma_s0(x,y)');
fprintf(fid,'%s\n','% Total scattering coefficient.');
fprintf(fid,'%s\n','f = (x>.3)*100;');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = sigma_sm(x,y,m)');
fprintf(fid,'%s\n','% Moments of scattering kernel.');
fprintf(fid,'%s\n','g = 0.85;');
fprintf(fid,'%s\n','f = (x>.3)*100*g^m;');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = initial(x,y,k)');
fprintf(fid,'%s\n','% Initial conditions (for (k-1)-st moment).');
fprintf(fid,'%s\n','f = 0;');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = source(x,y,t,k)');
fprintf(fid,'%s\n','% Source (for (k-1)-st moment).');
fprintf(fid,'%s\n',['f = feval(',func2str(space_beam),',x,y);']);
fprintf(fid,'%s\n','StarMAPmoments = [');
fprintf(fid,'%12.8f\n',StarMAPmoments);
fprintf(fid,'%s\n','];');
fprintf(fid,'%s\n','f = f*StarMAPmoments(k);');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function output(par,x,y,U,step)');
fprintf(fid,'%s\n','% Plotting routine.');
fprintf(fid,'%s\n','t = par.t_plot(step);');
fprintf(fid,'%s\n','clf, imagesc(x,y,U'') % 2D plot of approximation');
fprintf(fid,'%s\n','hold on, plot([1 1]*.3,y([1 end]),''k:''), hold off');
fprintf(fid,'%s\n','axis xy equal tight, title(sprintf(''%s with %s%d at t = %0.2f'',...');
fprintf(fid,'%s\n','    par.name,par.closure,par.n_mom,t))');
fprintf(fid,'%s\n','colormap jet(255); colorbar;');
fprintf(fid,'%s\n','caxis([-2 6]);');
fprintf(fid,'%s\n','drawnow');
fclose(fid);
fprintf(' Done.\n')

%========================================================================
% Functions
%========================================================================
function y = sph_cc(mu,phi,l,m)
% Complex conjugates of coefficients.
z = legendre(l,mu)'; ma = abs(m);
y = sqrt((2*l+1)/(4*pi).*factorial(l-ma)./factorial(l+ma)).*...
    (-1).^max(m,0).*exp(1i*m*phi).*z(ma+1);
