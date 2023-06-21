function starmap_create_water_phantom
%STARMAP_CREATE_WATER_PHANTOM
%   Creates the file STARMAP_EX_WATER_PHANTOM.M, which is an
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
%   Version 2.01
%   Copyright (c) 06/28/2022 Benjamin Seibold, Martin Frank, and
%                            Rujeko Chinomona
%   http://www.math.temple.edu/~seibold
%   https://www.scc.kit.edu/personen/martin.frank.php
%   https://rujekoc.github.io/
%
%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0,v2.01).
%
%   StaRMAP project website:
%   https://github.com/starmap-project

%   For license, see files LICENSE.txt or starmap_solver.m, as published on
%   https://github.com/starmap-project/starmap

%========================================================================
% Problem Parameters
%========================================================================
prob = struct(...
'name','Phantom Test',... % Name of example
'n_mom',13,... % Order of moment approximation
'ax',[0 1 0 1 0 1]*10,... % Coordinates of computational domain
'n',[200 1 1],... % Numbers of grid cells in each coordinate direction
'bc',[1 1 1],... % Type of boundary cond. (0 = periodic, 1 = extrapolation)
'E_plot',linspace(5,0,6)... % Output energies.
);

% The beam is a Dirac in angle, times a Gaussian in space.
phi_beam = 0;                               % Angle of beam w.r.t. x-axis.
mu_beam = 0;                           % Cosine of the angle w.r.t z-axis.
pos_beam = [2,5,5];
space_beam = @(x,y,z)normpdf(x,pos_beam(1),.1).*... % Spatial distribution
    ((pos_beam(2)-5)<y&y<(pos_beam(2)+5)).*...              % of beam.
    ((pos_beam(3)-5)<z&z<(pos_beam(3)+5));          

% Material parameters.                        
E_tab = [5e-05;6e-05;7e-05;8e-05; % Energy range: 5e-5 MeV ... 1000 MeV 
    9e-05;0.0001;0.000125;0.00015;0.000175;0.0002;0.00025;0.0003;0.00035;
    0.0004;0.00045;0.0005;0.0006;0.0007;0.0008;0.0009;0.001;0.00125;
    0.0015;0.00175;0.002;0.0025;0.003;0.0035;0.004;0.00450;0.00500;
    0.00600;0.00700;0.00800;0.00900;0.0100;0.0125;0.0150;0.0175;0.0200;
    0.0250;0.0300;0.0350;0.0400;0.0450;0.0500;0.0600;0.0700;0.0800;0.0900;
    0.100;0.125;0.150;0.175;0.200;0.250;0.300;0.350;0.400;0.450;0.500;
    0.600;0.700;0.800;0.900;1;1.25;1.50;1.75;2;2.50;3;3.50;4;4.50;5;6;7;8;
    9;10;12.5;15;17.5;20;25;30;35;40;45;50;60;70;80;90;100;125;150;175;
    200;250;300;350;400;450;500;600;700;800;900;1000];
StoppingPower = @(E) StoppingPowerElectrons(E,278,'tot'); % Stopping Power.
TransportCoef = @(E,m) TransportCoefElectronsP39(E,'tot',m); % n-th transport coefficient.

%========================================================================
% Compute Initial Condition
%========================================================================
n_mom = prob.n_mom;
n_sys = (n_mom+1)^2;

% Assemble transformation matrix
M = sparse(zeros(n_sys)); s = size(M);
for m = 2:n_mom+1
    i = 1:m-1; r = (m-1)^2+2*i;
    M(sub2ind(s,r-1,(m-1)^2+i)) = 1;
    M(sub2ind(s,r,(m-1)^2+i)) = -1i;
    M(sub2ind(s,r-1,m^2+1-i)) = (-1).^(m+i);
    M(sub2ind(s,r,m^2+1-i)) = (-1).^(m+i)*1i;
end
M = M/sqrt(2);
m = 1:1:n_mom+1;
M(sub2ind(s,m.^2,(m-1).^2+m)) = 1;

% Assemble untransformed moment vector
u = zeros(n_sys,1);
for l = 0:n_mom
    u(l^2+1:(l+1)^2) = sph_cc(mu_beam,phi_beam+pi,l,-l:1:l)';
end

% Transform to StarMAP variables
StarMAPmoments = M*u;

%========================================================================
% Compute Energy transformation.
%========================================================================
S_tab = StoppingPower(E_tab);
E_trans = zeros(1,size(E_tab,1)); % Transformed energy.
for i = 2:size(E_tab,1)
    E_trans(i) = E_trans(i-1)+(E_tab(i)-E_tab(i-1))/2*(1/S_tab(i)+1/S_tab(i-1));
end

%========================================================================
% Write Example File
%========================================================================
function_name = 'starmap_ex_water_phantom_auto';
fprintf('Writing example file %s.',[function_name,'.m']);
fid = fopen([function_name,'.m'],'w');
fprintf(fid,'function %s\n',function_name);
fprintf(fid,'%s%s\n','%',upper(function_name));
fprintf(fid,'%s\n','%   Example case for STARMAP_SOLVER, a second order staggered');
fprintf(fid,'%s\n','%   grid finite difference solver for linear hyperbolic moment');
fprintf(fid,'%s\n','%   approximations to radiative transfer in 3D geometry.');
fprintf(fid,'%s\n','%');
fprintf(fid,'%s%s%s\n','%   Created by the file ',mfilename,'.m');
fprintf(fid,'%s\n','%');
fprintf(fid,'%s\n','%   Version 2.01');
fprintf(fid,'%s\n','%   Copyright (c) 06/28/2022 Benjamin Seibold, Martin Frank, and');
fprintf(fid,'%s\n','%                            Rujeko Chinomona');
fprintf(fid,'%s\n','%   http://www.math.temple.edu/~seibold');
fprintf(fid,'%s\n','%   https://www.scc.kit.edu/personen/martin.frank.php');
fprintf(fid,'%s\n','%   https://rujekoc.github.io/');
fprintf(fid,'%s\n','%   ');
fprintf(fid,'%s\n','%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.01).');
fprintf(fid,'%s\n','%   ');
fprintf(fid,'%s\n','%   StaRMAP project website:');
fprintf(fid,'%s\n','%   https://github.com/starmap-project');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%   For license, see files LICENSE.txt or starmap_solver.m, as published on');
fprintf(fid,'%s\n','%   https://github.com/starmap-project/starmap');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Problem Parameters');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','prob = struct(...');
fprintf(fid,'%s\n','''name'',''Phantom Test'',... % name of example');
fprintf(fid,'%s\n','''closure'',''P'',... % type of closure (can be ''P'' or ''SP'')');
fprintf(fid,'%s\n',['''n_mom'',',num2str(prob.n_mom),',... % order of moment approximation']);
fprintf(fid,'%s\n','''sigma_s0'',@sigma_s0,... % isotropic scattering coefficient (def. below)');
fprintf(fid,'%s\n','''sigma_sm'',@sigma_sm,... % aniso. scattering coefficient (defined below)');
fprintf(fid,'%s\n','''ic'',@initial,... % initial condition');
fprintf(fid,'%s\n',['''ax'',[',num2str(prob.ax,'%g '),'],... % coordinates of computational domain']);
fprintf(fid,'%s\n',['''n'',[',num2str(prob.n,'%g '),'],... % numbers of grid cells in each coordinate direction']);
fprintf(fid,'%s\n',['''bc'',[',num2str(prob.bc,'%g '), '],... % type of boundary cond. (0 = periodic, 1 = extrapolation)']);
fprintf(fid,'%s\n',['''E_plot'',[',num2str(prob.E_plot,'%g '),'],... % output times']);
fprintf(fid,'%s\n','''output'',@output,... % output routine (defined below)');
fprintf(fid,'%s\n','''density'',@density...');
fprintf(fid,'%s\n',');');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Moment System Setup and Solver Execution');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','par = starmap_init(prob);     % Configure data structures for starmap solver');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','% Modify functions and run solver.');
fprintf(fid,'%s\n','E_CutOff = max(par.E_plot);');
fprintf(fid,'%s\n','par.t_plot = Energy2Time(par.E_plot,E_CutOff);');
fprintf(fid,'%s\n','par.sigma_s0 = @(x,y,z,t)par.sigma_s0(x,y,z,Time2Energy(t,E_CutOff));');
fprintf(fid,'%s\n','par.sigma_sm = @(x,y,z,m,t)par.sigma_sm(x,y,z,m,Time2Energy(t,E_CutOff));');
fprintf(fid,'%s\n','par.int_weight = @(m,t)StoppingPower(Time2Energy(t,E_CutOff)).*(m==1);');
fprintf(fid,'%s\n','figure');
fprintf(fid,'%s\n','solution = starmap_solver(par);');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','% Compute and plot depth dose'); 
fprintf(fid,'%s\n','x = solution(1).x; y = solution(1).y; z = solution(1).z;');
fprintf(fid,'%s\n','[X,Y,Z] = ndgrid(x,y,z);');
fprintf(fid,'%s\n','Rho = density(X,Y,Z,par);');
fprintf(fid,'%s\n','Dose = solution(1).Int./Rho;');
fprintf(fid,'%s\n',['figure, plot(x,Dose(:,',num2str(ceil(par.n(2)/2)),',',num2str(ceil(par.n(3)/2)),'))']);
fprintf(fid,'%s\n','xlabel(''x [cm]''), ylabel(''dose'')');
fprintf(fid,'%s\n','title([par.name,'': depth dose''])');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Problem Specific Functions');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','function f = sigma_s0(x,y,z,E)');
fprintf(fid,'%s\n','% Total scattering coefficient.');
fprintf(fid,'%s\n',['f = feval(',func2str(TransportCoef),',E,0);']);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = sigma_sm(x,y,z,m,E)');
fprintf(fid,'%s\n','% Moments of scattering kernel.');
fprintf(fid,'%s\n',['f = feval(',func2str(TransportCoef),',E,m);']);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = density(x,y,z,par)');
fprintf(fid,'%s\n','% Problem specific density function.');
fprintf(fid,'%s\n','f = 1;');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = initial(x,y,z,k)');
fprintf(fid,'%s\n','% Initial conditions (for (k-1)-st moment).');
fprintf(fid,'%s\n',['pos_beam = [',num2str(pos_beam,'%g '),'];']');
fprintf(fid,'%s\n',['f = feval(',func2str(space_beam),',x,y,z);']);
fprintf(fid,'%s\n','StarMAPmoments = ['); 
fprintf(fid,'%12.8f\n',StarMAPmoments);
fprintf(fid,'%s\n','];');
fprintf(fid,'%s\n','f = f*StarMAPmoments(k);');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = StoppingPower(E)');
fprintf(fid,'%s\n','% Stopping Power.');
fprintf(fid,'%s\n',['f = feval(',func2str(StoppingPower),',E);']);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function output(par,x,y,z,U,step)');
fprintf(fid,'%s\n','% Output function showing the progress.');
fprintf(fid,'%s\n','E = par.E_plot(step);');
fprintf(fid,'%s\n','fprintf(''Energy:%12.2fMeV\n'',E)');
fprintf(fid,'%s\n',['plot(x,U(:,',num2str(ceil(par.n(2)/2)),',',num2str(ceil(par.n(3)/2)),')), hold on']);
fprintf(fid,'%s\n','title([par.name,'': E = '',num2str(E),''MeV'']);');
fprintf(fid,'%s\n','xlabel(''x [cm]''), ylabel(''zeroth-moment'')');
fprintf(fid,'%s\n','drawnow');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Energy transformation.');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','function E = Time2Energy(t,E_CutOff)');
fprintf(fid,'%s\n','% Transformation: Time to energy.');
fprintf(fid,'%s\n','E = max(0,energyTansform(energyTansform(E_CutOff,0)-t'',1))'';');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function t = Energy2Time(E,E_CutOff)');
fprintf(fid,'%s\n','% Transformation: Energy to time.');
fprintf(fid,'%s\n','t = max(0,energyTansform(E_CutOff-E'',0))'';');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function TE = energyTansform(E,inv)');
fprintf(fid,'%s\n','% Transform the energy using linear interpolation.');
fprintf(fid,'%s\n','E_tab = [');
fprintf(fid,'%12.8f\n',E_tab);
fprintf(fid,'%s\n','];');
fprintf(fid,'%s\n','E_trans = [');
fprintf(fid,'%12.8f\n',E_trans);
fprintf(fid,'%s\n','];');
fprintf(fid,'%s\n','if inv==0');
fprintf(fid,'%s\n','    TE = interp1q(E_tab,E_trans,E);');
fprintf(fid,'%s\n','else');
fprintf(fid,'%s\n','    TE = interp1q(E_trans,E_tab,E);');
fprintf(fid,'%s\n','end');
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
