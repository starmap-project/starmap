function starmap_create_lung_cttest
%STARMAP_CREATE_LUNG_CTTEST
%   Creates an example test : starmap_2d_ex_lung_cttest_auto.m
%   which is an implementation of the Beam in 2D Patient CT
%   test problem from the paper 
%   "KiT-RT: An Extendable Framework for Radiative Transfer and Therapy."
%   https://doi.org/10.1145/3630001
% 
%   Version 2.01 (med)
%   Copyright (c) 09/27/2024 Benjamin Seibold, Martin Frank, and
%                            Rujeko Chinomona
%   http://www.math.temple.edu/~seibold
%   https://www.scc.kit.edu/personen/martin.frank.php
%   https://rujekoc.github.io/
%
%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0,v2.01),
%                 Pia Stammer (v2.01)
%
%   StaRMAP project website:
%   https://github.com/starmap-project

%   For license, see files LICENSE.txt or starmap_solver.m, as published on
%   https://github.com/starmap-project/starmap

%========================================================================
% Problem Parameters
%========================================================================
prob = struct(...
'name','Lung CT-Scan test',... % Name of example
'n_mom',13,... % Order of moment approximation
'ax',[0 1 0 1 0 1]*6,... % Coordinates of computational domain
'n',[40 40 1],... % Numbers of grid cells in each coordinate direction
'nangle',[30 60],...  % Num. grid points in polar angle, azimuthal angle 
'bc',[1 1 1],... % Type of boundary cond. (0 = periodic, 1 = extrapolation)
'E_plot',linspace(20,0,5)... % Output energies.
);

%========================================================================
% Initial Beam Parameters 
%========================================================================
mu_xy = [2.5,5.8]; sigma_xy = [0.1,0.1]; % spatial Gaussian parameters
mu_Omega2 = pi/2; sigma_Omega2 = 0.1; % angular Gaussian parameters
% Spatial Gaussian
space_beam = @(x,y) (1/(4*pi*sigma_xy(1)*sigma_xy(2)))*...
    exp(-((mu_xy(1) - x).^2)/(2*sigma_xy(1))).* ...
    exp(-((mu_xy(2) - y).^2)/(2*sigma_xy(2)));   
% Angular Gaussian
psi_angular = @(w,mu,sigma) 1/(sqrt(2*pi)*sigma)*exp((-(mu-w).^2)/(2*sigma^2));

%========================================================================
% Material parameters.   
%========================================================================
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
% Compute Initial Condition (Angular part) 
%========================================================================
n_mom = prob.n_mom;
n_sys = (n_mom+1)^2;
M = transformation_matrix(n_mom);  % Complex to real spherical harmonics

% Angular grids
ntheta = prob.nangle(1);     % number of grid points in polar angle 
nphi = prob.nangle(2);       % number of grid points in azimuthal angle 
theta = linspace(0, pi, ntheta);  % points in polar angle
phi = linspace(0, 2*pi, nphi);    % points in azimuthal angle 
[Theta_grid, Phi_grid] = meshgrid(theta, phi); % 2D angular grids
Mu_grid = cos(Theta_grid); % cosine of polar angle 
Omega2 = Theta_grid(:);
% Integration parameters 
dtheta = theta(2) - theta(1);
dphi = phi(2) - phi(1);
Weights = sin(Theta_grid(:))*dtheta*dphi;
% Angular Gaussian
G_angular = psi_angular(Omega2, mu_Omega2, sigma_Omega2).*Weights;
% Compute spherical harmonics
Y_lm_matrix = zeros(length(Theta_grid(:)),n_sys);
index = 1; 
for l = 0:n_mom
    for m = -l:l
        Y_lm = sph_cc(Mu_grid(:),Phi_grid(:),l,m); % returns complex conjugates
        Y_lm_matrix(:,index) = Y_lm;
        index = index + 1;
    end
end
u = Y_lm_matrix' * G_angular;       
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
function_name = 'starmap_2d_ex_lung_cttest_auto';
fprintf('Writing example file %s.',[function_name,'.m']);
fid = fopen([function_name,'.m'],'w');
fprintf(fid,'function %s\n',function_name);
fprintf(fid,'%s%s\n','%',upper(function_name));
fprintf(fid,'%s\n','%   Example case for STARMAP_SOLVER, a second order staggered');
fprintf(fid,'%s\n','%   grid finite difference solver for linear hyperbolic moment');
fprintf(fid,'%s\n','%   approximations to radiative transfer in 3D geometry.');
fprintf(fid,'%s\n','%   Implementation of the Beam in 2D Patient CT test problem from the paper:');
fprintf(fid,'%s\n','%   KiT-RT: An Extendable Framework for Radiative Transfer and Therapy.');
fprintf(fid,'%s\n','%   https://doi.org/10.1145/3630001');
fprintf(fid,'%s\n','%');
fprintf(fid,'%s%s%s\n','%   Created by the file ',mfilename,'.m');
fprintf(fid,'%s\n','%');
fprintf(fid,'%s\n','%   Version 2.01 (med)');
fprintf(fid,'%s\n','%   Copyright (c) 09/28/2024 Benjamin Seibold, Martin Frank, and');
fprintf(fid,'%s\n','%                            Rujeko Chinomona');
fprintf(fid,'%s\n','%   http://www.math.temple.edu/~seibold');
fprintf(fid,'%s\n','%   https://www.scc.kit.edu/personen/martin.frank.php');
fprintf(fid,'%s\n','%   https://rujekoc.github.io/');
fprintf(fid,'%s\n','%   ');
fprintf(fid,'%s\n','%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0,v2.01).');
fprintf(fid,'%s\n','%                     Pia Stammer (v2.01).');
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
fprintf(fid,'%s\n','''name'',''Lung CT-Scan test'',... % name of example');
fprintf(fid,'%s\n','''image_name'',''Lung_square.png'',... % name of image file');
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
fprintf(fid,'%s\n','par.tfinal = par.t_plot(end-1);');
fprintf(fid,'%s\n','par.sigma_s0 = @(x,y,z,t)par.sigma_s0(x,y,z,Time2Energy(t,E_CutOff));');
fprintf(fid,'%s\n','par.sigma_sm = @(x,y,z,m,t)par.sigma_sm(x,y,z,m,Time2Energy(t,E_CutOff));');
fprintf(fid,'%s\n','par.int_weight = @(m,t)StoppingPower(Time2Energy(t,E_CutOff)).*(m==1);');
fprintf(fid,'%s\n','par.image_matrix = processImage(par.image_name);');
fprintf(fid,'%s\n','figure');
fprintf(fid,'%s\n','solution = starmap_solver(par);');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','% Compute and plot depth dose'); 
fprintf(fid,'%s\n','x = solution(1).x; y = solution(1).y; z = solution(1).z;');
fprintf(fid,'%s\n','[X,Y,Z] = ndgrid(x,y,z);');
fprintf(fid,'%s\n','Rho = density(X,Y,Z,par);');
fprintf(fid,'%s\n','Dose = solution(1).Int./Rho;');
fprintf(fid,'%s\n',['figure, plot(x,Dose(:,',num2str(ceil(prob.n(2)/2)),',',num2str(ceil(prob.n(3)/2)),'))']);
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
fprintf(fid,'%s\n','[nx,ny] = size(par.image_matrix);');
fprintf(fid,'%s\n','xd = linspace(par.ax(1),par.ax(2),nx);');
fprintf(fid,'%s\n','yd = linspace(par.ax(3),par.ax(4),ny);');
fprintf(fid,'%s\n','[XD,YD] = meshgrid(xd,yd);');
fprintf(fid,'%s\n','f = interp2(XD,YD,par.image_matrix,x,y);');
fprintf(fid,'%s\n','% Scale densities to be between 0 (air) and 1.85 (bone)');
fprintf(fid,'%s\n','f=min(max(f * 1.85, 0.05), 1.85);');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = initial(x,y,z,k)');
fprintf(fid,'%s\n','% Initial conditions (for (k-1)-st moment).');
fprintf(fid,'%s\n',['sigma_xy = [',num2str(sigma_xy,'%g '),'];']');
fprintf(fid,'%s\n',['mu_xy = [',num2str(mu_xy,'%g '),'];']');
fprintf(fid,'%s\n',['f = feval(',func2str(space_beam),',x,y);']);
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
fprintf(fid,'%s\n',['plot(x,U(:,',num2str(ceil(prob.n(2)/2)),',',num2str(ceil(prob.n(3)/2)),')), hold on']);
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
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Process Image.');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','function I = processImage(image_name)');
fprintf(fid,'%s\n','% Read the image and convert to grayscale');
fprintf(fid,'%s\n','img = imread(image_name);');
fprintf(fid,'%s\n','if size(img, 3) == 3');
fprintf(fid,'%s\n','img = rgb2gray(img); % Convert RGB to grayscale');
fprintf(fid,'%s\n','end');
fprintf(fid,'%s\n','% Convert to double and scale values between 0 and 1');
fprintf(fid,'%s\n','I = double(img) / 255;');
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
function M = transformation_matrix(n_mom)
% Assemble transformation matrix
n_sys = (n_mom+1)^2;
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
