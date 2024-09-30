function starmap_create_lung_cttest
%STARMAP_CREATE_LUNG_CTTEST
%   Creates an example test : starmap_2d_ex_lung_cttest_auto.m
%   which is an implementation of the Beam in 2D Patient CT
%   test problem from the paper 
%   "KiT-RT: An Extendable Framework for Radiative Transfer and Therapy."
%   https://doi.org/10.1145/3630001
% 
%   Version 1.0-med
%   Copyright (c) 09/29/2024 Benjamin Seibold, Martin Frank, and
%                            Rujeko Chinomona
%   http://www.math.temple.edu/~seibold
%   https://www.scc.kit.edu/personen/martin.frank.php
%   https://rujekoc.github.io/
%   
%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0,v1.0-med).
%                     Pia Stammer (v1.0-med).
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
xslice = 2.5; yslice = 5;   % cross-section points
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
fprintf(fid,'%s\n','%   Version 1.0-med');
fprintf(fid,'%s\n','%   Copyright (c) 09/29/2024 Benjamin Seibold, Martin Frank, and');
fprintf(fid,'%s\n','%                            Rujeko Chinomona');
fprintf(fid,'%s\n','%   http://www.math.temple.edu/~seibold');
fprintf(fid,'%s\n','%   https://www.scc.kit.edu/personen/martin.frank.php');
fprintf(fid,'%s\n','%   https://rujekoc.github.io/');
fprintf(fid,'%s\n','%   ');
fprintf(fid,'%s\n','%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0,v1.0-med).');
fprintf(fid,'%s\n','%                     Pia Stammer (v1.0-med).');
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
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Dose computation and plotting');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Compute and plot depth dose'); 
fprintf(fid,'%s\n','x = solution(1).x; y = solution(1).y; z = solution(1).z;');
fprintf(fid,'%s\n','Dose = solution(1).Int./solution(1).Rho{1,1,1};');
fprintf(fid,'%s\n','Dose = Dose/max(Dose(:));  % normalized dose');
fprintf(fid,'%s\n',['xslice = ', num2str(xslice), '; yslice = ',num2str(yslice),';']);
fprintf(fid,'%s\n','% Store data to regenerate plots later');
fprintf(fid,'%s\n','ax = par.ax; save("DoseMatrix.mat","Dose","x","y","ax");');
fprintf(fid,'%s\n','ctImage = imread("Lung_square.png");  % Load ct scan');
fprintf(fid,'%s\n','if size(ctImage, 3) == 3,ctImage = rgb2gray(ctImage);end % Convert to gray');
fprintf(fid,'%s\n','[m,n] = size(ctImage); % size of image ');
fprintf(fid,'%s\n','xq = linspace(ax(1),ax(2),m); yq = linspace(ax(3),ax(4),n);');
fprintf(fid,'%s\n','[X,Y] = meshgrid(x,y);[XQ,YQ] = meshgrid(xq,yq);');
fprintf(fid,'%s\n','doseMatrix = interp2(X,Y,Dose,XQ,YQ);');
fprintf(fid,'%s\n','figure;');
fprintf(fid,'%s\n','ax1 = axes; imagesc(xq, yq, ctImage, ''Parent'', ax1); axis equal;');
fprintf(fid,'%s\n','set(gca, ''YTickLabel'', flip(get(gca, ''YTick'')));');
fprintf(fid,'%s\n','colormap(ax1, gray);  % Grayscale colormap for the image');
fprintf(fid,'%s\n','xlabel("x"); ylabel("y");');
fprintf(fid,'%s\n','title("Lung CT Scan with Dose Penetration");');
fprintf(fid,'%s\n','hold on;');
fprintf(fid,'%s\n','ax2 = axes;  % Create another set of axes for the contour plot');
fprintf(fid,'%s\n','hold(ax2, ''on'');');
fprintf(fid,'%s\n','contour(ax2, xq, yq, doseMatrix'', linspace(0, max(doseMatrix(:)),20), ''LineWidth'', 1.5);');
fprintf(fid,'%s\n','colormap(ax2, jet);  % Colored colormap for the contour plot');
fprintf(fid,'%s\n','colorbar(ax2);  % Add colorbar for the contour plot');
fprintf(fid,'%s\n','% Link the axes and set transparency for the second axes');
fprintf(fid,'%s\n','linkaxes([ax1, ax2]);');
fprintf(fid,'%s\n','ax2.Color = ''none'';  % Make ax2 transparent so ax1 can be seen');
fprintf(fid,'%s\n','ax2.XColor = ''none'';  % Hide x-axis for contour axes');
fprintf(fid,'%s\n','ax2.YColor = ''none'';  % Hide y-axis for contour axes');
fprintf(fid,'%s\n','ax2.Position = ax1.Position;  % Align the two axes');
fprintf(fid,'%s\n','% Add vertical and horizontal lines');
fprintf(fid,'%s\n','line([xslice, xslice], ylim, ''Color'', ''r'', ''LineWidth'', 2, ''LineStyle'', ''--'', ''Parent'', ax1);  % Vertical line');
fprintf(fid,'%s\n','line(xlim, [ax(4)-yslice, ax(4)-yslice], ''Color'', ''b'', ''LineWidth'', 2, ''LineStyle'', ''--'', ''Parent'', ax1);     % Horizontal line');
fprintf(fid,'%s\n','% Add labels for the lines');
fprintf(fid,'%s\n','text(xslice, ax(4), [''x = '', num2str(xslice), ''cm''], ''VerticalAlignment'', ''bottom'', ''HorizontalAlignment'', ''right'', ''Color'', ''r'', ''FontSize'', 12, ''Parent'', ax1);');
fprintf(fid,'%s\n','text(ax(2), ax(4)-yslice, [''y = '', num2str(yslice), ''cm''], ''VerticalAlignment'', ''bottom'', ''HorizontalAlignment'', ''left'', ''Color'', ''b'', ''FontSize'', 12, ''Parent'', ax1);');
fprintf(fid,'%s\n','hold off;');
fprintf(fid,'%s\n','% Dose along different slices');
fprintf(fid,'%s\n','dose_yslice = interp2(X,Y,Dose'',x,yslice);');
fprintf(fid,'%s\n','dose_xslice = interp2(X,Y,Dose'',xslice,y);');
fprintf(fid,'%s\n','figure;');
fprintf(fid,'%s\n','subplot(1,2,1), plot(y,dose_xslice,''LineWidth'',1.5)');
fprintf(fid,'%s\n','xlabel(''y [cm]''), ylabel(''normalized dose'')');
fprintf(fid,'%s\n','title([''Slice at x = '',num2str(xslice), '' cm''])');
fprintf(fid,'%s\n','subplot(1,2,2), plot(x,dose_yslice,''LineWidth'',1.5)');
fprintf(fid,'%s\n','xlabel(''x [cm]''), ylabel(''normalized dose'')');
fprintf(fid,'%s\n','title([''Slice at y = '',num2str(yslice), '' cm''])');
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
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Output function.');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','function output(par,x,y,z,U,step)');
fprintf(fid,'%s\n','% Output function showing the progress.');
fprintf(fid,'%s\n','E = par.E_plot(step);');
fprintf(fid,'%s\n','fprintf(''Energy:%12.2fMeV\n'',E)');
fprintf(fid,'%s\n','plot(x,U(:,ceil(par.n(2)/2),ceil(par.n(3)/2))'',''DisplayName'',[''E = '',num2str(E),''MeV''], ''LineWidth'',1.5)');
fprintf(fid,'%s\n','hold on; legend');
fprintf(fid,'%s\n','title([par.name,'': E = '',num2str(E),''MeV'']);');
fprintf(fid,'%s\n','xlabel(''x [cm]''), ylabel(''Radiative Intensity'')');
fprintf(fid,'%s\n','drawnow');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Energy transformation.');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','function E = Time2Energy(t,E_CutOff)');
fprintf(fid,'%s\n','% Transformation: Time to energy.');
fprintf(fid,'%s\n','E = max(0,energyTransform(energyTransform(E_CutOff,0)-t'',1))'';');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function t = Energy2Time(E,E_CutOff)');
fprintf(fid,'%s\n','% Transformation: Energy to time.');
fprintf(fid,'%s\n','t = max(0,energyTransform(E_CutOff-E'',0))'';');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function TE = energyTransform(E,inv)');
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
