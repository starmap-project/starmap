function starmap_2d_ex_lung_cttest_auto
%STARMAP_2D_EX_LUNG_CTTEST_AUTO
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 3D geometry.
%   Implementation of the Beam in 2D Patient CT test problem from the paper:
%   KiT-RT: An Extendable Framework for Radiative Transfer and Therapy.
%   https://doi.org/10.1145/3630001
%
%   Created by the file starmap_create_lung_cttest.m
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
'name','Lung CT-Scan test',... % name of example
'image_name','Lung_square.png',... % name of image file
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',13,... % order of moment approximation
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'sigma_sm',@sigma_sm,... % aniso. scattering coefficient (defined below)
'ic',@initial,... % initial condition
'ax',[0 6 0 6 0 6],... % coordinates of computational domain
'n',[40 40  1],... % numbers of grid cells in each coordinate direction
'bc',[1 1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
'E_plot',[20 15 10  5  0],... % output times
'output',@output,... % output routine (defined below)
'density',@density...
);

%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
par = starmap_init(prob);     % Configure data structures for starmap solver

% Modify functions and run solver.
E_CutOff = max(par.E_plot);
par.t_plot = Energy2Time(par.E_plot,E_CutOff);
par.tfinal = par.t_plot(end-1);
par.sigma_s0 = @(x,y,z,t)par.sigma_s0(x,y,z,Time2Energy(t,E_CutOff));
par.sigma_sm = @(x,y,z,m,t)par.sigma_sm(x,y,z,m,Time2Energy(t,E_CutOff));
par.int_weight = @(m,t)StoppingPower(Time2Energy(t,E_CutOff)).*(m==1);
par.image_matrix = processImage(par.image_name);
figure
solution = starmap_solver(par);

%========================================================================
% Dose computation and plotting
%========================================================================
% Compute and plot depth dose
x = solution(1).x; y = solution(1).y; z = solution(1).z;
Dose = solution(1).Int./solution(1).Rho{1,1,1};
Dose = Dose/max(Dose(:));  % normalized dose
xslice = 2.5; yslice = 5;
% Store data to regenerate plots later
ax = par.ax; save("DoseMatrix.mat","Dose","x","y","ax");
ctImage = imread("Lung_square.png");  % Load ct scan
if size(ctImage, 3) == 3,ctImage = rgb2gray(ctImage);end % Convert to gray
[m,n] = size(ctImage); % size of image 
xq = linspace(ax(1),ax(2),m); yq = linspace(ax(3),ax(4),n);
[X,Y] = meshgrid(x,y);[XQ,YQ] = meshgrid(xq,yq);
doseMatrix = interp2(X,Y,Dose,XQ,YQ);
figure;
ax1 = axes; imagesc(xq, yq, ctImage, 'Parent', ax1); axis equal;
set(gca, 'YTickLabel', flip(get(gca, 'YTick')));
colormap(ax1, gray);  % Grayscale colormap for the image
xlabel("x"); ylabel("y");
title("Lung CT Scan with Dose Penetration");
hold on;
ax2 = axes;  % Create another set of axes for the contour plot
hold(ax2, 'on');
contour(ax2, xq, yq, doseMatrix', linspace(0, max(doseMatrix(:)),20), 'LineWidth', 1.5);
colormap(ax2, jet);  % Colored colormap for the contour plot
colorbar(ax2);  % Add colorbar for the contour plot
% Link the axes and set transparency for the second axes
linkaxes([ax1, ax2]);
ax2.Color = 'none';  % Make ax2 transparent so ax1 can be seen
ax2.XColor = 'none';  % Hide x-axis for contour axes
ax2.YColor = 'none';  % Hide y-axis for contour axes
ax2.Position = ax1.Position;  % Align the two axes
% Add vertical and horizontal lines
line([xslice, xslice], ylim, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--', 'Parent', ax1);  % Vertical line
line(xlim, [ax(4)-yslice, ax(4)-yslice], 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--', 'Parent', ax1);     % Horizontal line
% Add labels for the lines
text(xslice, ax(4), ['x = ', num2str(xslice), 'cm'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r', 'FontSize', 12, 'Parent', ax1);
text(ax(2), ax(4)-yslice, ['y = ', num2str(yslice), 'cm'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'b', 'FontSize', 12, 'Parent', ax1);
hold off;
% Dose along different slices
dose_yslice = interp2(X,Y,Dose',x,yslice);
dose_xslice = interp2(X,Y,Dose',xslice,y);
figure;
subplot(1,2,1), plot(y,dose_xslice,'LineWidth',1.5)
xlabel('y [cm]'), ylabel('normalized dose')
title(['Slice at x = ',num2str(xslice), ' cm'])
subplot(1,2,2), plot(x,dose_yslice,'LineWidth',1.5)
xlabel('x [cm]'), ylabel('normalized dose')
title(['Slice at y = ',num2str(yslice), ' cm'])
%========================================================================
% Problem Specific Functions
%========================================================================
function f = sigma_s0(x,y,z,E)
% Total scattering coefficient.
f = feval(@(E,m)TransportCoefElectronsP39(E,'tot',m),E,0);

function f = sigma_sm(x,y,z,m,E)
% Moments of scattering kernel.
f = feval(@(E,m)TransportCoefElectronsP39(E,'tot',m),E,m);

function f = density(x,y,z,par)
% Problem specific density function.
[nx,ny] = size(par.image_matrix);
xd = linspace(par.ax(1),par.ax(2),nx);
yd = linspace(par.ax(3),par.ax(4),ny);
[XD,YD] = meshgrid(xd,yd);
f = interp2(XD,YD,par.image_matrix,x,y);
% Scale densities to be between 0 (air) and 1.85 (bone)
f=min(max(f * 1.85, 0.05), 1.85);

function f = initial(x,y,z,k)
% Initial conditions (for (k-1)-st moment).
sigma_xy = [0.1 0.1];
mu_xy = [2.5 5.8];
f = feval(@(x,y)(1/(4*pi*sigma_xy(1)*sigma_xy(2)))*exp(-((mu_xy(1)-x).^2)/(2*sigma_xy(1))).*exp(-((mu_xy(2)-y).^2)/(2*sigma_xy(2))),x,y);
StarMAPmoments = [
  1.79350525
  0.05177404
 -0.00000000
  3.10644222
  0.01929504
 -0.00000000
  0.03859009
 -0.00000000
  4.01039966
  0.00416821
  0.00000000
  0.01020998
  0.00000000
  0.03228678
 -0.00000000
  4.74516887
  0.00063158
 -0.00000000
  0.00178637
  0.00000000
  0.00668400
 -0.00000000
  0.02835781
 -0.00000000
  5.38051575
  0.00007360
 -0.00000000
  0.00023275
 -0.00000000
  0.00098746
 -0.00000000
  0.00483752
  0.00000000
  0.02559777
 -0.00000000
  5.94838397
  0.00000696
  0.00000000
  0.00002412
 -0.00000000
  0.00011315
 -0.00000000
  0.00061977
  0.00000000
  0.00371863
 -0.00000000
  0.02351871
 -0.00000000
  6.46657514
  0.00000055
  0.00000000
  0.00000207
  0.00000000
  0.00001058
 -0.00000000
  0.00006348
 -0.00000000
  0.00042105
 -0.00000000
  0.00297729
  0.00000000
  0.02187852
 -0.00000000
  6.94621596
  0.00000004
 -0.00000000
  0.00000015
  0.00000000
  0.00000083
  0.00000000
  0.00000541
 -0.00000000
  0.00003901
 -0.00000000
  0.00030221
 -0.00000000
  0.00245514
  0.00000000
  0.02054114
 -0.00000000
  7.39481158
  0.00000000
  0.00000000
  0.00000001
 -0.00000000
  0.00000006
  0.00000000
  0.00000039
  0.00000000
  0.00000306
 -0.00000000
  0.00002558
 -0.00000000
  0.00022591
  0.00000000
  0.00207053
  0.00000000
  0.01942325
  0.00000000
  7.81770814
  0.00000000
 -0.00000000
  0.00000000
  0.00000000
  0.00000000
 -0.00000000
  0.00000003
  0.00000000
  0.00000021
  0.00000000
  0.00000186
 -0.00000000
  0.00001761
 -0.00000000
  0.00017428
  0.00000000
  0.00177733
  0.00000000
  0.01847054
 -0.00000000
  8.21887357
  0.00000000
 -0.00000000
  0.00000000
 -0.00000000
  0.00000000
  0.00000000
  0.00000000
 -0.00000000
  0.00000001
  0.00000000
  0.00000012
  0.00000000
  0.00000119
 -0.00000000
  0.00001259
 -0.00000000
  0.00013788
  0.00000000
  0.00154765
  0.00000000
  0.01764588
 -0.00000000
  8.60134901
  0.00000000
  0.00000000
  0.00000000
 -0.00000000
  0.00000000
 -0.00000000
  0.00000000
  0.00000000
  0.00000000
 -0.00000000
  0.00000001
  0.00000000
  0.00000007
  0.00000000
  0.00000080
 -0.00000000
  0.00000928
 -0.00000000
  0.00011134
  0.00000000
  0.00136368
 -0.00000000
  0.01692287
 -0.00000000
  8.96752625
  0.00000000
  0.00000000
  0.00000000
  0.00000000
  0.00000000
 -0.00000000
  0.00000000
 -0.00000000
  0.00000000
  0.00000000
  0.00000000
 -0.00000000
  0.00000000
  0.00000000
  0.00000004
  0.00000000
  0.00000055
 -0.00000000
  0.00000702
 -0.00000000
  0.00009148
 -0.00000000
  0.00121360
  0.00000000
  0.01628218
 -0.00000000
  9.31932665
];
f = f*StarMAPmoments(k);

function f = StoppingPower(E)
% Stopping Power.
f = feval(@(E)StoppingPowerElectrons(E,278,'tot'),E);

%========================================================================
% Output function.
%========================================================================
function output(par,x,y,z,U,step)
% Output function showing the progress.
E = par.E_plot(step);
fprintf('Energy:%12.2fMeV\n',E)
plot(x,U(:,ceil(par.n(2)/2),ceil(par.n(3)/2))','DisplayName',['E = ',num2str(E),'MeV'], 'LineWidth',1.5)
hold on; legend
title([par.name,': E = ',num2str(E),'MeV']);
xlabel('x [cm]'), ylabel('Radiative Intensity')
drawnow

%========================================================================
% Energy transformation.
%========================================================================
function E = Time2Energy(t,E_CutOff)
% Transformation: Time to energy.
E = max(0,energyTransform(energyTransform(E_CutOff,0)-t',1))';

function t = Energy2Time(E,E_CutOff)
% Transformation: Energy to time.
t = max(0,energyTransform(E_CutOff-E',0))';

function TE = energyTransform(E,inv)
% Transform the energy using linear interpolation.
E_tab = [
  0.00005000
  0.00006000
  0.00007000
  0.00008000
  0.00009000
  0.00010000
  0.00012500
  0.00015000
  0.00017500
  0.00020000
  0.00025000
  0.00030000
  0.00035000
  0.00040000
  0.00045000
  0.00050000
  0.00060000
  0.00070000
  0.00080000
  0.00090000
  0.00100000
  0.00125000
  0.00150000
  0.00175000
  0.00200000
  0.00250000
  0.00300000
  0.00350000
  0.00400000
  0.00450000
  0.00500000
  0.00600000
  0.00700000
  0.00800000
  0.00900000
  0.01000000
  0.01250000
  0.01500000
  0.01750000
  0.02000000
  0.02500000
  0.03000000
  0.03500000
  0.04000000
  0.04500000
  0.05000000
  0.06000000
  0.07000000
  0.08000000
  0.09000000
  0.10000000
  0.12500000
  0.15000000
  0.17500000
  0.20000000
  0.25000000
  0.30000000
  0.35000000
  0.40000000
  0.45000000
  0.50000000
  0.60000000
  0.70000000
  0.80000000
  0.90000000
  1.00000000
  1.25000000
  1.50000000
  1.75000000
  2.00000000
  2.50000000
  3.00000000
  3.50000000
  4.00000000
  4.50000000
  5.00000000
  6.00000000
  7.00000000
  8.00000000
  9.00000000
 10.00000000
 12.50000000
 15.00000000
 17.50000000
 20.00000000
 25.00000000
 30.00000000
 35.00000000
 40.00000000
 45.00000000
 50.00000000
 60.00000000
 70.00000000
 80.00000000
 90.00000000
100.00000000
125.00000000
150.00000000
175.00000000
200.00000000
250.00000000
300.00000000
350.00000000
400.00000000
450.00000000
500.00000000
600.00000000
700.00000000
800.00000000
900.00000000
1000.00000000
];
E_trans = [
  0.00000000
  0.00000003
  0.00000005
  0.00000008
  0.00000010
  0.00000012
  0.00000019
  0.00000025
  0.00000032
  0.00000040
  0.00000056
  0.00000074
  0.00000094
  0.00000115
  0.00000139
  0.00000164
  0.00000219
  0.00000281
  0.00000349
  0.00000423
  0.00000503
  0.00000727
  0.00000982
  0.00001268
  0.00001583
  0.00002302
  0.00003133
  0.00004073
  0.00005118
  0.00006265
  0.00007514
  0.00010303
  0.00013474
  0.00017015
  0.00020916
  0.00025168
  0.00037280
  0.00051443
  0.00067574
  0.00085601
  0.00127074
  0.00175467
  0.00230417
  0.00291610
  0.00358765
  0.00431628
  0.00593463
  0.00775569
  0.00976489
  0.01194933
  0.01429742
  0.02080739
  0.02813854
  0.03617646
  0.04482673
  0.06363014
  0.08408959
  0.10583774
  0.12859951
  0.15216623
  0.17638099
  0.22626057
  0.27760375
  0.32995349
  0.38300115
  0.43653241
  0.57147911
  0.70718625
  0.84291833
  0.97830062
  1.24728607
  1.51351565
  1.77682635
  2.03725247
  2.29489423
  2.54987043
  3.05235322
  3.54549379
  4.03000023
  4.50646156
  4.97537287
  6.11760947
  7.22029542
  8.28734912
  9.32182484
 11.30353595
 13.17999838
 14.96335835
 16.66326130
 18.28755912
 19.84292382
 22.77184827
 25.48526839
 28.01267457
 30.37776924
 32.60004307
 37.64405937
 42.07735627
 46.03152811
 49.59970748
 55.85787485
 61.19700674
 65.85324342
 69.98187149
 73.69030148
 77.05609581
 82.99736493
 88.10216461
 92.57711656
 96.56071719
100.15035659
];
if inv==0
    TE = interp1q(E_tab,E_trans,E);
else
    TE = interp1q(E_trans,E_tab,E);
end
%========================================================================
% Process Image.
%========================================================================
function I = processImage(image_name)
% Read the image and convert to grayscale
img = imread(image_name);
if size(img, 3) == 3
img = rgb2gray(img); % Convert RGB to grayscale
end
% Convert to double and scale values between 0 and 1
I = double(img) / 255;
