function starmap_2d_ex_l2norm
%STARMAP_2D_EX_L2NORM
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: Expanding Gaussian in void. The temporal
%   evolution of the L2 norm of the numerical solution is
%   plotted, to test how well numerical scheme conserves the
%   discrete L2 norm.
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
'name','Expanding Gaussian Test',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',5,... % order of moment approximation
'ic',@initial,... % initial conditions (defined below)
'ax',[-1 1 -1 1],... % coordinates of computational domain
'n',[1 1]*50,... % numbers of grid cells in each coordinate direction
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',linspace(0,0.5,25),... % output times
'output',@output,... % problem-specific output routine (defined below)
'mom_output','all'... % output all moments to plotting routine
);
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
par = starmap_init(prob); % Configure data structures for starmap solver.
clf
for k = 1:3             % Run test case with three different resolutions.
    starmap_solver(par);                                    % Run solver.
    prob.n = prob.n*2;                               % Double resolution.
    par = starmap_init(prob); %Configure data structures for starmap solver.
end

%========================================================================
% Problem Specific Functions
%========================================================================
function f = initial(x,y,k)
% Initial conditions (for (k-1)-st moment).
t = 1e-2;                % Pseudo-time for smoothing initial Dirac delta.
f = 1/(4*pi*t)*exp(-(x.^2+y.^2)/(4*t))*(k==1);                % Gaussian.

function par = output(par,x,y,U,step)
% Plotting routine, showing the 2D results and time-evolution of
% discrete energy.
t = par.t_plot(step);                                   % Time of output.
phi_max = max(max(U{1}))+eps;                       % Axis/color scaling.
h = [x{1}(2)-x{1}(1),y{1}(2)-y{1}(1)];                 % Grid resolution.
sol_norm = sqrt(sum(cell2mat(cellfun(@(u)L2(u,par.prob.n,h),U,'Un',0))));
if step==1, par.sol_norm = sol_norm;   % Create vector of solution norms.
else, par.sol_norm(step) = sol_norm; end%Update vector of solution norms.
subplot(1,2,1)
imagesc(x{1},y{1},U{1}')                      % 2D plot of approximation.
axis xy equal tight, caxis([-1 1]*phi_max)
title(sprintf('%s with %s%d at t = %0.2f',par.name,par.closure,par.n_mom,t))
colormap jet(255)
subplot(1,2,2)
title('Evolution of Discrete Energy')
xlabel('t'), ylabel('Relative L^2 norm of numerical solution')
h_graph = findobj(gca,'Marker','.');                  % Handle of graphs.
h_text = findobj(gca,'Type','Text');                   % Handle of texts.
if step>1                                      % After first step, remove
    delete(h_graph(1)), delete(h_text(1))      % previous graph and text.
end
n_graphs = length(findobj(gca));                      % Number of graphs.
c = mod(n_graphs*5,7)/6;                        % Assign color for graph.
hold on
plot(par.t_plot(1:step),par.sol_norm/par.sol_norm(1),'.-',...% Plot time-
   'Color',[c 0 1-c],'MarkerSize',12,'LineWidth',1)%evolution of L2 norm.
text(par.t_plot(step),par.sol_norm(step)/par.sol_norm(1),...  % Plot mesh
    sprintf(' h = %g',h(1)),'Color',[c 0 1-c])                    % size.
hold off, set(gca,'XLim',[0 par.t_plot(end-1)])
drawnow

function v = L2(U,n,h)
% Compute squared L2 norm of a solution component, considering the
% staggered grid nature.
s = size(U)==n;                   % Check which components are staggered.
w1 = ones(1,n(1)+1-s(1)); w1([1 end]) = s(1)/2; % Weights in x-direction.
w2 = ones(1,n(2)+1-s(2)); w2([1 end]) = s(2)/2; % Weights in y-direction.
v = prod(h)*w1*U.^2*w2';                % Square of L2 norm of component.
