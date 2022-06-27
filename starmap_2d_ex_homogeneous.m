function starmap_2d_ex_homogeneous
%STARMAP_2D_EX_HOMOGENEOUS
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: Initial condition spreading outwards in domain
%   with homogeneous material coefficients. 
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
'name','Homogeneous Material',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'timestep',1,... % parabolic time step
'ic',@initial,... % initial conditions (defined below)
'ax',[-1 1 -1 1],... % coordinates of computational domain
'n',[1 1]*100,... % numbers of grid cells in each coordinate direction
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
't_plot',0.7,... % output times
'output',@output... % problem-specific output routine (defined below)
);

% Define test cases
n_mom = [1 3 5];                      % List of moment orders to compute.
sigma = {@small,@large};     % List of functions for material parameters.
top_zoom = [[-1 1]*.2,[4 5.1]*1e-3];    % Zoom window at top of solution.
save_figure = 0; % Flag whether outputs should be saved into image files.

%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
for k = 1:numel(sigma)       % Loop over different material coefficients.
    for j = n_mom                              % Loop over moment orders.
        prob.sigma_a = sigma{k}; prob.sigma_s0 = sigma{k};% Absorp.&scatter.
        for i = 1:numel(prob.t_plot), prob.figh(i) = figure; end %Open figs.
        % Reference solution
        prob.run = 1; prob.n_mom = 9;
        par = starmap_init(prob); % Configure data structures for starmap solver
        starmap_solver(par)                                 % Run solver.
        % P_N closure
        prob.run = 2; prob.n_mom = j;
        par = starmap_init(prob);     % Configure data structures for starmap solver
        starmap_solver(par)                                 % Run solver.

        if k==2&&any(j==[3 5]), axis(top_zoom), end%Activate zoom window.
        if save_figure
            filename = sprintf('fig_homogeneous_sigma%d_N%d',k,j);
            set(gcf,'position',[10 50 560 500])
            set(gcf,'paperpositionmode','auto')
            print(gcf,'-dpng',filename,'-r580')
        end
    end
end

%========================================================================
% Problem Specific Functions
%========================================================================
function f = initial(x,y,k)
% Initial conditions (for (k-1)-st moment).
f = exp(-((x.^2+y.^2)/.25^2).^2)*(k==1);

function f = small(~,~)
f = 3;

function f = large(~,~)
f = 6;

function output(par,x,y,U,step)
% Plotting routine.
width = [5,2,3,3,1.5];                                     % Line widths.
color = [.7 .7 .8;0 .2 0;1 0 0;0 0 1;.7 0 .7];             % Line colors.
style = {'-','-','--','-.','-'};                           % Line styles.
figure(par.figh(step))
if par.run>1, hold on, end
u = interp2(x,y,U',x,y*0);                   % Eval. sol. along line y=0.
plot(x,u,style{par.run},'color',color(par.run,:),...
            'linewidth',width(par.run))
title(sprintf('%s (\\sigma_a=%g, \\sigma_s=%g) at t = %0.2f',...
    par.name,par.sigma_a(0,0,0),par.sigma_s0(0,0,0),par.t_plot(step)))
xlabel('x'), ylabel('\psi_0')
[~,~,~,text] = legend; text{end+1} = sprintf('%s%d',par.closure,par.n_mom); 
legend(text)
axis tight, set(gca,'xlim',[-1 1])
ax = get(gca,'ylim'); set(gca,'ylim',ax+[-1 1]*diff(ax)*.02)
drawnow
