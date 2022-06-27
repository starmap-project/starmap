function starmap_2d_ex_gauss_filter_convergence
%STARMAP_2D_EX_GAUSS_CONVERGENCE
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: This test case uses a Gaussian as initial
%   condition and computes the PN and FPN solution for several
%   moment orders approximation and filter orders. Then, the
%   L^2 error is approximated using a highly resolved PN
%   solution. This test case demonstrates that the convergence
%   order of the FPN equations is limited by the order of the
%   filter, c.f. [Frank, Hauck, Kuepper, Convergence of filtered
%   spherical harmonic equations for radiation transport (2014)].
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
'name','Gauss Test',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',33,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'ic',@initial,... % initial conditions (defined below)
'ax',[-1 1 -1 1]*0.6,... % coordinates of computational domain
'n',[75 75],... % numbers of grid cells in each coordinate direction
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
'tfinal',0.4,... % final time (initial time is t=0)
't_plot',[],... % output times
'output',[],... % problem-specific output routine (defined below)
'f_eff',10 ... % filter effective opacity
);
orders = [3 5 9 17];                    % Orders of moment approximation.
f_order = [2 4 8 16];                             % Orders of the filter.
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
par = starmap_init(prob); % Configure data structures for starmap solver.

fprintf('Warning: This may take some minutes...\n')
fprintf('Running %s%d...\n',par.closure,par.n_mom)
solution_true = starmap_solver(par);        % Compute reference solution.
% Compute PN/FPN solutions and L2 error terms.
E = zeros(length(f_order)+1,length(orders)); R = E;
for k = 0:length(f_order)              % Loop over various filter orders.
    if k == 0, par.f_position = '';
    else
        prob.filter_order = f_order(k);     % Define order of the filter.
        prob.f_position = 'substep';            % Define filter position.
        fprintf('Running FPN with filter order %d:\n',prob.filter_order);
        prob.filterfunction = @exp_filter;      % Filter function defined
    end                                                          % below.
    for l = 1:length(orders)   % Loop over order of moment approximation.
        prob.n_mom = orders(l);
        par = starmap_init(prob); %Configure data structures for starmap solver
        fprintf('Running %s%s%d...\n',par.closure_mod,par.closure,par.n_mom);
        solution = starmap_solver(par);                     % Run solver.
        for m = 1:length(solution)
            R(k+1,l) = R(k+1,l) + ...
                sum(sum((solution(m).U-solution_true(m).U).^2));
        end
        E(k+1,l) = R(k+1,l);
        for m = length(solution)+1:length(solution_true)
            E(k+1,l) = E(k+1,l)+sum(sum((solution_true(m).U).^2));
        end
        R(k+1,l) = sqrt(R(k+1,l))/prob.n(1);   % Scaled L^2 error without
        E(k+1,l) = sqrt(E(k+1,l))/prob.n(1); % and with projection error.
    end
end
OrderE = -diff(log(E),[],2)./(ones(size(E,1),1)*diff(log(orders)));
OrderR = -diff(log(R),[],2)./(ones(size(E,1),1)*diff(log(orders)));
cname = cellstr([num2str([orders(1:end-1);orders(2:end)]','E(%d,%d)');...
    num2str([orders(1:end-1);orders(2:end)]','R(%d,%d)')]);
cformat = cellstr(repmat('bank',2*size(E,1)-2,1));
rname = cellstr(['      no filter';...
    num2str((f_order)','filter order %d')]);
text = sprintf(['%s: The terms E(N1,N2) and R(N1,N2) are the ' ...
    'convergence rates when going from moment order N1 to N2 for the ' ...
    'total and projected error.'],par.name);
output_table([OrderE,OrderR],cname,cformat',rname,text)

%========================================================================
% Problem Specific Functions
%========================================================================
function f = initial(x,y,k)
% Initial conditions (for (k-1)-st moment).
t = 3.2e-3; % pseudo-time for smoothing initial Dirac delta
f = 1/(4*pi*t)*exp(-(x.^2+y.^2)/(4*t))*(k==1);

function f = sigma_a(x,y)
% Absorption coefficient.
f = 0;

function f = sigma_s0(x,y)
% Isotropic scattering coefficient.
f = 0;

function f = exp_filter(par,dt,k)
% Filter function.
p = par.filter_order;                              % Order of the filter.
f_b = -par.f_eff/log(exp(log(eps)*(par.n_mom/(par.n_mom+1))^p));
f = exp(log(eps)*(k/(par.n_mom+1)).^p).^(dt*f_b);   % Exponential filter.

function output_table(data,cname,cformat,rname,text)
% Plot a table from data.
clf, uitable('Units','normalized','Position',[0 0.3 1 0.7],...
    'Data',data,'ColumnName',cname,'ColumnFormat',cformat,...
    'RowName',rname);
annotation('textbox',[0 0 1 0.3],'String',text)
