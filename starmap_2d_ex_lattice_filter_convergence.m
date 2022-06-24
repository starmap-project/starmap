function starmap_2d_ex_lattice_filter_convergence
%STARMAP_2D_EX_LATTICE_FILTER_CONVERGENCE
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Example: Lattice/Checkerboard geometry, as defined in
%   [Brunner, Holloway, JCP 210 (2005) 386-399].
%   Computation up to time t=2.8.
%   PN and FPN solutions are computed for several moment order 
%   approximations and filter orders. The L^2 error is
%   approximated using a highly resolved PN solution. This
%   test case demonstrates that the convergence order of the
%   PN and FPN equations is limited by the regularity of the
%   solution in the spatial and angular variable, c.f.
%   [Frank, Hauck, Kuepper, Convergence of filtered spherical
%   harmonic equations for radiation transport (2014)].
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
'name','Lattice Test',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',33,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'source',@source,... % source term (defined below)
'ax',[0 7 0 7],... % coordinates of computational domain
'n',[125 125],... % numbers of grid cells in each coordinate direction
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
'tfinal',2.8,... % final time (initial time is t=0)
't_plot',[],... % output times
'output',[],... % problem-specific output routine
'f_eff',10 ... % filter effective opacity
);
orders = [3 5 9 17];                    % Orders of moment approximation.
f_order = [2 4 8 16];                             % Orders of the filter.
%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
par = starmap_init(prob);  % Configure data structures for starmap solver

fprintf('Warning: This may take some minutes...\n')
fprintf('Running %s%d...\n',par.closure,par.n_mom)
solution_true = starmap_solver(par);    % Compute reference solution.

% Compute estimates from reference solution.
n_sys = (prob.n_mom+1)*(prob.n_mom+2)/2;     % Number of system components.
moment_order = ceil(sqrt(2*(1:n_sys)+1/4)-3/2);           % Moment order.
Moments = zeros(1,prob.n_mom+1);           % Moments of the same order and
DiffXMoments = Moments;                               % spatial gradient.
for i = 1:length(solution_true)
    index = moment_order(i)+1;
    Moments(index) = Moments(index)+sum(sum((solution_true(i).U).^2));
    DiffXMoments(index) = DiffXMoments(index) +...
        sum(sum(diff(solution_true(i).U,[],1).^2)) +...
        sum(sum(diff(solution_true(i).U,[],2).^2));
end
Moments = sqrt(Moments)/prob.n(1);
DiffXMoments = sqrt(DiffXMoments)/prob.n(1);
output_moments(par,Moments,DiffXMoments)                % Plot estimates.
OrderMoments = -diff(log(Moments([orders;orders+1])),[],2)./...
    diff(log([orders-1;orders]),[],2);
OrderDiffXMoments = -diff(log(DiffXMoments([orders;orders+1])),[],2)./...
    diff(log([orders-1;orders]),[],2);
cname = {'B(N1,N2)|even','D(N1,N2)|even','B(N1,N2)|odd','D(N1,N2)|odd'};
rname = cellstr(num2str([orders(1:end-1)-1;orders(2:end)-1;...
    orders(1:end-1);orders(2:end)]','(N1,N2) = (%d,%d) and (%d,%d)'));
text = sprintf('%s: The terms B(N1,N2) and D(N1,N2) are the decay rates of the moments and their differentials when going from N1 to N2 (computed with P%d).',par.name,par.n_mom);
figure
output_table([OrderMoments',OrderDiffXMoments'],cname,[],rname,text)
% Compute PN/FPN solutions and L2 error terms.
E = zeros(length(f_order)+1,length(orders)); R = E;
for k = 0:length(f_order)              % Loop over various filter orders.  
    if k == 0, prob.f_position = '';
    else, prob.filter_order = f_order(k);     % Define order of the filter.
        prob.f_position = 'substep';             % Define filter position.
        fprintf('Running FPN with filter order %d:\n',prob.filter_order);
        prob.filterfunction = @exp_filter;       % Filter function defined
    end                                                          % below.
    for l = 1:length(orders)   % Loop over order of moment approximation.
        prob.n_mom = orders(l);
        par = starmap_init(prob); % Configure data structures for starmap solver
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
        R(k+1,l) = sqrt(R(k+1,l))/par.n(1);    % Scaled L^2 error without
        E(k+1,l) = sqrt(E(k+1,l))/par.n(1);  % and with projection error.
    end
end
OrderE = -diff(log(E),[],2)./(ones(size(E,1),1)*diff(log(orders)));
OrderR = -diff(log(R),[],2)./(ones(size(E,1),1)*diff(log(orders)));
cname = cellstr([num2str([orders(1:end-1);orders(2:end)]','E(%d,%d)');...
    num2str([orders(1:end-1);orders(2:end)]','R(%d,%d)')]);
cformat = cellstr(repmat('bank',2*size(E,1)-2,1));
rname = cellstr(['      no filter';...
    num2str((f_order)','filter order %d')]);
text = sprintf('%s: The terms E(N1,N2) and R(N1,N2) are the convergence rates when going from moment order N1 to N2 for the total and projected error.',par.name);
figure, output_table([OrderE,OrderR],cname,cformat',rname,text)

%========================================================================
% Problem Specific Functions
%========================================================================
function f = sigma_a(x,y)
% Absorption coefficient.
cx = ceil(x); cy = ceil(y);
g = (ceil((cx+cy)/2)*2==(cx+cy)).*(1<cx&cx<7&1<cy&cy-2*abs(cx-4)<4);
f = (1-g)*0+g*10;

function f = sigma_s0(x,y)
% Isotropic scattering coefficient.
cx = ceil(x); cy = ceil(y);
g = (ceil((cx+cy)/2)*2==(cx+cy)).*(1<cx&cx<7&1<cy&cy-2*abs(cx-4)<4);
f = (1-g)*1+g*0;

function f = source(x,y)
% Radiation source (only for zeroth moment).
f = 3<x&x<4&3<y&y<4;

function f = exp_filter(par,dt,k)
% Filter function; exponential filter.
p = par.filter_order;                              % Order of the filter.
f_b =  -par.f_eff/log(exp(log(eps)*(par.n_mom/(par.n_mom+1))^p));
f = exp(log(eps)*(k/(par.n_mom+1)).^p).^(dt*f_b);

function output_table(data,cname,cformat,rname,text)
% Plot table from data.
clf, uitable('Units','normalized','Position',[0 0.3 1 0.7],...
    'Data',data,'ColumnName',cname,'ColumnFormat',cformat,...
    'RowName',rname);
annotation('textbox',[0 0 1 0.3],'String',text)

function output_moments(par,B,D)
% Plot sequences B and D and indicate even/odd sequence elements.
figure, subplot(1,2,1)
loglog(1:2:par.n_mom,B(2:2:par.n_mom+1),'*','Color',[0 0 1])
hold on, loglog(2:2:par.n_mom,B(3:2:par.n_mom+1),'x','Color',[1 0 0])
legend({'B_N (N odd)','B_N (N even)'},'Location','southwest')
xlabel({'degree N'}), xlim([1 par.n_mom])
ylabel({'B_N (L2 norm of the moments of degree N)'}), ylim([1e-6 1e-1])
subplot(1,2,2)
loglog(1:2:par.n_mom,D(2:2:par.n_mom+1),'*','Color',[0 0 1])
hold on, loglog(2:2:par.n_mom,D(3:2:par.n_mom+1),'x','Color',[1 0 0])
legend({'D_N (N odd)','D_N (N even)'},'Location','southwest')
xlabel({'degree N'}), xlim([1 par.n_mom])
ylabel({'D_N (L2 norm of the differentials of the moments of degree N)'})
ylim([1e-6 1e-1])
