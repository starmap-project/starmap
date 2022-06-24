function starmap_2d_create_hemisphere_filter_convergence
%STARMAP_2D_CREATE_HEMISPHERE_CONVERGENCE
%   Creates the file STARMAP_2D_EX_HEMISPHERE_CONVERGENCE_AUTO.M,
%   which is an example case for STARMAP_SOLVER, a second
%   order staggered grid finite difference solver for linear
%   hyperbolic moment approximations to radiative transfer
%   in 1D-3D geometry.
%
%   An initial characteristic function in angle, which is a
%   hemisphere pointing in a certain direction, and a spatial
%   distribution can be specified. The routine generates an
%   initial condition and writes default functions for all
%   other problem parameters that can later be changed in the
%   generated example file. This routine is specific to the PN
%   and FPN equations. In the example file PN and FPN solutions
%   are computed for several moment order approximations and
%   filter orders. The L^2 error is approximated using a highly
%   resolved PN solution. This test case demonstrates that the
%   convergence order of the PN and FPN equations is limited by
%   the regularity of the solution in the angular variable,
%   c.f. [Frank, Hauck, Kuepper, Convergence of filtered
%   spherical harmonic equations for radiation transport (2014)].
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
'name','Hemisphere Test',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',33,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'ax',[-0.6  0.6 -0.6  0.6],... % coordinates of computational domain
'n',[75 75],... % numbers of grid cells in each coordinate direction
'bc',[0 0],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
'tfinal',0.3,... % final time (initial time is t=0)
't_plot',[],... % output times
'output',[],... % problem-specific output routine
'f_eff',10 ... % filter effective opacity
);
orders = [3 5 9 17];                    % Orders of moment approximation.
f_order = [2 4 8 16];                             % Orders of the filter.
% The source is a characteristic function on a hemisphere in angle,
% times a Gaussian in space.
phi_beam = 0;                         % Main Angle of beam w.r.t. x-axis.
space_beam = @(x,y) 1/(4*pi*3.2e-3)*...            % Spatial distribution
    exp(-(x.^2+y.^2)/(4*3.2e-3));                              % of beam.
%========================================================================
% Compute Initial Condition
%========================================================================
mu_beam = 0;
n_mom = prob.n_mom;
n_sys = (n_mom+1)*(n_mom+2)/2;

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
    k = -l;
    for index = l*(l+1)/2+1:(l+1)*(l+2)/2
        u(index) = CoeffHemisphere(l,k,mu_beam,phi_beam);
        k = k+2;
    end
end

% Transform to StarMAP variables
StarMAPmoments = M*u;

%========================================================================
% Write Example File
%========================================================================
function_name = 'starmap_2d_ex_hemisphere_filter_convergence_auto';
fprintf('Writing example file %s.',[function_name,'.m']);
fid = fopen([function_name,'.m'],'w');
fprintf(fid,'function %s\n',function_name);
fprintf(fid,'%s%s\n','%',upper(function_name));
fprintf(fid,'%s\n','%   Example case for STARMAP_SOLVER, a second order staggered');
fprintf(fid,'%s\n','%   grid finite difference solver for linear hyperbolic moment');
fprintf(fid,'%s\n','%   approximations to radiative transfer in 1D-3D geometry.');
fprintf(fid,'%s\n','%');
fprintf(fid,'%s%s%s\n','%   Created by the file ',mfilename,'.m');
fprintf(fid,'%s\n','%');
fprintf(fid,'%s\n','%   Version 2.0');
fprintf(fid,'%s\n','%   Copyright (c) 06/28/2022 Benjamin Seibold, Martin Frank, and');
fprintf(fid,'%s\n','%                            Rujeko Chinomona');
fprintf(fid,'%s\n','%   http://www.math.temple.edu/~seibold');
fprintf(fid,'%s\n','%   https://www.scc.kit.edu/personen/martin.frank.php');
fprintf(fid,'%s\n','%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0)');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%   For license, see file starmap_solver.m, as published on');
fprintf(fid,'%s\n','%   https://github.com/starmap-project/starmap');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Problem Parameters');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','prob = struct(...');
fprintf(fid,'%s\n','''name'',''Hemisphere Test'',... % name of example');
fprintf(fid,'%s\n','''closure'',''P'',... % type of closure (can be ''P'' or ''SP'')');
fprintf(fid,'%s\n',['''n_mom'',',num2str(prob.n_mom,'%g'),',... % order of moment approximation']);
fprintf(fid,'%s\n','''sigma_a'',@sigma_a,... % absorption coefficient (defined below)');
fprintf(fid,'%s\n','''sigma_s0'',@sigma_s0,... % isotropic scattering coefficient (def. below)');
fprintf(fid,'%s\n','''source'',@source,... % source term (defined below)');
fprintf(fid,'%s\n',['''ax'',[',num2str(prob.ax,'%g '),'],... % coordinates of computational domain']);
fprintf(fid,'%s\n',['''n'',[',num2str(prob.n,'%g '),'],... % numbers of grid cells in each coordinate direction']);
fprintf(fid,'%s\n',['''bc'',[',num2str(prob.bc,'%g '), '],... % type of boundary cond. (0 = periodic, 1 = extrapolation)']);
fprintf(fid,'%s\n',['''tfinal'',',num2str(prob.tfinal),',... % final time (initial time is t=0)']);
fprintf(fid,'%s\n',['''t_plot'',[',num2str(prob.t_plot,'%g '),'],... % output times']);
fprintf(fid,'%s\n','''output'',[],... % output routine');
fprintf(fid,'%s\n',['''f_eff'',',num2str(prob.f_eff),' ... % filter effective opacity']);
fprintf(fid,'%s\n',');');
fprintf(fid,'%s\n',['orders = [',num2str(orders,'%g '),'];                  % Orders of moment approximation.']);
fprintf(fid,'%s\n',['f_order = [',num2str(f_order,'%g '),'];                            % Order of the filter.']);
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Moment System Setup and Solver Execution');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','par = starmap_init(prob);     % Configure data structures for starmap solver');
fprintf(fid,'%s\n','fprintf(''Warning: This may take some minutes...\n'')');
fprintf(fid,'%s\n','fprintf(''Running %s%d...\n'',par.closure,par.n_mom)');
fprintf(fid,'%s\n','solution_true = starmap_solver(par);        % Compute reference solution.');
fprintf(fid,'%s\n','% Compute estimates from reference solution.');
fprintf(fid,'%s\n','n_sys = (par.n_mom+1)*(par.n_mom+2)/2;     % Number of system components.');
fprintf(fid,'%s\n','moment_order = ceil(sqrt(2*(1:n_sys)+1/4)-3/2);           % Moment order.');
fprintf(fid,'%s\n','Moments = zeros(1,par.n_mom+1);           % Moments of the same order and');
fprintf(fid,'%s\n','DiffXMoments = Moments;                               % spatial gradient.');
fprintf(fid,'%s\n','for i = 1:length(solution_true)');
fprintf(fid,'%s\n','    index = moment_order(i)+1;');
fprintf(fid,'%s\n','    Moments(index) = Moments(index)+sum(sum((solution_true(i).U).^2));');
fprintf(fid,'%s\n','    DiffXMoments(index) = DiffXMoments(index) +...');
fprintf(fid,'%s\n','        sum(sum(diff(solution_true(i).U,[],1).^2)) +...');
fprintf(fid,'%s\n','        sum(sum(diff(solution_true(i).U,[],2).^2));');
fprintf(fid,'%s\n','end');
fprintf(fid,'%s\n','Moments = sqrt(Moments)/prob.n(1);');
fprintf(fid,'%s\n','DiffXMoments = sqrt(DiffXMoments)/prob.n(1);');
fprintf(fid,'%s\n','output_moments(par,Moments,DiffXMoments)                % Plot estimates.');
fprintf(fid,'%s\n','OrderMoments = -diff(log(Moments([orders;orders+1])),[],2)./...');
fprintf(fid,'%s\n','    diff(log([orders-1;orders]),[],2);');
fprintf(fid,'%s\n','OrderDiffXMoments = -diff(log(DiffXMoments([orders;orders+1])),[],2)./...');
fprintf(fid,'%s\n','    diff(log([orders-1;orders]),[],2);');
fprintf(fid,'%s\n','cname = {''B(N1,N2)|even'',''B(N1,N2)|odd'',''D(N1,N2)|even'',''D(N1,N2)|odd''};');
fprintf(fid,'%s\n','rname = cellstr(num2str([orders(1:end-1)-1;orders(2:end)-1;...');
fprintf(fid,'%s\n','    orders(1:end-1);orders(2:end)]'',''(N1,N2) = (%d,%d) and (%d,%d)''));');
fprintf(fid,'%s\n','text = sprintf(''%s: The terms B(N1,N2) and D(N1,N2) are the decay rates of the moments and their differentials when going from N1 to N2 (computed with P%d).'',par.name,par.n_mom);');
fprintf(fid,'%s\n','figure, output_table([OrderMoments'',OrderDiffXMoments''],cname,[],rname,text)');
fprintf(fid,'%s\n','% Compute PN/FPN solutions and L2 error terms.');
fprintf(fid,'%s\n','E = zeros(length(f_order)+1,length(orders)); R = E;');
fprintf(fid,'%s\n','for k = 0:length(f_order)              % Loop over various filter orders.');
fprintf(fid,'%s\n','    if k == 0, prob.f_position = '''';');
fprintf(fid,'%s\n','    else prob.filter_order = f_order(k);     % Define order of the filter.');
fprintf(fid,'%s\n','        prob.f_position = ''substep'';             % Define filter position.');
fprintf(fid,'%s\n','        fprintf(''Running FPN with filter order %d:\n'',prob.filter_order);');
fprintf(fid,'%s\n','        prob.filterfunction = @exp_filter;       % Filter function defined');
fprintf(fid,'%s\n','    end                                                          % below.');
fprintf(fid,'%s\n','    for l = 1:length(orders)   % Loop over order of moment approximation.');
fprintf(fid,'%s\n','        prob.n_mom = orders(l);');
fprintf(fid,'%s\n','        par = starmap_init(prob);     % Configure data structures for starmap solver');
fprintf(fid,'%s\n','        fprintf(''Running %s%d...\n'',par.closure,par.n_mom);');
fprintf(fid,'%s\n','        solution = starmap_solver(par);                     % Run solver.');
fprintf(fid,'%s\n','        for m = 1:length(solution)');
fprintf(fid,'%s\n','            R(k+1,l) = R(k+1,l) + ...');
fprintf(fid,'%s\n','                sum(sum((solution(m).U-solution_true(m).U).^2));');
fprintf(fid,'%s\n','        end');
fprintf(fid,'%s\n','        E(k+1,l) = R(k+1,l);');
fprintf(fid,'%s\n','        for m = length(solution)+1:length(solution_true)');
fprintf(fid,'%s\n','            E(k+1,l) = E(k+1,l)+sum(sum((solution_true(m).U).^2));');
fprintf(fid,'%s\n','        end');
fprintf(fid,'%s\n','        R(k+1,l) = sqrt(R(k+1,l))/par.n(1);    % Scaled L^2 error without');
fprintf(fid,'%s\n','        E(k+1,l) = sqrt(E(k+1,l))/par.n(1);  % and with projection error.');
fprintf(fid,'%s\n','    end');
fprintf(fid,'%s\n','end');
fprintf(fid,'%s\n','OrderE = -diff(log(E),[],2)./(ones(size(E,1),1)*diff(log(orders)));');
fprintf(fid,'%s\n','OrderR = -diff(log(R),[],2)./(ones(size(E,1),1)*diff(log(orders)));');
fprintf(fid,'%s\n','cname = cellstr([num2str([orders(1:end-1);orders(2:end)]'',''E(%d,%d)'');...');
fprintf(fid,'%s\n','    num2str([orders(1:end-1);orders(2:end)]'',''R(%d,%d)'')]);');
fprintf(fid,'%s\n','cformat = cellstr(repmat(''bank'',2*size(E,1)-2,1));');
fprintf(fid,'%s\n','rname = cellstr([''      no filter'';...');
fprintf(fid,'%s\n','    num2str((f_order)'',''filter order %d'')]);');
fprintf(fid,'%s\n','text = sprintf(''%s: The terms E(N1,N2) and R(N1,N2) are the convergence rates when going from moment order N1 to N2 for the total and projected error.'',par.name);');
fprintf(fid,'%s\n','figure, output_table([OrderE,OrderR],cname,cformat'',rname,text)');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Problem Specific Functions');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','function f = sigma_a(x,y)');
fprintf(fid,'%s\n','% Absorption coefficient.');
fprintf(fid,'%s\n','f = 0;');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = sigma_s0(x,y)');
fprintf(fid,'%s\n','% Total scattering coefficient.');
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
fprintf(fid,'%s\n','function f = exp_filter(par,dt,k)');
fprintf(fid,'%s\n','% Filter function; exponential filter.');
fprintf(fid,'%s\n','p = par.filter_order;                              % Order of the filter.');
fprintf(fid,'%s\n','f_b =  -par.f_eff/log(exp(log(eps)*(par.n_mom/(par.n_mom+1))^p));');
fprintf(fid,'%s\n','f = exp(log(eps)*(k/(par.n_mom+1)).^p).^(dt*f_b);');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function output_table(data,cname,cformat,rname,text)');
fprintf(fid,'%s\n','% Plot table from data.');
fprintf(fid,'%s\n','clf, uitable(''Units'',''normalized'',''Position'',[0 0.3 1 0.7],...');
fprintf(fid,'%s\n','    ''Data'',data,''ColumnName'',cname,''ColumnFormat'',cformat,...');
fprintf(fid,'%s\n','    ''RowName'',rname);');
fprintf(fid,'%s\n','annotation(''textbox'',[0 0 1 0.3],''String'',text)');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function output_moments(par,B,D)');
fprintf(fid,'%s\n','% Plot sequences B and D and indicate even/odd sequence elements.');
fprintf(fid,'%s\n','figure, loglog(1:2:par.n_mom,B(2:2:par.n_mom+1),''*'',''Color'',[0 0 1])');
fprintf(fid,'%s\n','hold on, loglog(2:2:par.n_mom,B(3:2:par.n_mom+1),''o'',''Color'',[0 0 0.8])');
fprintf(fid,'%s\n','loglog(1:2:par.n_mom,D(2:2:par.n_mom+1),''+'',''Color'',[1 0 0])');
fprintf(fid,'%s\n','loglog(2:2:par.n_mom,D(3:2:par.n_mom+1),''x'',''Color'',[0.8 0 0])');
fprintf(fid,'%s\n','legend({''B_N (N odd)'',''B_N (N even)'',''D_N (N odd)'',''D_N (N even)''},''Location'',''SouthWest'')');
fprintf(fid,'%s\n','xlabel({''degree N''}), xlim([1 par.n_mom])');
fprintf(fid,'%s\n','ylabel({''B_N and D_N (L2 norm of the moments and differentials of degree N)''})');
fclose(fid);
fprintf(' Done.\n')

%========================================================================
% Functions
%========================================================================
function y = CoeffHemisphere(l,k,mu,phi)
% Expansion coefficients for a hemisphere in angle pointing in
% x-direction.
if l==0, y = sqrt(pi);                        % Define zeros coefficient.
elseif mod(l,2) && mod(k,2)   % Compute all other non-zero coefficients.
    % Compute integral from 0 to 1 of the legendre polynomial of order l.
    z1 = legendre(l-1,0);
    z2 = legendre(l+1,0);
    c1 = (z1(1)-z2(1))/(2*l+1);
    % Compute normalization factor of the spherical harmonics.
    c2 = 1;
    ak = abs(k);
    for i = ak:-1:-ak+1
        c2 = c2*sqrt(l+i);
    end
    c2 = 1/c2*sqrt((2*l+1)/(4*pi));
    % Compute associated Legendre polynomial at 0.
    z3 = legendre(l,mu)*exp(1i*k*phi);
    c3 = z3(ak+1);
    % Compute Expansion coefficient from c1, c2, c3.
    y = 2*pi*(-1)^min(0,k)*c1*c2*c3;
else                                % Set all other coefficients to zero.
    y=0;
end
