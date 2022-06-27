function starmap_2d_create_mms
%STARMAP_2D_CREATE_MMS
%   Creates the file STARMAP_2D_EX_MMS_AUTO.M, which is an
%   example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometry.
%
%   Uses the method of manufactured solutions (MMS) to
%   generate test cases for the PN or the SPN equations, for
%   any moment order N, by prescribing a solution and then
%   symbolically computing the source that generates this
%   solution.
%
%   Note: Running this file requires the Symbolic Toolbox.
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
'name','Manufactured Solution Test',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',3,... % order of moment approximation
'tfinal',1,... % final time (initial time is t=0)
'source',@source,... % source term (defined below)
'ic',@initial,... % initial condition
'ax',[0 1 0 1],... % coordinates of computational domain
'bc',[1 1],... % type of boundary cond. (0 = periodic, 1 = extrapolation)
...                % Careful: manufactured solution must satisfy the b.c.
't_plot',[],... % output times (none here)
'output',[]... % output routine (none here)
);

points = 2.^(4:7);      % Set number of grid points for convergence test.

%========================================================================
% Manufacture Solution
%========================================================================
syms t x y solution
% Assign variables
n_mom = prob.n_mom;
switch prob.closure
    case 'P',  [mom_order,~,Mx,My] = feval('starmap_closure_pn',n_mom,2);
    case 'SP', [mom_order,~,Mx,My] = feval('starmap_closure_spn',n_mom,2);
end

n_sys = length(Mx);
u = zeros(n_sys,1); scatterm = u;
u = sym(u); scatterm = sym(scatterm);
absorb = sigma_a(x,y,t);
scatter0 = sigma_sm(x,y,t,0);
for k = 1:n_sys
    u(k) = solution(x,y,t,k);
    scatterm(k) = sigma_sm(x,y,t,mom_order(k));
end

% Compute source
q = simplify(diff(u,t) + Mx*diff(u,x) + My*diff(u,y) + ...
    (absorb+scatter0)*u - scatterm.*u);

%========================================================================
% Write Example File
%========================================================================
function_name = 'starmap_2d_ex_mms_auto';
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
fprintf(fid,'%s\n','''name'',''Manufactured Solution Test'',... % name of example');
fprintf(fid,'%s\n',['''closure'',''',prob.closure,''',... % type of closure (can be ''P'' or ''SP'')']);
fprintf(fid,'%s\n',['''n_mom'',',num2str(prob.n_mom),',... % order of moment approximation']);
fprintf(fid,'%s\n',['''tfinal'',',num2str(prob.tfinal),',... % final time (initial time is t=0)']);
fprintf(fid,'%s\n','''sigma_a'',@sigma_a,... % absorption coefficient (defined below)');
fprintf(fid,'%s\n','''sigma_s0'',@sigma_s0,... % isotropic scattering coefficient (def. below)');
fprintf(fid,'%s\n','''sigma_sm'',@sigma_sm,... % aniso. scattering coefficient (defined below)');
fprintf(fid,'%s\n','''source'',@source,... % source term (defined below)');
fprintf(fid,'%s\n','''ic'',@initial,... % initial condition');
fprintf(fid,'%s\n',['''ax'',[',num2str(prob.ax,'%g '),'],... % coordinates of computational domain']);
fprintf(fid,'%s\n',['''bc'',[',num2str(prob.bc,'%g '), '],... % type of boundary cond. (0 = periodic, 1 = extrapolation)']);
fprintf(fid,'%s\n','''t_plot'',[],... % output times (none here)');
fprintf(fid,'%s\n','''output'',[]... % output routine (none here)');
fprintf(fid,'%s\n',');');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Moment System Setup and Solver Execution');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','par = starmap_init(prob);     % Configure data structures for starmap solver');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Run solver');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n',['points = [',num2str(points,'%d '), '];']);
fprintf(fid,'%s\n','E_1 = zeros(par.n_mom+1,length(points)); E_2 = E_1; E_inf = E_1;');
fprintf(fid,'%s\n','for k = 1:length(points)            % Loop over various grid resolutions.');
fprintf(fid,'%s\n','    prob.n = [1 1]*points(k);                     % Numbers of grid cells.');
fprintf(fid,'%s\n','    par = starmap_init(prob);     % Configure data structures for starmap solver');
fprintf(fid,'%s\n','    sol = starmap_solver(par);                              % Run solver.');
fprintf(fid,'%s\n','    for j = 1:length(sol)                % Loop over solution components.');
fprintf(fid,'%s\n','        [X,Y] = ndgrid(sol(j).x,sol(j).y);% Grid on which solution lives.');
fprintf(fid,'%s\n','        Utrue = solution(X,Y,j,par.tfinal);     % Evaluate true solution.');
fprintf(fid,'%s\n','        D = abs(sol(j).U-Utrue);  % Difference between num. and true sol.');
fprintf(fid,'%s\n','        ind = par.mom_order(j)+1;                    % Moment order index');
fprintf(fid,'%s\n','        E_1(ind,k) = E_1(ind,k)+sum(sum(D))/points(k)^2;%Scaled L1 error.');
fprintf(fid,'%s\n','        E_2(ind,k) = E_2(ind,k)+sum(sum(D.^2))/points(k)^2;%Sc. L2 error.');
fprintf(fid,'%s\n','        E_inf(ind,k) = max(E_inf(ind,k),max(max(D)));        % Max error.');
fprintf(fid,'%s\n','    end');
fprintf(fid,'%s\n','end');
fprintf(fid,'%s\n','E_2 = sqrt(E_2);');
fprintf(fid,'%s\n','h = 1./points; axh = [min(h) max(h)];                        % Axes range');
fprintf(fid,'%s\n','f = (axh(2)/axh(1))^.15; axh = axh.*[1/f f];         % for resolutions h.');
fprintf(fid,'%s\n','E = [E_inf;E_1;E_2];                               % All errors combined.');
fprintf(fid,'%s\n','cmin = min(min(E).*points.^2)/1.5;    % Lower bounding line with slope 2.');
fprintf(fid,'%s\n','cmax = max(max(E).*points.^2)*1.5;    % Upper bounding line with slope 2.');
fprintf(fid,'%s\n','axE = axh''.^2*[cmin,cmax];                       % Error values of lines.');
fprintf(fid,'%s\n','cols = ceil(sqrt(par.n_mom+1));              % Number of subplot columns.');
fprintf(fid,'%s\n','rows = ceil((par.n_mom+1)/cols);                % Number of subplot rows.');
fprintf(fid,'%s\n','clf');
fprintf(fid,'%s\n','for j = 1:par.n_mom+1                          % Loop over moment orders.');
fprintf(fid,'%s\n','    subplot(rows,cols,j)');
fprintf(fid,'%s\n','    loglog(h,E_1(j,:),''.'',h,E_2(j,:),''x'',h,E_inf(j,:),''o'',...%Plot errors');
fprintf(fid,'%s\n','        axh,axE,''k-'',''Linewidth'',1)     % and reference lines of slope 2.');
fprintf(fid,'%s\n','    axis([axh axE([1 4])])');
fprintf(fid,'%s\n','    xlabel(''resolution h'')');
fprintf(fid,'%s\n','    ylabel(sprintf(''error in moments of order %d'',j-1))');
fprintf(fid,'%s\n','    set(gca,''xtick'',fliplr(h),''xticklabel'',...    % Create axes labels on');
fprintf(fid,'%s\n','        cellfun(@(s)cat(2,''1/'',num2str(s)),...  % h-axis that are inverse');
fprintf(fid,'%s\n','        num2cell(fliplr(points)),''UniformOutput'',false))   % powers of 2.');
fprintf(fid,'%s\n','    legend(''L^1 error'',''L^2 error'',''L^\infty error'',''slope 2'',''Location'',''NorthWest'')');
fprintf(fid,'%s\n','    title(sprintf(''%s with %s%d'',par.name,par.closure,par.n_mom))');
fprintf(fid,'%s\n','end');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','% Problem Specific Functions');
fprintf(fid,'%s\n','%========================================================================');
fprintf(fid,'%s\n','function f = initial(x,y,k)');
fprintf(fid,'%s\n','switch k');
for k = 1:n_sys
    fh_u = matlabFunction(u(k),'vars',[x,y,t]);
    fprintf(fid,'%s\n',['case ',num2str(k),', f = feval(',func2str(fh_u),',x,y,0);']);
end
fprintf(fid,'%s\n','end');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = source(x,y,t,k)');
fprintf(fid,'%s\n','switch k');
for k = 1:n_sys
    fh_q = matlabFunction(q(k),'vars',[x,y,t]);
    fprintf(fid,'%s\n',['case ',num2str(k),', f = feval(',func2str(fh_q),',x,y,t);']);
end
fprintf(fid,'%s\n','end');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = solution(x,y,k,t)');
fprintf(fid,'%s\n','switch k');
for k = 1:n_sys
    fh_u = matlabFunction(u(k),'vars',[x,y,t]);
    fprintf(fid,'%s\n',['case ',num2str(k),', f = feval(',func2str(fh_u),',x,y,t);']);
end
fprintf(fid,'%s\n','end');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = sigma_a(x,y,t)');
fh_a = matlabFunction(absorb,'vars',[x,y,t]);
fprintf(fid,'%s\n',['f = feval(',func2str(fh_a),',x,y,t);']);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = sigma_s0(x,y,t)');
fh_s0 = matlabFunction(scatter0,'vars',[x,y,t]);
fprintf(fid,'%s\n',['f = feval(',func2str(fh_s0),',x,y,t);']);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','function f = sigma_sm(x,y,m,t)');
fprintf(fid,'%s\n','switch m');
for m = 0:n_mom
    scatterm(1) = sigma_sm(x,y,t,m);
    fh_sm = matlabFunction(scatterm(1),'vars',[x,y,t]);
    fprintf(fid,'%s\n',['case ',num2str(m),...
        ', f = feval(',func2str(fh_sm),',x,y,t);']);
end
fprintf(fid,'%s\n','end');
fclose(fid);
fprintf(' Done.\n')

%========================================================================
% Problem Specific Functions
%========================================================================
function f = sigma_a(x,y,t)
% Absorption coefficient.
syms f x y t
f = cos(2*pi*y)*t;

function f = sigma_sm(x,y,t,m)
% Moments of scattering kernel.
syms f x y t
g = 0.9;
f = g^m;

function u = solution(x,y,t,k)
% Manufactured solution
syms u x y t
switch k
    case 1, u = exp(-t)*(sin(2*pi*(x)))^2;
    otherwise, u = 0;
end
