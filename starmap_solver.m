function output = starmap_solver(par)
%STARMAP_SOLVER
%   A second order staggered grid finite difference solver for
%   linear hyperbolic moment approximations (such as P_N and
%   SP_N) to the equations of radiative transfer in 1D-3D
%   geometry. Needs to be called from an example file by
%   STARMAP_SOLVER(PAR), where PAR defines the problem
%   parameters and system matrices.
%   OUTPUT = STARMAP_SOLVER(PAR) writes the final state of the
%   zeroth moment and its spatial coordinates into the struct
%   OUTPUT.
%
%   Version 2.0
%   Copyright (c) 06/28/2022 Benjamin Seibold, Martin Frank, and
%                            Rujeko Chinomona
%   http://www.math.temple.edu/~seibold
%   https://www.scc.kit.edu/personen/martin.frank.php
%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0),
%
%   StaRMAP project website:
%   https://github.com/starmap-project

%========================================================================
% License
%========================================================================
%   Permission is hereby granted, free of charge, to any person obtaining
%   a copy of this software and associated documentation files (the
%   "Software"), to deal in the Software without restriction, including
%   without limitation the rights to use, copy, modify, merge, publish,
%   distribute, sublicense, and/or sell copies of the Software, and to
%   permit persons to whom the Software is furnished to do so, subject to
%   the following conditions:
%   * The above copyright notice and this permission notice shall be
%     included in all copies or substantial portions of the Software.
%   * Appropriate references shall be given to the name of the Software
%     (StaRMAP), the website of the Software
%     (https://github.com/starmap-project/starmap)
%     and the appropriate research article(s) corresponding to the
%     Software, in particular:
%     B. Seibold, M. Frank, ACM Transactions on Mathematical Software,
%     41(1) [4]. DOI: 10.1145/2590808 (http://arxiv.org/abs/1211.2205)
%
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
%   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
%   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
%   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.


%========================================================================
% Set parameters and necessary flags
%========================================================================
n_sys = length(par.Mz);                    % Number of system components.
if isfield(par,'mom_order')
    f_mom_order = par.mom_order;          % Create copy of par.mom_order.
end
flag_filter = isfield(par,'filterfunction')&&isfield(par,'f_position')...
    &&isfield(par,'mom_order');
if flag_filter, f_mom_order = par.mom_order;
else, par.filterfunction = @ones; par.f_position = '';
end
flag_ani_scatter = isfield(par,'sigma_sm')&&isfield(par,'mom_order');
if ~flag_ani_scatter           % If anisotropic part of scattering kernel
    par.mom_order = [0 ones(1,n_sys-1)]; par.sigma_sm = @zero;      % not
end                    % available, set up data for isotropic scattering.
if strcmp(par.mom_output,'all'), par.mom_output = 1:n_sys; end

%========================================================================
% Initialization
%========================================================================
% Create index matrices and vectors for grid type
Ix = cellfun(@find,num2cell(par.Mx',1),'Un',0);         % Grid dependency
Iy = cellfun(@find,num2cell(par.My',1),'Un',0);                % indices.
Iz = cellfun(@find,num2cell(par.Mz',1),'Un',0);
c111 = 1; c122 = []; c121 = []; c112 = [];        % First component is on
c211 = []; c222 = []; c221 = []; c212 = [];                  %(1,1)-grid.
while length([c111,c122,c121,c112,c211,c222,c221,c212])<n_sys % Determine
    c211 = unique(vertcat(c211',Ix{c111},Iy{c221},Iz{c212})');  % indices
    c121 = unique(vertcat(c121',Ix{c221},Iy{c111},Iz{c122})');  % of grid
    c112 = unique(vertcat(c112',Ix{c212},Iy{c122},Iz{c111})');    % types
    c221 = unique(vertcat(c221',Ix{c121},Iy{c211},Iz{c222})');     % from
    c212 = unique(vertcat(c212',Ix{c112},Iy{c222},Iz{c211})');   % system
    c122 = unique(vertcat(c122',Ix{c222},Iy{c112},Iz{c121})');% matrices.
    c111 = unique(vertcat(c111',Ix{c211},Iy{c121},Iz{c112})');
    c222 = unique(vertcat(c222',Ix{c122},Iy{c212},Iz{c221})');
end
if isempty(c221), c221 = []; end
if isempty(c222), c222 = []; end

gtx = ones(1,n_sys); gty = gtx; gtz = gtx;                  % Grid types
gtx([c211,c212,c221,c222]) = 2;                              % in x,y,z.
gty([c121,c122,c221,c222]) = 2;
gtz([c112,c122,c212,c222]) = 2;

% Grid size and time step
h = (par.ax([2 4 6])-par.ax([1 3 5]))./par.n;                % Grid size.
dt_hyp = min(h)/abs(eigs(par.Mz,1,'lm'));            % Maximal time step.
if par.n(1)~=1 % Adapt dt to 3D test case.
    dt_hyp = dt_hyp/3;
elseif par.n(2)~=1 % Adapt dt to 2D test case.
    dt_hyp = dt_hyp/2;
end
M = -1i*(par.Mx/h(1)+par.My/h(2)+par.Mz/h(3)); % Scaled matrices.
time_dependent = [nargin(par.sigma_a)>3,nargin(par.sigma_s0)>3,...
            nargin(par.sigma_sm)>4];

% Create staggered grids and coefficients
x{1} = par.ax(1)+h(1)/2:h(1):par.ax(2)-h(1)/2;              % Cell center
y{1} = par.ax(3)+h(2)/2:h(2):par.ax(4)-h(2)/2;             % coordinates.
z{1} = par.ax(5)+h(3)/2:h(3):par.ax(6)-h(3)/2;
x{2} = par.ax(1)+h(1)*(1-par.bc(1)):h(1):par.ax(2);           % Staggered
y{2} = par.ax(3)+h(2)*(1-par.bc(2)):h(2):par.ax(4);        % coordinates.
z{2} = par.ax(5)+h(3)*(1-par.bc(3)):h(3):par.ax(6);

X = cell(2,2,2); Y = X; Z = X;                         % Staggered grids.
[X(:,:,1),Y(:,:,1),Z(:,:,1)] = cellfun(@ndgrid,[x;x]',[y;y],...
    {z{1},z{1};z{1},z{1}},'Un',0);
[X(:,:,2),Y(:,:,2),Z(:,:,2)] = cellfun(@ndgrid,[x;x]',[y;y],...
    {z{2},z{2};z{2},z{2}},'Un',0);

n1 = par.n; n2 = size(X{2,2,2});              % Sizes of staggered grids.
extendx = {[1:par.bc(1),1:n1(1),par.bc(1)*(n1(1)-1)+1],...    % Extension
    [n2(1)*(1:1-par.bc(1)),1:n2(1)]};                       % indices for
extendy = {[1:par.bc(2),1:n1(2),par.bc(2)*(n1(2)-1)+1],...     % boundary
    [n2(2)*(1:1-par.bc(2)),1:n2(2)]};                        % conditions.
extendz = {[1:par.bc(3),1:n1(3),par.bc(3)*(n1(3)-1)+1],...
    [n2(3)*(1:1-par.bc(3)),1:n2(3)]};
sg = [gtx;gty;gtz;par.mom_order]';             % Construct staggered grid
sg = unique(sg(2:end,:),'rows')';            % indices with moment order.

% Initialize field variables
sS = cell(2,2,2,max(par.mom_order));                 % Scattering fields.
sT = sS; ET = sS;
U = cell(1,n_sys); dxU = U; dyU = U; dzU = U; Q = U;          % Unknowns.
for j = 1:n_sys                         % If initial condition only three
    U{j} = X{gtx(j),gty(j),gtz(j)}*0+...                 % arguments, set
        capargs(par.ic,X{gtx(j),gty(j),gtz(j)},...   % all higher moments
        Y{gtx(j),gty(j),gtz(j)},Z{gtx(j),gty(j),gtz(j)},j)*... % to zero.
        (nargin(par.ic)>3||j==1);
    Q{j} = X{gtx(j),gty(j),gtz(j)}*0; %Initialize source terms with zeros.
end

% Initialize filtering variables
if flag_filter
    if strcmp(par.f_position,'full')
        filterTerm = cell(2,2,2,max(f_mom_order));    % Full filter term.
    elseif strcmp(par.f_position,'substep')
        filterTerm = cell(2,2,2,max(f_mom_order),2); %Substep filter term.
    end
    f_sg = [gtx;gty;gtz;f_mom_order]';      % Staggered grid indices with
    f_sg = unique(f_sg(2:end,:),'rows')';                 % moment order.
end

%========================================================================
% Time Loop
%========================================================================
par.t_plot = [par.t_plot inf]; plot_count = 1;
cputime = zeros(1,3);
t = 0; dt = 0; step_count = 0;
while t<par.tfinal      % Loop until final time is reached (or exceeded).
    % Evaluate material parameters and determine maximal time step
    tic
    tleft = (par.tfinal-t)/par.cn;  % Time remaining, scaled with 1/CFL#.
    for k = double(t>0):1      % Loop twice for t=0, otherwise only once.
        if ~any(time_dependent)&&t>0, break, end
        flag_check_stab = k>0&&par.timestep>0;    % If not first time and
                                                        % parabolic step.
        if flag_check_stab
            % Determine maximum admissible time step
            dt = max(dt,dt_hyp);  % Initial guess at least hyperbolic dt.
            decay = ones(1,n_sys)*min(min(sA{1,1,1})); %Smallest sigma_a.
            for j = 2:n_sys        % Smallest total cross-section values.
                decay(j) = min(min(sT{gtx(j),gty(j),...
                    gtz(j),par.mom_order(j)}));
            end
            la = growth_factor(par,M,[c111,c221,c212,c122],...
                [c211,c121,c112,c222],decay,dt);
            dt_bracket = [1 1]*dt; la_bracket = [1 1]*la;
            la0 = 1+1e-8;  % Allow eigenv. bound of 1+a bit (->no decay).
            while ~diff(la_bracket<la0)...%While eigs on same side of la0
                &&(dt_bracket(1)<=tleft...    % and (t_final not exceeded
                ||la_bracket(1)>la0)                % or dt(1) unstable).
                ind = (la_bracket(1)<la0)+1;% dt(1) unstable:1; stable:2.
                dt_bracket = dt_bracket(ind)*[1 2]*2^(ind-2); % Scale dt.
                la_bracket = la_bracket([2 1]); % Overwrite new time step
                la_bracket(ind) = growth_factor(par,M,...  % with new ei-
                    [c111,c221,c212,c122],...              %genvalue.
                    [c211,c121,c112,c222],decay,dt_bracket(ind));
                if min(la_bracket)<1, la0 = 1; end% Lower eig bound to 1.
            end
            while dt_bracket(1)/dt_bracket(2)<0.98...   % Bisection iter-
                &&dt_bracket(1)<=tleft   % ation, if not reaching tfinal.
                dt_mid = sqrt(prod(dt_bracket)); % Geom. mid point of dt.
                la_mid = growth_factor(par,M,...          % Eigenvalue at
                    [c111,c221,c212,c122],...             % midpoint dt.
                    [c211,c121,c112,c222],decay,dt_mid);
                ind = 1+(la_mid>=la0); % Index in bracket to be replaced.
                dt_bracket(ind) = dt_mid; la_bracket(ind) = la_mid;
            end
        else, dt_bracket = dt_hyp;             % Use hyperbolic time step.
        end
        if k>0                                  % If not very first time.
            dt = min(dt_bracket(1),tleft)*par.cn;      %Actual time step.
            sA0 = sA; sT0 = sT;   % Copies of parameters at current time.
            if ~any(time_dependent)% If autonom.
                sA1 = sA{1,1,1}; sT1 = sT; break % use current parameters
            end                                          % and exit loop.
        end
        % Evaluate material parameters at t+dt
        evaluate = time_dependent|(dt==0&[1 1 flag_ani_scatter]);
        if evaluate(1)                            % Absorption: evaluate,
            sA = cellfun(@(x,y,z)capargs(...      % if time-dependent or
                par.sigma_a,x,y,z,t+dt),X,Y,Z,'Un',0); % first time step.
        end
        if evaluate(2)        % Scattering, isotropic component: evaluate
            for j = 1:size(sg,2)  % if time-dependent or first time step.
                sS{sg(1,j),sg(2,j),sg(3,j),sg(4,j)} = ...
                    capargs(par.sigma_s0,X{sg(1,j),sg(2,j),sg(3,j)},...
                    Y{sg(1,j),sg(2,j),sg(3,j)},...
                    Z{sg(1,j),sg(2,j),sg(3,j)},t+dt);
            end
        end
        if evaluate(3)     % Scattering, anisotropic components: evaluate
            for j = 1:size(sg,2)  % if time-dependent or first time step.
                sS{sg(1,j),sg(2,j),sg(3,j),sg(4,j)} = ...
                    sS{sg(1,j),sg(2,j),sg(3,j),sg(4,j)}-...
                    capargs(par.sigma_sm,X{sg(1,j),sg(2,j),sg(3,j)}, ...
                    Y{sg(1,j),sg(2,j),sg(3,j)}, ...
                    Z{sg(1,j),sg(2,j),sg(3,j)},sg(4,j),t+dt);
            end
        end
        if any(evaluate) % Total sigma_t^m = sigma_a+sigma_s^0-sigma_s^m.
            for j = 1:size(sg,2)             % (Re)compute, if any of the
                sT{sg(1,j),sg(2,j),sg(3,j),sg(4,j)} = ...         % above
                    sA{sg(1,j),sg(2,j),sg(3,j)}+...       %(re)evaluated.
                    sS{sg(1,j),sg(2,j),sg(3,j),sg(4,j)};
            end
        end
        % Determine material parameters used to conduct time step
        if k>0                                  % If not very first time.
            sA1 = (sA0{1,1}+sA{1,1})/2; sT1 = sT;               % Average
            for j = 1:size(sg,2)                               % material
                sT1{sg(1,j),sg(2,j),sg(3,j),sg(4,j)} = ...   % parameters
                    (sT0{sg(1,j),sg(2,j),sg(3,j),sg(4,j)}+....     % at t
                    sT{sg(1,j),sg(2,j),sg(3,j),sg(4,j)})/2;   % and t+dt.
            end
        end
        if flag_check_stab % If stability investigations were done above.
            % Check whether avg. material parameters are stable with dt
            decay = ones(1,n_sys)*min(min(sA1));      % Smallest sigma_a.
            for j = 2:n_sys        % Smallest total cross-section values.
                decay(j) = min(min(sT1{gtx(j),gty(j),gtz(j),...
                    par.mom_order(j)}));
            end
            if ~growth_factor(par,M,...% If step with averaged parameters
                [c111,c221,c212,c122],[c211,c121,c112,c222],...      % is
                decay,dt)<la0                       % unstable, fall back
                sA1 = sA0{1,1,1}; sT1 = sT0;   % to parameters at time t.
            end
        end
    end
    % Compute decay terms
    if t==0||any(time_dependent)                 % Recompute decay terms.
        EA = expm1div(-sA1*dt/2);         % Decay term for zeroth moment.
        for j = 1:size(sg,2)
            ET{sg(1,j),sg(2,j),sg(3,j),sg(4,j)} = ... % Decay of sigma_t.
                expm1div(-sT1{sg(1,j),sg(2,j),sg(3,j),sg(4,j)}*dt/2);
        end
    end
    % Time increment
    t = t+dt; step_count = step_count+1;
    % Evaluate source
    for j = par.source_ind     % Source: evaluate only active components.
        Q{j} = capargs(par.source,X{gtx(j),gty(j),gtz(j)},...  % Evaluate
            Y{gtx(j),gty(j),gtz(j)},Z{gtx(j),gty(j),gtz(j)},...  % source
            t-dt/2,j);                            % at half-time of step.
    end
    % Evaluate filter function
    if (t==dt||nargin(par.filterfunction)>6)&&flag_filter  % Evaluate, if
        for j = 1:size(f_sg,2)       % first time step or time-dependent.
            if strcmp(par.f_position,'full')
                filterTerm{f_sg(1,j),f_sg(2,j),f_sg(3,j),f_sg(4,j)} = ...
                    capargs(par.filterfunction,par,dt,f_sg(4,j),...
                    X{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Y{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Z{f_sg(1,j),f_sg(2,j),f_sg(3,j)},t);
            elseif strcmp(par.f_position,'substep')
                filterTerm{f_sg(1,j),f_sg(2,j),f_sg(3,j),f_sg(4,j),1} = ...
                    capargs(par.filterfunction,par,dt,f_sg(4,j),...
                    X{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Y{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Z{f_sg(1,j),f_sg(2,j),f_sg(3,j)},t-dt/2);
                filterTerm{f_sg(1,j),f_sg(2,j),f_sg(3,j),f_sg(4,j),2} = ...
                    capargs(par.filterfunction,par,dt,f_sg(4,j),...
                    X{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Y{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Z{f_sg(1,j),f_sg(2,j),f_sg(3,j)},t);
            end
        end
    end
    cputime(2) = cputime(2)+toc;
    % Update solution
    U0 = U(par.mom_output);                     % Store current solution.
    k = 1; tic
    for step = [1 2 1]                   % Strang splitting of sub-steps.
        switch step
        case 1         % Half step with (111), (221), (212), (122) grids.
        for j = c111                         % Compute update from (111).
            if ~isempty(Ix{j})
                dxU{j} = diff(U{j}(extendx{1},:,:),[],1)/h(1);
            end
            if ~isempty(Iy{j})
                dyU{j} = diff(U{j}(:,extendy{1},:),[],2)/h(2);
            end
            dzU{j} = diff(U{j}(:,:,extendz{1}),[],3)/h(3);
        end
        for j = c221                         % Compute update from (221).
            if ~isempty(Ix{j})
                dxU{j} = diff(U{j}(extendx{2},:,:),[],1)/h(1);
            end
            if ~isempty(Iy{j})
                dyU{j} = diff(U{j}(:,extendy{2},:),[],2)/h(2);
            end
            dzU{j} = diff(U{j}(:,:,extendz{1}),[],3)/h(3);
        end
        for j = c212                         % Compute update from (212).
            if ~isempty(Ix{j})
                dxU{j} = diff(U{j}(extendx{2},:,:),[],1)/h(1);
            end
            if ~isempty(Iy{j})
                dyU{j} = diff(U{j}(:,extendy{1},:),[],2)/h(2);
            end
            dzU{j} = diff(U{j}(:,:,extendz{2}),[],3)/h(3);
        end
        for j = c122                         % Compute update from (122).
            if ~isempty(Ix{j})
                dxU{j} = diff(U{j}(extendx{1},:,:),[],1)/h(1);
            end
            if ~isempty(Iy{j})
                dyU{j} = diff(U{j}(:,extendy{2},:),[],2)/h(2);
            end
            dzU{j} = diff(U{j}(:,:,extendz{2}),[],3)/h(3);
        end
        for j = [c211 c121 c112 c222]       % Update components on grids.
            W = -sumcell([dxU(Ix{j}),dyU(Iy{j}),dzU(Iz{j})],...
                [par.Mx(j,Ix{j}),par.My(j,Iy{j}),par.Mz(j,Iz{j})]);
            U{j} = U{j}+dt/2*(W+Q{j}-...
                sT1{gtx(j),gty(j),gtz(j),par.mom_order(j)}.*U{j}).*...
                ET{gtx(j),gty(j),gtz(j),par.mom_order(j)};
            if strcmp(par.f_position,'substep')       % Multiplication by
                U{j} = U{j}.*filterTerm{gtx(j),gty(j),...       % substep
                    gtz(j),f_mom_order(j),k};              % filter term.
            end
        end
        case 2         % Half step with (211), (121), (112), (222) grids.
        for j = c222                         % Compute update from (222).
            if ~isempty(Ix{j})
                dxU{j} = diff(U{j}(extendx{2},:,:),[],1)/h(1);
            end
            if ~isempty(Iy{j})
                dyU{j} = diff(U{j}(:,extendy{2},:),[],2)/h(2);
            end
            dzU{j} = diff(U{j}(:,:,extendz{2}),[],3)/h(3);
        end
        for j = c112                         % Compute update from (112).
            if ~isempty(Ix{j})
                dxU{j} = diff(U{j}(extendx{1},:,:),[],1)/h(1);
            end
            if ~isempty(Iy{j})
                dyU{j} = diff(U{j}(:,extendy{1},:),[],2)/h(2);
            end
            dzU{j} = diff(U{j}(:,:,extendz{2}),[],3)/h(3);
        end
        for j = c121                         % Compute update from (121).
            if ~isempty(Ix{j})
                dxU{j} = diff(U{j}(extendx{1},:,:),[],1)/h(1);
            end
            if ~isempty(Iy{j})
                dyU{j} = diff(U{j}(:,extendy{2},:),[],2)/h(2);
            end
            dzU{j} = diff(U{j}(:,:,extendz{1}),[],3)/h(3);
        end
        for j = c211                         % Compute update from (211).
            if ~isempty(Ix{j})
                dxU{j} = diff(U{j}(extendx{2},:,:),[],1)/h(1);
            end
            if ~isempty(Iy{j})
                dyU{j} = diff(U{j}(:,extendy{1},:),[],2)/h(2);
            end
            dzU{j} = diff(U{j}(:,:,extendz{1}),[],3)/h(3);
        end
        for j = [c111 c221 c212 c122]       % Update components on grids.
            W = -sumcell([dxU(Ix{j}),dyU(Iy{j}),dzU(Iz{j})],...
                [par.Mx(j,Ix{j}),par.My(j,Iy{j}),par.Mz(j,Iz{j})]);
            for k = 1:2                         % Perform two half-steps.
                if j==1        % Zeroth moment decays by absorption only.
                    U{j} = U{j}+dt/2*(W+Q{j}-sA1.*U{j}).*EA;
                else                  % All other moments decay normally.
                    U{j} = U{j}+dt/2*(W+Q{j}-...
                        sT1{gtx(j),gty(j),gtz(j),par.mom_order(j)}.*...
                        U{j}).*ET{gtx(j),gty(j),gtz(j),par.mom_order(j)};
                    if strcmp(par.f_position,'substep')  % Multiplication
                        U{j} = U{j}.*...                     % by substep
                            filterTerm{gtx(j),gty(j),gtz(j),...  % filter
                            f_mom_order(j),k};                    % term.
                    end
                end
            end
        end
        end
    end
    if strcmp(par.f_position,'full') % Multiplication by full filter term
        for j= 2:n_sys                         % at the end of time step.
            U{j} = U{j}.*filterTerm{gtx(j),gty(j),gtz(j),f_mom_order(j)};
        end
    end
    cputime(1) = cputime(1)+toc;
    % Plotting
    tic
    while t>=par.t_plot(plot_count)-1e-14  % If current time has exceeded
        lambda = (par.t_plot(plot_count)-t+dt)/dt; %plotting time, define
        Uplot = cellfun(@(x,y)(1-lambda)*x+lambda*y,...    % solution via
            U0,U(par.mom_output),'Un',0); % linear interpolation in time.
        xplot = x(gtx(par.mom_output));                 % Assign grids at
        yplot = y(gty(par.mom_output));              % outputted moments.
        zplot = z(gtz(par.mom_output));
        switch par.dim  % Configure output according to problem dimension
            case 1
                for j = 1:length(Uplot)                  % Reshape for 1D
                    Uplot{j} = reshape(Uplot{j},length(zplot{j}),[]);
                end
                % If only a single moment plotted, remove cell structure.
                if length(par.mom_output)==1
                    zplot = zplot{:}; Uplot = Uplot{:};
                end
                if nargout(par.output)             % Call output routine,
                    par = par.output(par,zplot,Uplot,plot_count); % allo-
                else                          % wing for it to modify the
                    par.output(par,zplot,Uplot,plot_count)  % struct par.
                end
            case 2
                for j = 1:length(Uplot)                 % Reshape for 2D
                    Uplot{j} = reshape(Uplot{j},length(yplot{j}),[]);
                end
                % If only a single moment plotted, remove cell structure.
                if length(par.mom_output)==1
                    yplot = yplot{:}; zplot = zplot{:}; Uplot = Uplot{:};
                end
                if nargout(par.output)             % Call output routine,
                    par = par.output(par,yplot,zplot,Uplot,plot_count);
                else               % allowing for it to modify the struct
                    par.output(par,yplot,zplot,Uplot,plot_count)   % par.
                end
            case 3
                % If only a single moment plotted, remove cell structure.
                if length(par.mom_output)==1
                    xplot = xplot{:}; yplot = yplot{:};...
                        zplot = zplot{:}; Uplot = Uplot{:};
                end
                if nargout(par.output)             % Call output routine,
                    par = par.output(par,xplot,yplot,...       % allowing
                        zplot,Uplot,plot_count);       % for it to modify
                else                                    % the struct par.
                    par.output(par,xplot,yplot,zplot,Uplot,plot_count)
                end
        end
        plot_count = plot_count+1;
    end
    cputime(3) = cputime(3)+toc;
end
% Output problem overview and statistics
fprintf('%s with %s%s%d\n %d moments, %dx%dx%d grid, %0.0f time steps\n',...
    par.name,par.closure_mod,...
    par.closure,par.n_mom,n_sys,par.n,step_count) % Display test
cputime = reshape([cputime;cputime/sum(cputime)*1e2],1,[]);   % case info
fprintf(['CPU-times\n advection:%12.2fs%5.0f%%\n ',...   % and CPU times.
    'absorb/scatter:%7.2fs%5.0f%%\n plotting:%13.2fs%5.0f%%\n'],cputime)

% Output final state of solution
if nargout
    switch par.dim
        case 1
            for j = 1:length(U)
                U{j} = reshape(U{j},length(z{gtz(j)}),[]);
            end
            output = struct('x',z(gtz),'U',U);
        case 2
            for j = 1:length(U)
                U{j} = reshape(U{j},length(y{gty(j)}),[]);
            end
            output = struct('x',y(gty),'y',z(gtz),'U',U);
        otherwise
           output = struct('x',x(gtx),'y',y(gty),'z',z(gtz),'U',U);
    end
end


%========================================================================
% Technical Functions
%========================================================================
function la = growth_factor(par,M,c_even,c_odd,decay,dt)
% Calculate largest eigenvalue of von-Neumann growth factor matrix.
G = diag(dt*expm1div(-decay*dt/2))*(M+diag(-decay/2));
I = eye(length(decay));
G_e = G; G_e(c_odd,:) = 0; G_e = G_e+I;    %Growth factor of even update.
G_o = G; G_o(c_even,:) = 0; G_o = G_o+I;   % Growth factor of odd update.
G = G_o*G_e*G_e*G_o;            % Growth factor matrix of full time step.
la = max(abs(eig(G)));      % Largest eigenvalue of growth factor matrix.
if par.timestep==2, la = la*exp(decay(1)*dt); end   % Relative stability.

function S = sumcell(A,w)
% Add vector of cells A, weighted with vector w.
w = full(w);
S = A{1}*w(1); for j = 2:length(w), S = S+A{j}*w(j); end

function z = capargs(fct,varargin)
% Call function fct with as many arguments as it requires (at least 1),
% and ignore further arguments.
narg = max(nargin(fct),1);
z = fct(varargin{1:narg});

function f = zero(varargin)
% Zero function.
f = zeros(size(varargin{1}));

function y = expm1div(x)
% Function (exp(x)-1)/x that is accurate for x close to zero.
y = 1+x*.5+x.^2/6;
ind = abs(x)>2e-4;
y(ind) = (exp(x(ind))-1)./x(ind);
