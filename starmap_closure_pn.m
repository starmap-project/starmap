function [moment_order,Mx,My,Mz] = starmap_closure_pn(n_mom,dim)
%STARMAP_CLOSURE_PN
%   Creates P_N moment system matrices, to be used by
%   STARMAP_SOLVER, a second order staggered grid finite
%   difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 1D-3D geometries.
%   In 1D : Moment matrix stored in Mz. Mx and My are zero matrices.
%   In 2D : Moment matrices My (for x-coordinate) and Mz (for y-coordinate)
%           Mx is a zero matrix
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

switch dim
    case 1
    n_sys = n_mom+1;
    % Moment matrix
    subdiag = zeros(n_sys,1);
    superdiag = zeros(n_sys,1);
    for k=0:n_sys-2
        subdiag(k+1) = (k+1)/(2*k+3);
        superdiag(k+2) = (k+1)/(2*k+1);
    end
    Mz = spdiags([subdiag superdiag],[-1 1],n_sys,n_sys);
    Mx = sparse(zeros(n_sys));
    My = sparse(zeros(n_sys));
    % Order of moments
    moment_order = n_mom+1;
    case 2
    n_sys = (n_mom+1)*(n_mom+2)/2; % number of system components
    Mx = sparse(zeros(n_sys)); My = Mx; s = size(Mx);
    for m = 1:n_mom
        i = 1:m; p = m*(m-1)/2+i;
        v = D(m,-m+2*(ceil(i/2)-1));
        Mx(sub2ind(s,p,p+m)) = v;
        My(sub2ind(s,p,p+m-(-1).^i)) = -(-1).^i.*v;
        i = 1:m-1; p = m*(m-1)/2+i;
        v = F(m,-m+2+2*(ceil(i/2)-1));
        Mx(sub2ind(s,p,p+m+2)) = -v;
        My(sub2ind(s,p-(-1).^i,p+m+2)) = (-1).^i.*v;
    end
    m = 1:2:n_mom; i = m.*(m+1)/2;
    Mx(i,:) = sqrt(2)*Mx(i,:); My(i,:) = sqrt(2)*My(i,:);
    m = 2:2:n_mom; i = (m+1).*(m+2)/2;
    Mx(:,i) = sqrt(2)*Mx(:,i); My(:,i) = sqrt(2)*My(:,i);
    Mx = full((Mx+Mx')/2); My = full((My+My')/2);
    % Adopt starmap solver ordering
    Mz = My; My=Mx; Mx = sparse(zeros(n_sys));
    % Order of moments
    moment_order = ceil(sqrt(2*(1:n_sys)+1/4)-3/2);
    case 3
    n_sys = (n_mom+1)^2; % number of system components in 3D
    Mx = sparse(zeros(n_sys)); My = Mx; Mz = Mx; s = size(Mx);
    for m = 1:n_mom
        i = 1:2*m-1; p = (m-1)^2+i;
        v = F(m,m+1-ceil(i/2));
        Mx(sub2ind(s,p,p+2*m-1)) = v;
        My(sub2ind(s,p,p+2*m-1-(-1).^i)) = -(-1).^i.*v;
        v = B(m,m-ceil(i/2));
        Mz(sub2ind(s,p,p+2*m+1)) = v;
        i = 1:2*m-3; p = (m-1)^2+i;
        v = D(m,m-1-ceil(i/2));
        Mx(sub2ind(s,p,p+2*m+3)) = -v;
        if m>1, p(end) = p(end)+1; i(end) = i(end)+1; end
        My(sub2ind(s,p,p+2*m+3-(-1).^i)) = -(-1).^i.*v;
    end
    m = 1:n_mom; i = m.^2;
    Mx(i,:) = sqrt(2)*Mx(i,:); My(i,:) = sqrt(2)*My(i,:);
    m = 2:n_mom; i = (m+1).^2;
    Mx(:,i) = sqrt(2)*Mx(:,i); My(:,i) = sqrt(2)*My(:,i);
    Mx = full((Mx+Mx')/2); My = full((My+My')/2); % symmetry and factor 1/2
    Mz = full(Mz+Mz');
    % Order of moments in 3D
    moment_order = ceil(sqrt((1:n_sys))-1);
end

function y = D(l,m)
y = sqrt((l-m).*(l-m-1)/(2*l+1)/(2*l-1));

function y = F(l,m)
y = sqrt((l+m).*(l+m-1)/(2*l+1)/(2*l-1));

function y = B(l,m)
y = sqrt((l-m).*(l+m)/(2*l+1)/(2*l-1));
