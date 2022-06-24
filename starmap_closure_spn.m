function [moment_order,Mx,My,Mz] = starmap_closure_spn(n_mom,dim)
  %STARMAP_CLOSURE_SPN
  %   Creates SP_N moment system matrices, to be used by
  %   STARMAP_SOLVER, a second order staggered grid finite
  %   difference solver for linear hyperbolic moment
  %   approximations to radiative transfer in 1D-3D geometry.
  %   In 1D : Moment matrix stored in Mz. Mx and My are zero matrices.
  %   In 2D : Moment matrices My (for x-coord) and Mz (for y-coord)
  %           Mx is a zero matrix
  %
  %   The time-dependent SPN equations encoded here are
  %   the ones given in
  %   [Olbrant, Larsen, Frank, Seibold, JCP 238 (2013) 315-336].
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
    nstar = ceil((n_mom+1)/2);
    i = 0:nstar-1;
    k = 2*i.*(2*i-1)./((4*i+1).*(4*i-1));
    l = 4*i.^2./((4*i+1).*(4*i-1))+...
    (2*i+1).^2./((4*i+1).*(4*i+3)).*(i<nstar-1|floor(n_mom/2)*2<n_mom);
    m = 2*(2*i+1).*(i+1)./((4*i+1).*(4*i+3));
    % Create matrix for x-direction
    Mx = zeros(nstar*3);
    Mx(nstar+(1:nstar)*2-1,1:nstar) = ...
    diag(k(2:end),-1)+diag(l)+diag(m(1:end-1),1);
    Mx(1:nstar,nstar+(1:nstar)*2-1) = eye(nstar);
    % Create matrix for y-direction
    i = 1:nstar*2; i = [1:nstar,nstar+i+(ceil(i/2)*2~=i)*2-1];
    My = Mx(i,i);
    % Adopt starmap_solver ordering
    Mz = My; My=Mx; Mx = sparse(zeros(size(Mz)));
    % Order of moments
    moment_order = [0:2:n_mom reshape([1;1]*((0:2:n_mom)+1),1,[])];

    case 3
    nstar = ceil((n_mom+1)/2);
    i = 0:nstar-1;
    k = 2*i.*(2*i-1)./((4*i+1).*(4*i-1));
    l = 4*i.^2./((4*i+1).*(4*i-1))+...
    (2*i+1).^2./((4*i+1).*(4*i+3)).*(i<nstar-1|floor(n_mom/2)*2<n_mom);
    m = 2*(2*i+1).*(i+1)./((4*i+1).*(4*i+3));
    % Create matrix for x-direction
    Mx = zeros(nstar*4);
    Mx(nstar+(1:nstar)*3-2,1:nstar) = ...
    diag(k(2:end),-1)+diag(l)+diag(m(1:end-1),1);
    Mx(1:nstar,nstar+(1:nstar)*3-2) = eye(nstar);
    % Create matrix for y-direction
    i = 1:nstar*3; i = [1:nstar,nstar+i+mod(i+1,3)-1];
    My = Mx(i,i);
    % Create matrix for z-direction
    i = 1:nstar*3; i = [1:nstar,nstar+i+2*(1-mod(i+2,3))];
    Mz = Mx(i,i);
    % Order of moments
    moment_order = [0:2:n_mom reshape([1;1;1]*((0:2:n_mom)+1),1,[])];

  end
