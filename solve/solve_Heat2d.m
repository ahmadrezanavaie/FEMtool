function [ u_h, varargout ] = solve_Heat2d( problem, method, output )
%SOLVE_HEAT2D solves a 2-d heat problem by means of the Finite Element
%Method (FEM). The differential problem reads:
%
% find u : \Omega x ]0,T] --> R s.t.
%
%     	du/dt - grad( \mu div(u) ) = f(x,t),     in \Omega x ]0,T]
%                           u(.,0) = u_0,        in \Omega
%                           u(x,t) = g(x,t),     on \Gamma_D x ]0,T]
%                   \mu du(x,t)/dn = h(x,t),     on \Gamma_N x ]0,T]
%
%
%The computational domain is defined in the file geo_2d.m through a
%parameterization from a rectangular parametric domain (usually the unit
%square) into the physical domain in R^2:
%
%        _________________
%       |                 |
%       |       (4)       |
%       |                 |
%       |                 |       F
%       | (1)         (2) |     -----> \Omega
%       |                 |
%  p2   |                 |
%  ^    |       (3)       |
%  |    |_________________|
%  |
% -|-----> p1
%
%The boundary in the computational domain is defined by the mapping
%
%       \Gamma = F( (1) U (2) U (3) U (4) ).
%
%On each part of the boundary, the boundary conditions are imposed by three
%parameters:
%
%   - index of the boundary part
%   - relative limits of the boundary condition type
%   - function handle.
%
%Let us consider an example.
% Example BC:
%   \Gamma_1 = F( (1) )
%
%         Dirichlet BC      Neumann BC
%            g_1(x)           h_1(x)
%       |--------------|-------------------|
%       0             0.4                 1.0     (relative length)
%
%   The two boundary conditions are imposed by
%       - dir_sides = [1];      % boundary index
%       - dir_lim = [ 0; 0.4];  % relative limits
%       - dir_fun{1} = g_1;     % function handle
%       - neum_sides = [1];     % boundary index
%       - neum_lim = [ 0.4; 1]; % relative limits
%       - neum_fun{1} = h_1;    % function handle.
%   
%
% -------------------------------------------------------------------------
%
%   u_h = solve_Heat2d( problem, method, output );
%   [ u_h, data ] = solve_Heat2d( problem, method, output );
%   [ u_h, data, mesh, space ] = solve_Heat2d( problem, method, output );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% problem           (struct)
%   .geo            (see output geo_2d.m)
%   .BC             (struct) See example above
%       .dir_sides  (1 x nDir double) Dirichlet boundary indices
%       .dir_lim    (2 x nDir double) Dirichlet limits
%       .dir_fun    (1 x nDir cell array) Dirichlet function handles
%       .neum_sides (1 x nNeum double) Neumann boundary indices
%       (if .neum_sides non-empty)
%       .neum_lim   (2 x nNeum double) Neumann limits
%       .neum_fun   (1 x nNeum cell array) Neumann function handles
%   .mu             (function handle) or (1 x #mesh elements double)
%                   Diffusion coefficient
%   .sigma          (function handle) or (1 x #mesh elements double)
%                   Storativity coefficient
%   .f              (function handle) or (1 x #mesh elements double)
%                   Reaction term
%   .u_0            (function handle) For the initial solution
%   .Tmax           (1 x 1 double) Simulation ending time
% method
%   .nElDir         (2 x 1 double) Number of mesh elements in each
%                   parametric direction
%   .nQRDir         (2 x 1 double) Number of quadrature nodes per element
%                   in each parametric direction
%   .rBFDir         (2 x 1 double) Order of the piecewise polynomial basis
%                   functions in each parametric direction
%   .lMax           (1 x 1 integer) Maximal number of loops
%   .dt             (1 x 1 double) Time step (initial or fixed)
%   .adapt_dt       (struct)
%       .flag       (1 x 1 logical) Use fixed (=0) or adaptive (=1) time
%                   step selection
%       (if flag == 1)
%       .facmin     (1 x 1 double) Max decrease factor of time step
%       .facmax     (1 x 1 double) Max increase factor of time step
%       .fac        (1 x 1 double) Factor for time step selection
%       .tol        (1 x 1 double) Tolerance for error control
%   .timeScheme     (string)
%       'EE': For Euler Explicit (no time step adaption possible)
%       'EI': For Euler Implicite
%       'CN': For Crank Nicolson
%       'RKC2': For 2nd order damped Runge Kutta Chebyshev
%       'RIIA': For Radau IIA Runge Kutta
%       'GL4': For 4th order Gauss-Legendre Runge Kutta
%   .eta            (1 x 1 double) Damping factor (for 'RKC2' only)
%   .discIn         (string) Initial function discretization
%       'IP': For interpolation
%       'L2': For L^2 projection on the discrete space
%       'VW': If u_0 is a vector of weights of size (nDof x 1 double)
% output
%   .screen_flag    (1 x 1 logical) Flag for writing percentage of time
%                   integration on the screen
%   .vtk_flag       (1 x 1 logical) Flag for writing .vtk files of solution
%   (if vtk_flag == 1)
%   .vtk_name       (string) Name of output file(s) (f.e: /path/filename_ )
%   .vtk_freq       (1 x 1 integer) Frequence of output (usually = 1)
%   (optional)
%   .lmA_flag       (1 x 1 logical) Flag for computing the eigenvalue of
%                   largest magnitude of the stiffness matrix A
%   .smA_flag       (1 x 1 logical) Flag for computing the eigenvalue of
%                   smallest magnitude of the stiffness matrix A
%
%
% OUTPUT
% ------
% u_h               (nDof x N+1 double) Dof-weights of the entire history of
%                   the numerical solution
% (optional)
% data              (struct)
%   .dtHist         (2 x N+1 double) First row contains the history of the
%                   time steps, while the second indicates its acceptance
%                   (1) or rejection (0)
%   .lambda         (1 x 1 double) Estimation of the spectral radius of the
%                   differential problem (maximum time step for the Euler
%                   explicit scheme is then 2/lambda)
%   (if output.lmA_flag == 1 or output.smA_flag == 1 resp.)
%   .lmA            (1 x 1 double) Eigenvalue of largest magnitude of A
%   .smA            (1 x 1 double) Eigenvalue of smallest magnitude of A
% mesh              (see mesh_2d.m
% space             (see space_2d.m)
%
%
% -------------------------------------------------------------------------
%   Authors: Christoph Jaeggli, Julien Straubhaar and Philippe Rendard
%   Year: 2015
%   Institut: University of Neuchatel
%
%   This program is free software, you can redistribute it and/or modify
%   it.
%
%   Copyright (C) 2015 Christoph Jaeggli

%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


% Extract the fields from the data structures into local variables
prob_names = fieldnames(problem);
for iField = 1:numel(prob_names)
    eval( [ prob_names{iField}, ' = problem.(prob_names{iField});' ]);
end
meth_names = fieldnames(method);
for iField = 1:numel(meth_names)
    eval( [ meth_names{iField}, ' = method.(meth_names{iField});' ]);
end

% Construct mesh and space structures
mesh = mesh_2d(geo,nElDir,nQRDir);
space = space_2d(mesh,rBFDir);

% Assemble matrices
op.A = opGradUGradV(space, space, mesh, mu, geo);
op.M = opUV(  space, space, mesh, sigma, geo  );
Mlumped = sum(op.M,2);
if nargin(f) == 1
    op.f = opFV(space,mesh,f,geo);
end
lambda = abs(eigs( bsxfun(@times,-1./Mlumped,op.A),1,'lm' ));
if strcmp(timeScheme,'EE') && dt >= 2/lambda
    warning('Explicit time integration scheme %s is possibly unstable for selected time step',timeScheme)
end

% Compute initial solution weights
switch discIn
    case 'L2' % L2 projection
        u_h = op.M\opFV(space,mesh,u_0,geo);
    case 'IP' % Interpolation weights
        u_h = u_0(space_dofCoord(space,geo))';
    case 'VW' % Vector weights
        u_h = u_0;
    otherwise
        error('Unknown discretization scheme: %s',discIn)
end
if output.vtk_flag
    matlab2vtk( u_h, space, geo, 0, output.vtk_name)
end

% Boundary information
% Dirichlet
if ~isempty(BC.dir_sides)
    BC.dir = bndry_info(geo,mesh,space,BC.dir_sides,BC.dir_lim,...
        BC.dir_fun,'dof',1,'projfv',1,'mass',1);
else
    BC.dir.dof = [];
end
% Neuman
if ~isempty(BC.neum_sides)
    BC.neum = bndry_info(geo,mesh,space,BC.neum_sides,BC.neum_lim,...
        BC.neum_fun,'projfv',1);
else
    BC.neum = [];
end


% Time Integration
% Set up initial data
FE.geo = geo; FE.mesh = mesh; FE.space = space;
F_0 = rhsF( 0, f, op, FE ) + rhsDir( 0, op.A, FE.space, BC ) + rhsNeum( 0, FE.space, BC );
F_old = F_0;
RK.t = 0;
RK.dt = dt;
if strcmp(timeScheme,'RKC2')
    stage = 1;
    RK.s = max(ceil(sqrt(RK.dt*lambda/0.6)),2);
    RK.eta = eta;
else
    stage = 0;
end

% Time steps n = 1,2,...
n = 1; loop = 1;
dtHist = [];
while RK.t < Tmax && loop < lMax
    dtHist = [dtHist,[RK.dt;0]];
    
    % Time Integration
    eval( [ '[u_h(:,n+1), outData] = intHeat', timeScheme, '( u_h(:,n), F_old, f, FE, RK, op, BC);' ] );
    if adapt_dt.flag
        if strcmp(timeScheme,'EE')
            error('Time step selection scheme is not working for Euler explicit time integration')
        end
        adapt_dt.u_1p = outData.u_1p;
        adapt_dt.orderp = outData.orderp;
    end
    F_old = outData.F_1;
    
    % Error Estimation and Time Step Selection
    adapt_dt.u_0 = u_h(:,n); adapt_dt.u_1 = u_h(:,n+1);
    [RK, acc] = timeStepSel(RK,Tmax,adapt_dt);
    if stage == 1
        RK.s = max(ceil(sqrt(RK.dt*lambda/0.6)),2);
    end
    % Fprint Percentage of Elapsed Time
    if  output.screen_flag
        if loop > 1
            fprintf('\b\b\b\b\b\b')
        end
        fprintf(sprintf('%5.1f',round(1000*RK.t/Tmax)/10))
        fprintf('%%')
    end
    
    % Update
    loop = loop+1;
    dtHist(2,end) = acc;
    if acc
        if output.vtk_flag && mod(n,output.vtk_freq) == 0
            matlab2vtk( u_h(:,n+1), space, geo, n, output.vtk_name)
        end
        n = n + 1;
    end
end

data = struct('dtHist',dtHist,'lambda',lambda);

if isfield(output,'lmA_flag')
    if output.lmA_flag
        intDof = setdiff(1:space.nDof,BC.dir.dof);
        data.lmA = eigs(op.A(intDof,intDof),1,'lm');
    end
end
if isfield(output,'smA_flag')
    if output.smA_flag
        intDof = setdiff(1:space.nDof,BC.dir.dof);
        data.smA = eigs(op.A(intDof,intDof),1,'sm');
    end
end
if nargout == 2
    varargout{1} = data;
elseif nargout == 4
    varargout{1} = data;
    varargout{2} = mesh;
    varargout{3} = space;
end


end

%%%%%%%%%%%%%%%   Generate Boundary Geo/Mesh/Space    %%%%%%%%%%%%%%%%%
function [info] = bndry_info(geo,mesh,space,sides,limits,funs,varargin)
%BNDRY_INFO generates the complete boundary information consisting of
%   - geometry
%   - mesh
%   - space
%and optionally
%   - dof indices and sort map
%   - projection weights of the function ( (int fv) on the boundary)
%   - mass matrix on the boundary
%
%
% -------------------------------------------------------------------------
%
%   [info] = bndry_info(geo,mesh,space,sides,limits,funs);
%   [info] = bndry_info(geo,mesh,space,sides,limits,funs, 'name', flag, ...);
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% geo               (see geo_2d.m)
% mesh              (see mesh_2d.m)
% space             (see space_2d.m)
% sides             (1 x N integer) Side indices
% limits            (2 x N double) Relative limits (0<=a1 < a2 <=1)
% funs              (1 x N cell) Function handles
% (optional)
% Pairwise optional inputs are
%
%   'Name'   | flag |
% -------------------------------------------------------------------------
%   'dof'    | 0/1  |  - Compute all dof indices associated to the
%            |      |    boundary condition
%            |      |  - Gives the sorting map (output of unique(dofs)
%            |      |    for the dofs placed all in one big vector
%   -----------------------------------------------------------------------
%   'projfv' | 0/1  |  - Computes the weights of int(fv) on the boundary
%   -----------------------------------------------------------------------
%   'mass'   | 0/1  |  - Computes the mass matrix int(uv) on the boundary
%   -----------------------------------------------------------------------
%
%
% OUTPUT
% ------
% info              (sturct)
%   .geo            (1 x N cell) (see geo_2d:geo_boundary.m)
%   .mesh           (1 x N cell) (see mesh_2d:mesh_boundary.m)
%   .space          (1 x N cell) (see space_2d:space_boundary.m)
% (depending on the optional input)
%   .dof            (1 x nDofBndry integer) Dof indices of boundary
%                   condition
%   .sortMap        (1 x nDofBndry+M integer) Sort map for boundary weights
%   .projfv         (1 x N cell) Weights of int(fv) on the boundary
%   .mass           (1 x N cell) Mass matrices int(uv) on the boundary

info = [];

% Additional inputs for the output specification
nArgIn0 = -nargin('bndry_info')-1;
if mod(nargin-nArgIn0,2) == 1
    error('bndry_info: Additional input arguments have to come pairwise')
end
dofFlag = 0; projFlag = 0; massFlag = 0;
for iIn = 1:2:(nargin-nArgIn0-1)
    switch varargin{iIn}
        case 'dof'
            dofFlag = varargin{iIn+1};
            dof = [];
        case 'projfv'
            projFlag = varargin{iIn+1};
        case 'mass'
            massFlag = varargin{iIn+1};
    end
end

% Create complete boundary information
for iSide = 1:length(sides)
    info.geo{iSide} = geo_boundary( geo, sides(iSide) );
    info.mesh{iSide} = mesh_boundary( mesh, sides(iSide), limits(:,iSide) );
    info.space{iSide} = space_boundary( space, info.mesh{iSide}, sides(iSide) );
    
    if dofFlag
        dof = [dof; info.space{iSide}.dof'];
    end
    if projFlag
        if nargin(funs{iSide}) == 2
            fun = @(x)( funs{iSide}(x,0) );
%             warning('Boundary information evaluated at t=0')
        else
            fun = funs{iSide};
        end
        info.projfv{iSide} = opFV( info.space{iSide},...
            info.mesh{iSide}, fun, info.geo{iSide} );
    end
    if massFlag
        info.mass{iSide} = opUV( info.space{iSide}, info.space{iSide},...
            info.mesh{iSide}, ones(1,prod(info.mesh{iSide}.nElDir)), info.geo{iSide});
    end
end

if dofFlag
    [info.dof, info.sortMap, ~] = unique(dof);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%   Error Control and Time Step Selection   %%%%%%%%%%%%%%%%%%
function [RK, acc] = timeStepSel(RK,Tmax,adapt_dt)
%TIMESTEPSEL is an automatic time step selection using an error control
%algorithm for the acceptance or rejection of the previous time step based
%on a lower order approximation of the numerical solution.
%
%
% -------------------------------------------------------------------------
%
%   [RK, acc] = timeStepSel(RK,Tmax,adapt_dt);
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% RK, Tmax          (see main function input)
% adapt_dt          (struct)
%	.facmin         (1 x 1 double) Max decrease factor of time step
%   .facmax         (1 x 1 double) Max increase factor of time step
%	.fac            (1 x 1 double) Factor for time step selection
% 	.tol            (1 x 1 double) Tolerance for error control
%   .u_0            (nDof x 1 double) Solution at t_0
%   .u_1            (nDof x 1 double) Solution at t_1
%   .u_1p           (nDof x 1 double) Lower order approximation of the
%                   solution at at t_1
%   .orderp         (1 x 1 integer) Lower order
%
%
% OUTPUT
% ------
% RK                (sturct) As before, but with adapted time step .dt
% acc               (1 x 1 logical) Accept previous time step or not


if adapt_dt.flag
    rtol = adapt_dt.tol;
    atol = rtol*1e-3;
    facmax = adapt_dt.facmax;
    facmin = adapt_dt.facmin;
    fac = adapt_dt.fac;
    
    sc = sqrt(length(adapt_dt.u_0))*( atol + rtol*max(abs(adapt_dt.u_0),abs(adapt_dt.u_1)) );
    err = norm( (adapt_dt.u_1-adapt_dt.u_1p)./sc );
    
    if err == 0
        acc = 1;
        RK.t = RK.t + RK.dt;
        RK.dt = facmax*RK.dt;
    elseif err < 1
        acc = 1;
        RK.t = RK.t + RK.dt;
        RK.dt = RK.dt*min(facmax,max(facmin,fac*nthroot(1/err,adapt_dt.orderp+1)));
    else
        acc = 0;
        RK.dt = RK.dt*min(facmax,max(facmin,fac*nthroot(1/err,adapt_dt.orderp+1)));
    end
    RK.dt = min(RK.dt,Tmax-RK.t);
else
    acc = 1;
    RK.t = RK.t + RK.dt;
    RK.dt = min(RK.dt,Tmax-RK.t);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   OLD   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     switch timeScheme
%         case 'EE'
%             u_h(:,n+1) = intEE( u_h(:,n), F_old, FE, RK, op, BC);
%             F_old = opFV(FE.space,FE.mesh,@(x)(f(x,RK.t+RK.dt)),FE.geo) +...
%                 rhsDir( RK.t+RK.dt, op.A, FE.space, BC ) +...
%                 rhsNeum( RK.t+RK.dt, FE.space, BC );
%         case 'EI'
%             [u_h(:,n+1), ~] = intEI( u_h(:,n), f, FE, RK, op, BC);
%         case 'CN'
%             [u_h(:,n+1), outData] = intCN( u_h(:,n), F_old, f, FE, RK, op, BC);
%             if adapt_dt.flag
%                 adapt_dt.u_1p = intEE( u_h(:,n), F_old, FE, RK, op, BC);
%                 adapt_dt.orderp = 1;
%             end
%             F_old = outData.F_1;
%         case 'RKC2'
%             RK.s = max(ceil(sqrt(RK.dt*lambda/0.6)),2);
%             RK.eta = eta;
%             [u_h(:,n+1), outData] = intRKC2( u_h(:,n), F_old, f, FE, RK, op, BC);
%             if adapt_dt.flag
%                 adapt_dt.u_1p = intEE( u_h(:,n), F_old, FE, RK, op, BC);
%                 adapt_dt.orderp = 1;
%             end
%             F_old = outData.F_1;
%         case 'RIIA'
%             [u_h(:,n+1), outData] = intRIIA( u_h(:,n), F_old, f, FE, RK, op, BC);
%             if adapt_dt.flag
%                 adapt_dt.u_1p = outData.u_1p;
%                 adapt_dt.orderp = outData.orderp;
%             end
%             F_old = outData.F_1;
%         otherwise
%             error('Unknown time-scheme %s',timeScheme)
%     end
