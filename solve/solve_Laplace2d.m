function [ u_h, varargout] = solve_Laplace2d( problem, method )
%SOLVE_LAPLACE2D solves a 2-d Laplace problem by means of the Finite
%Element Method. The differential problem reads:
%
%       - div( mu*grad( u ) ) = f   in \Omega
%                           u = g   on \Gamma_D
%                    mu*du/dn = h   on \Gamma_N
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
%   u_h = solve_Laplace2D(problem, method);
%   [u_h, mesh, space] = solve_Laplace2D(problem, method);
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
%   .f              (function handle) or (1 x #mesh elements double)
%                   Forcing term
% method
%   .nElDir         (2 x 1 double) Number of mesh elements in each
%                   parametric direction
%   .nQRDir         (2 x 1 double) Number of quadrature nodes per element
%                   in each parametric direction
%   .rBFDir         (2 x 1 double) Order of the piecewise polynomial basis
%                   functions in each parametric direction
%
%
% OUTPUT
% ------
% u_h               ( nDof x 1 double) The computed degrees of freedom
% (optional)
% mesh              (see output of mesh_2d.m)
% space             (see output of space_2d.m)
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
BC_names = fieldnames(BC);
for iField = 1:numel(BC_names)
    eval( [ BC_names{iField}, ' = BC.(BC_names{iField});' ]);
end
meth_names = fieldnames(method);
for iField = 1:numel(meth_names)
    eval( [ meth_names{iField}, ' = method.(meth_names{iField});' ]);
end

% Construct mesh and space structures
mesh = mesh_2d(geo,nElDir,nQRDir);
space = space_2d(mesh,rBFDir);
if nargout == 3
    varargout{1} = mesh;
    varargout{2} = space;
end

% Assemble matrices
A = opGradUGradV(space, space, mesh, mu, geo);
rhs = opFV(space,mesh,f,geo);

% Boundary conditions
% Neumann
for iNeum = 1:length(neum_sides)
    geo_bndry = geo_boundary( geo, neum_sides(iNeum) );
    mesh_bndry = mesh_boundary( mesh, neum_sides(iNeum), neum_lim(:,iNeum) );
    space_bndry = space_boundary( space, mesh_bndry, neum_sides(iNeum));
    
    rhs(space_bndry.dof) = rhs(space_bndry.dof) + opFV( space_bndry, mesh_bndry, neum_fun{iNeum}, geo_bndry );
end
% Dirichlet
dirDof = [];
dirWeights = [];
for iDir = 1:length(dir_sides)
    geo_bndry = geo_boundary( geo, dir_sides(iDir) );
    mesh_bndry = mesh_boundary( mesh, dir_sides(iDir), dir_lim(:,iDir) );
    space_bndry = space_boundary( space, mesh_bndry, dir_sides(iDir));
    
    M_bndry = opUV(space_bndry,space_bndry,mesh_bndry,ones(1,prod(mesh_bndry.nElDir)),geo_bndry);
    l2Proj_bndry = M_bndry\opFV( space_bndry, mesh_bndry, dir_fun{iDir}, geo_bndry );
    
    dirDof = [dirDof; space_bndry.dof'];
    dirWeights = [dirWeights; l2Proj_bndry];
end

% Priority to DIRICHLET BC (on the eventual intersection)
[dirDof, indRestr, ~] = unique(dirDof);
dirWeights = dirWeights(indRestr);

intDof = setdiff(1:space.nDof,dirDof);
rhs(intDof) = rhs(intDof) - A(intDof,dirDof)*dirWeights;

A(dirDof,:) = 0;
A(:,dirDof) = 0;
A(dirDof,dirDof) = eye(length(dirDof));
rhs(dirDof) = dirWeights;

% Solve the linear problem
u_h = A\rhs;

end

