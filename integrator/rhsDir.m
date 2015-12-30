function [ rhs, varargout ] = rhsDir( t, op, space, BC )
%RHSDIR computes the right hand side weights of the Dirichlet boundary
%conditions defined in BC. If the conditions are stationary, previously
%computed weights are re-used.
%
%
% -------------------------------------------------------------------------
%
%   rhs = rhsDir( t, A, space, BC );
%   [ rhs, wDir ] = rhsDir( t, A, space, BC );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% t                 (1 x 1 double) Time t
% op                (struct)
%   .A              (nDof x nDof double) Stiffness matrix
%   .M              (nDof x nDof double) Mass matrix
% space             (see output space_2d.m/space_3d.m)
% BC                (struct)
%	.dir_sides      (1 x nDir double) Dirichlet boundary indices
%	.dir_lim        (2 x nDir double) Dirichlet limits
%	.dir_fun        (1 x nDir cell array) Dirichlet function handles
%   .dir            (see solve_Heat2d:bndry_info)
%
%
% OUTPUT
% ------
% rhs               ( nDof x 1 double) Right hand side weights associated
%                   to the Dirichlet boundary conditions
% (optional)
% wDir              ( nDirDof x 1 double) Projection weights of the
%                   Dirichlet boundary condition functions
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

rhs = zeros(space.nDof,1);
wDir = [];
if ~isempty(BC.dir_sides)
    for iSide = 1:length(BC.dir_sides)
        if nargin(BC.dir_fun{iSide}) == 2
            fun = @(x)( BC.dir_fun{iSide}(x,t) );
            proj = opFV( BC.dir.space{iSide},...
                BC.dir.mesh{iSide}, fun, BC.dir.geo{iSide} );
        else
            proj = BC.dir.projfv{iSide};
        end
        weights = BC.dir.mass{iSide}\proj;
        wDir = [wDir; weights];
    end
    wDir = wDir(BC.dir.sortMap);

    dirDof = BC.dir.dof;
    intDof = setdiff(1:space.nDof,dirDof);


    rhs(intDof) = -op.A(intDof,dirDof)*wDir;
end

if nargout == 2
    varargout{1} = wDir;
end

end

