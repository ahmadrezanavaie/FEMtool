function [ rhs ] = rhsNeum( t, space, BC )
%RHSNEUM computes the right hand side weights of the Neumann boundary
%conditions defined in BC. If the conditions are stationary, previously
%computed weights are re-used.
%
%
% -------------------------------------------------------------------------
%
%   rhs = rhsNeum( t, space, BC );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% t                 (1 x 1 double) Time t
% space             (see output space_2d.m/space_3d.m)
% BC                (struct)
%	.neum_sides     (1 x nNeum double) Neumann boundary indices
%   (if .neum_sides non-empty)
%   .neum_lim       (2 x nNeum double) Neumann limits
%   .neum_fun       (1 x nNeum cell array) Neumann function handles
%   .neum           (see solve_Heat2d:bndry_info)
%
%
% OUTPUT
% ------
% rhs               ( nDof x 1 double) Right hand side weights associated
%                   to the Neumann boundary conditions
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
for iNeum = 1:length(BC.neum_sides)
    if nargin(BC.neum_fun{iNeum}) == 2
        fun = @(x)( BC.neum_fun{iNeum}(x,t) );
        rhs(BC.neum.space{iNeum}.dof) = rhs(BC.neum.space{iNeum}.dof) +...
            opFV( BC.neum.space{iNeum}, BC.neum.mesh{iNeum}, fun, BC.neum.geo{iNeum} );
    else
        rhs(BC.neum.space{iNeum}.dof) = rhs(BC.neum.space{iNeum}.dof) + BC.neum.projfv{iNeum};
    end
end

end

