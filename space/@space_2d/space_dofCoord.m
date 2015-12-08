function [ dofCoord ] = space_dofCoord( space, geo )
%SPACE_DOFCOORD Is a method of the class space_2d and computes
%coordinates of the nodes associated to the degrees of freedom in the
%physical domain.
%
%
% -------------------------------------------------------------------------
%
%   dofCoord = space_dofCoord( space, geo );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% space         (see output space_2d.m)
% geo           (see output geo_2d.m)
%
%
% OUTPUT
% ------
% dofCoord      (2 x space.nDof double) Physical coordinates of the dof 
%               nodes
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

origin = geo.paramDom.origin;
sizeDir = geo.paramDom.sizeDir;

xx = origin(1) + linspace(0,1,space.nDofDir(1))*sizeDir(1);
yy = origin(2) + linspace(0,1,space.nDofDir(2))*sizeDir(2);

xx = reshape(repmat(xx',1,space.nDofDir(2)),1,space.nDof);
yy = reshape(repmat(yy,space.nDofDir(1),1),1,space.nDof);

dofCoord = geo.map([xx;yy]);


end

