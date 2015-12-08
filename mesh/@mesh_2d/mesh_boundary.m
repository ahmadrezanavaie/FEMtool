function [ mesh_bndry ] = mesh_boundary( mesh, iSide, varargin )
%MESH_BOUNDARY is a method of the class mesh_2d that generates a 1d mesh of
%a part of the boundary of a 2d parametric domain.

%     2D Parametric Mesh
%
%            s4
%      ___ ___ ___ ___ 
%     |   |   |   |   |                     1D Bndry Mesh
%     |___|___|___|___|
%     |   |   |   |   |          --->     .___.___.___.
%  s1 |___|___|___|___| s2
%     |   |   |   |   |
% p2  |___|___|___|___|
%  ^
%  |-> p1    s3
%
%
% -------------------------------------------------------------------------
%
%   mesh_bndry = mesh_boundary( mesh, iSide );
%   mesh_bndry = mesh_boundary( mesh, iSide, limits );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% mesh              (see output of mesh_2d.m)
% iSide             (1 x 1 integer) Indicating the identification number of
%                   the boundary space 
% (optional)
% limits            (2 x 1 double) Two numbers 0 <= c1 < c2 <= 1 indicating
%                   the subinterval of the boundary of which the mesh
%                   should be extracted (by default c1 = 0, c2 = 1)
%
%
% OUTPUT
% ------
% mesh_bndry        (1 x 1 struct)
%   .nElDir         (1 x 1 double) Number of elements
%   .nQNElDir       (2 x 1 double) Number of quadrature nodes per element 
%   .qNDir          (1 x 1 cell-array)
%       .qNDir{1}   (nQNElDir x nElDir double) Contains the quadr. nodes
%   .qWDir          (1 x 2 cell-array) 
%       .qWDir{1}   (nQNElDir x nElDir double) Contains the quadr. weights
%   .indBndryEl     (1 x nElDir double) Contains the subset of element
%                   identification numbers of the boundary elements
%                   contained in the 1d mesh
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


nArgIn0 = -nargin('mesh_boundary')-1;

if nargin-nArgIn0 == 1
    limits = varargin{1};
elseif nargin-nArgIn0 == 0
    limits = [0;1];
else
    error('mesh_2d/mesh_boundary: Wrong number of input arguments')
end

if iSide == 1 || iSide == 2
    iDir = 2;
elseif iSide == 3 || iSide == 4
    iDir = 1;
else
    error('mesh_2d/mesh_boundary: Wrong side number index')
end

nodes = linspace(0,1,mesh.nElDir(iDir)+1);
[~, iNodeMin] = min(abs(nodes-limits(1)));
[~, iNodeMax] = min(abs(nodes-limits(2)));
bndryEl = iNodeMin:(iNodeMax-1);

if isempty(bndryEl)
    error('mesh_2d/mesh_boundary: Defined boundary subdomain is empty')
end

mesh_bndry.nElDir = length(bndryEl);
mesh_bndry.nQNElDir = mesh.nQNElDir(iDir);

mesh_bndry.qNDir{1} = mesh.qNDir{iDir}(:,bndryEl);
mesh_bndry.qWDir{1} = mesh.qWDir{iDir}(:,bndryEl);

mesh_bndry.indBndryEl = bndryEl;

end

