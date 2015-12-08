function [ space_bndry ] = space_boundary( space, mesh_bndry, iSide )
%SPACE_BOUNDARY is a method of the class space_2d that generates the finite
%element function space on the 1d parametric boundary domain. Depending on
%the order of the piecewise polynomial basis functions some
%degrees of freedom may be added to the vertices in the mesh. The local
%degrees of freedom are enumerated as follows:
%
%         order = 1                   order = 2
%         ---------                   ---------
%
%       1___________2               1_____2_____3
%
%
%         order = 3                   order = 4
%         ---------                   ---------
%
%       1___2___3___4               1__2__3__4__5
%
%
% -------------------------------------------------------------------------
%
%   space_bndry = space_boundary( space, mesh, iSide );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% space             (see output of space_2d.m)
% mesh_bndry        (see output of mesh_2d/mesh_boundary)
% iSide             
%
%
% OUTPUT
% ------
% space             (1 x 1 struct)
%   .nDim           (1 x 1 double) Dimension of the boundary domain
%   .nDofElDir      (1 x 1 double) Number of degrees of freedom per element
%   .nDof           (1 x 1 double) Number of degrees of freedom
%   .bFDir          (1 x 1 cell-array)
%       .bFDir{1}   (nDofElDir x mesh_bndry.nQNElDir x mesh_bndry.nElDir
%                   double ) For each element in it contains all the
%                   non-zero basis functions evaluated at all the
%                   quadrature nodes  
%   .bFDerDir       (1 x 1 cell-array)
%       .bFDerDir{1}(nDofElDir x mesh_bndry.nQNElDir x mesh_bndry.nElDir
%                   double ) For each element in it contains the
%                   derivatives of all the non-zero basis functions
%                   evaluated at all the quadrature nodes  
%   .loc2Glob       (nDofElDir x mesh_bndry.nElDir double) For each
%                   element this local-to-global map indicates the global
%                   dof indices (in the boundary space)
%   .dof            (1 x nDof double) A list of all dof identification
%                   number in the 2 dimensional space associated to the
%                   boundary space
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
%   Copyright (C) 2015 Christoph Jäggli

%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

if iSide == 1
    iDir = 2;
    firstDof = ((mesh_bndry.indBndryEl(1)-1)*space.rBFDir(2))*space.nDofDir(1) + 1;
    incr = space.nDofDir(1);
elseif iSide == 2
    iDir = 2;
    firstDof = ((mesh_bndry.indBndryEl(1)-1)*space.rBFDir(1) + 1)*space.nDofDir(1);
    incr = space.nDofDir(1);
elseif iSide == 3
    iDir = 1;
    firstDof = (mesh_bndry.indBndryEl(1)-1)*space.rBFDir(1) + 1;
    incr = 1;
elseif iSide == 4
    iDir = 1;
    firstDof = (space.nDofDir(2)-1)*space.nDofDir(1) + (mesh_bndry.indBndryEl(1)-1)*space.rBFDir(1) + 1;
    incr = 1;
else
    error('space_2d/space_boundary: Wrong side number index')
end

space_bndry.nDim = space.nDim-1;
space_bndry.nDofElDir = space.nDofElDir(iDir);
space_bndry.nDof = mesh_bndry.nElDir*space.rBFDir(iDir) + 1;
space_bndry.bFDir{1} = space.bFDir{iDir}(:,:,mesh_bndry.indBndryEl);
space_bndry.bFDerDir{1} = space.bFDerDir{iDir}(:,:,mesh_bndry.indBndryEl);

indBf = bsxfun(@plus,(1:space_bndry.nDofElDir)',0:(space_bndry.nDofElDir-1):((space_bndry.nDofElDir-1)*(mesh_bndry.nElDir-1)));
indBf = reshape(indBf,space_bndry.nDofElDir,mesh_bndry.nElDir,1);
space_bndry.loc2Glob = indBf;

space_bndry.dof = firstDof + incr*(0:(space_bndry.nDof-1));

end

