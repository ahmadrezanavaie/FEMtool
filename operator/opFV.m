function [ rhs ] = opFV( spV, mesh, f, geo )
%OPFV assembles the right-hand side vector (rhs)_(i) where 
%
%   (rhs)_(i) = \int (f v_i)dx.
%
% -------------------------------------------------------------------------
%
%   rhs = opFV( spV, mesh, f, geo );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% spV               (see output of space_2d.m/space_3d.m)
% mesh              (see output of mesh_2d.m/mesh_3d.m)
% f                 (function handle) or (1 x mesh.nElTot double)
%                   representing the source function
% geo               (see output of geo_2d.m/geo_3d.m --> or define your own)
%
%
% OUTPUT
% ------
% rhs               (spV.nDof x 1 double) Assembled rhs vector
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

nDim = spV.nDim;
nEl = prod(mesh.nElDir);
nQNEl = prod(mesh.nQNElDir);

% REGROUPE F, QUADRATURE WEIGHTS and DETERMINANT of the JACOBIAN matrix
% such that
% size(...) = [#QN per El x #Elements]
pQN = zeros(nDim,nQNEl*nEl);
qWDir = zeros(nDim,nQNEl*nEl);
for iDim = 1:nDim
    indResh = ones(1,2*nDim);
    indResh([iDim, nDim + iDim]) = [mesh.nQNElDir(iDim), mesh.nElDir(iDim)];
    indRep = [mesh.nQNElDir', mesh.nElDir'];
    indRep([iDim, nDim + iDim]) = [1, 1];
    
    qN = reshape(mesh.qNDir{iDim},indResh);
    qN = repmat(qN,indRep);
    pQN(iDim,:) = reshape(qN,1,numel(qN));
    
    qW = reshape(mesh.qWDir{iDim},indResh);
    qW = repmat(qW,indRep);
    qWDir(iDim,:) = reshape(qW,1,numel(qW));
end
weights = reshape(prod(qWDir,1),nQNEl,nEl);
detJac = reshape(geo.detJac(pQN),nQNEl,nEl);

if isa(f, 'function_handle')
    f = reshape(f(geo.map(pQN)),nQNEl,nEl);
elseif size(f,1) == prod(mesh.nElDir) && size(f,2) == 1 ||...
        size(f,1) == 1 && size(f,2) == prod(mesh.nElDir)
    f = repmat(reshape(f,1,prod(mesh.nElDir)),prod(mesh.nQNElDir),1);
else
    error('Bad definition of f (function handle or vector of length #elements)')
end
clear pQN qN qW;


% REGROUP FUNCTIONS such that
% size(...) = [#shape function per el, #QN per el, #Elements]
%
% Remember that in space_2d and space_3d the basis functions are defined
% as:
%
% 2d:
%   f(x) = p1(x1)*p2(x2)
%
% 3d:
%   f(x) = p1(x1)*p2(x2)*p3(x3).
nBfElV = prod(spV.nDofElDir);

indReshV = ones(1,3*nDim);
indReshV([1, nDim+1, 2*nDim+1]) = [spV.nDofElDir(1),mesh.nQNElDir(1),mesh.nElDir(1)]; 
bFunV = reshape(spV.bFDir{1},indReshV);

for iDim = 2:nDim
    indReshV = ones(1,3*nDim);
    indReshV([iDim, nDim+iDim, 2*nDim+iDim]) = [spV.nDofElDir(iDim),mesh.nQNElDir(iDim),mesh.nElDir(iDim)]; 
    bFunV = bsxfun(@times, bFunV, reshape(spV.bFDir{iDim},indReshV));
end
bFunV = reshape(bFunV,nBfElV,nQNEl,nEl);


% LOOP through the elements of the mesh
loc2GlobV = spV.loc2Glob;

rhs = zeros(spV.nDof,1);
for iEl = 1:prod(mesh.nElDir)
    bFunV_iEl = bFunV(:,:,iEl);
    fBFunV_iEl = bsxfun(@times,bFunV_iEl,reshape(f(:,iEl),1,nQNEl));
    qWDetJac_iEl = reshape(weights(:,iEl).*detJac(:,iEl),1,nQNEl);
    
    vals_iEl = sum(bsxfun(@times,fBFunV_iEl,qWDetJac_iEl),2);
    
    rhs(loc2GlobV(:,iEl)) = rhs(loc2GlobV(:,iEl)) + vals_iEl;
end

end

