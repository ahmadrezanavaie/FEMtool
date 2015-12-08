function [ varargout ] = opUV(  spU, spV, mesh, sigma, geo  )
%OPUV assembles the mass matrix (M)_(i,j) where 
%
%   (M)_(i,j) = \int (sigma u_jv_i)dx.
%
% -------------------------------------------------------------------------
%
%   M = opUV( spU, spV, mesh, sigma, geo );
%   [rows, cols, vals] = opUV( spU, spV, mesh, sigma, geo );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% spU,spV           (see output of space_2d.m/space_3d.m)
% mesh              (see output of mesh_2d.m/mesh_3d.m)
% sigma             (function handle) or (1 x mesh.nElTot double)
%                   representing the reaction coefficient
% geo               (see output of geo_2d.m/geo_3d.m --> or define your own)
%
%
% OUTPUT
% ------
% M                 (spU.nDof x spV.nDof double) Assembled mass matrix
%                   (sparse)
% rows              (nnz x 1 double) (repeated) Row indices of the nonzero
%                   entries
% cols              (nnz x 1 double) (repeated) Columnd indices of the
%                   nonzero entries
% vals              (nnz x 1 double) Values corresponding to the row/column
%                   indices
%
%
% Remark:
% -------
% 'help sparse' gives:
% 'S = sparse(i,j,s,m,n) uses vectors i, j, and s to generate an m-by-n
% sparse matrix such that S(i(k),j(k)) = s(k). Vectors i, j, and s are all
% the same length. ANY ELEMENTS OF s THAT HAVE DUPLICATE VALUES OF i AND j
% ARE ADDED TOGETHER.
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

if spU.nDim ~= spV.nDim
    error('Physical dimensions of the domains of the function spaces have to be equal')
else
    nDim = spU.nDim;
end

nEl = prod(mesh.nElDir);
nQNEl = prod(mesh.nQNElDir);

% REGROUPE SIGMA, QUADRATURE WEIGHTS and DETERMINANT of the JACOBIAN matrix
% such that
% size(...) = [#QN per El (in TOTAL not per direction) x #Elements (total)]
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

if isa(sigma, 'function_handle')
    sigma = reshape(sigma(geo.map(pQN)),nQNEl,nEl);
elseif size(sigma,1) == prod(mesh.nElDir) && size(sigma,2) == 1 ||...
        size(sigma,1) == 1 && size(sigma,2) == prod(mesh.nElDir)
    sigma = repmat(reshape(sigma,1,prod(mesh.nElDir)),prod(mesh.nQNElDir),1);
else
    error('Bad definition of mu (function handle or vector of length #elements)')
end


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
nBfElU = prod(spU.nDofElDir);
nBfElV = prod(spV.nDofElDir);

indReshU = ones(1,3*nDim);
indReshU([1, nDim+1, 2*nDim+1]) = [spU.nDofElDir(1),mesh.nQNElDir(1),mesh.nElDir(1)]; 
bFunU = reshape(spU.bFDir{1},indReshU);

indReshV = ones(1,3*nDim);
indReshV([1, nDim+1, 2*nDim+1]) = [spV.nDofElDir(1),mesh.nQNElDir(1),mesh.nElDir(1)]; 
bFunV = reshape(spV.bFDir{1},indReshV);

for iDim = 2:nDim
    indReshU = ones(1,3*nDim);
    indReshU([iDim, nDim+iDim, 2*nDim+iDim]) = [spU.nDofElDir(iDim),mesh.nQNElDir(iDim),mesh.nElDir(iDim)]; 
    bFunU = bsxfun(@times, bFunU, reshape(spU.bFDir{iDim},indReshU));
    
    indReshV = ones(1,3*nDim);
    indReshV([iDim, nDim+iDim, 2*nDim+iDim]) = [spV.nDofElDir(iDim),mesh.nQNElDir(iDim),mesh.nElDir(iDim)]; 
    bFunV = bsxfun(@times, bFunV, reshape(spV.bFDir{iDim},indReshV));
end
bFunU = reshape(bFunU,nBfElU,nQNEl,nEl);
bFunV = reshape(bFunV,nBfElV,nQNEl,nEl);


% LOOP through the elements of the mesh
rows = zeros(prod(mesh.nElDir)*nBfElU*nBfElV,1);
cols = zeros(prod(mesh.nElDir)*nBfElU*nBfElV,1);
vals = zeros(prod(mesh.nElDir)*nBfElU*nBfElV,1);

loc2GlobU = spU.loc2Glob;
loc2GlobV = spV.loc2Glob;

nCounter = 0;
for iEl = 1:prod(mesh.nElDir)
    bFunU_iEl = reshape(bFunU(:,:,iEl),nBfElU,1,nQNEl);
    bFunV_iEl = reshape(bFunV(:,:,iEl),1,nBfElV,nQNEl);
    
    bFunUBFunV_iEl = bsxfun(@times,bFunU_iEl,bFunV_iEl);
    sigmaBFunUBFunV_iEl = bsxfun(@times,bFunUBFunV_iEl,reshape(sigma(:,iEl),1,1,nQNEl));
    qWDetJac_iEl = reshape(weights(:,iEl).*abs(detJac(:,iEl)),1,1,nQNEl);
    
    vals_iEl = sum(bsxfun(@times,sigmaBFunUBFunV_iEl,qWDetJac_iEl),3);
    
    rows_iEl = bsxfun(@times,loc2GlobV(:,iEl),ones(1,nBfElU));
    cols_iEl = bsxfun(@times,loc2GlobU(:,iEl),ones(1,nBfElV));
    rows(nCounter+(1:numel(rows_iEl))) = reshape(rows_iEl',numel(rows_iEl),1);
    cols(nCounter+(1:numel(cols_iEl))) = reshape(cols_iEl,numel(cols_iEl),1);
    vals(nCounter+(1:numel(rows_iEl))) = reshape(vals_iEl,numel(rows_iEl),1);
    nCounter = nCounter+numel(rows_iEl);
end

if nargout == 1
    varargout{1} = sparse(rows(1:nCounter),cols(1:nCounter),vals(1:nCounter),spV.nDof,spU.nDof);
elseif (nargout == 3)
    varargout{1} = rows(1:nCounter);
    varargout{2} = cols(1:nCounter);
    varargout{3} = vals(1:nCounter);
else
    error ('opUV: Wrong number of output arguments')
end

end

