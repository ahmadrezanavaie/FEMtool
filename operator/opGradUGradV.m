function [ varargout ] = opGradUGradV( spU, spV, mesh, mu, geo )
%OPGRADUGRADV assembles the stiffness matrix (A)_(i,j) where 
%
%   (A)_(i,j) = \int (mu*grad u_j.grad v_i)dx.
%
% -------------------------------------------------------------------------
%
%   A = opGradUGradV( spU, spV, mesh, mu, geo );
%   [rows, cols, vals] = opGradUGradV( spU, spV, mesh, mu, geo );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% spU,spV           (see output of space_2d.m/space_3d.m)
% mesh              (see output of mesh_2d.m/mesh_3d.m)
% mu                (function handle) or (1 x mesh.nElTot double)
%                   representing the diffusion coefficient
% geo               (see output of geo_2d.m/geo_3d.m --> or define your own)
%
%
% OUTPUT
% ------
% A                 (spU.nDof x spV.nDof double) Assembled stiffness matrix
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

% REGROUPE MU, QUADRATURE WEIGHTS and DETERMINANT of the JACOBIAN matrix
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

if isa(mu, 'function_handle')
    mu = reshape(mu(geo.map(pQN)),nQNEl,nEl);
elseif size(mu,1) == prod(mesh.nElDir) && size(mu,2) == 1 ||...
        size(mu,1) == 1 && size(mu,2) == prod(mesh.nElDir)
    mu = repmat(reshape(mu,1,prod(mesh.nElDir)),prod(mesh.nQNElDir),1);
else
    error('Bad definition of mu (function handle or vector of length #elements)')
end

% REGROUP JACINV such that 
% size(...) = [#dim, #dim, #QN per el, #Elements]
JacInv = reshape(geo.JacInvT(geo.map(pQN)),nDim,nDim,nQNEl,nEl);
clear pQN qN qW;

% REGROUP GRADIENTS of basis functions such that
% size(...) = [#dim, #shape function per el, #QN per el, #Elements]
%
% Remember that in space_2d and space_3d the basis functions are defined
% as:
%
% 2d:
%   f(x) = p1(x1)*p2(x2)
%
% 3d:
%   f(x) = p1(x1)*p2(x2)*p3(x3),
%
% therefore the gradients are
%
% 2d:
%                 / p1'(x1)*p2(x2) \
%   grad(f(x)) = |                  |
%                 \ p1(x1)*p2'(x2) /
%
% 3d:
%                 / p1'(x1)*p2(x2)*p3(x3) \
%   grad(f(x)) = |  p1(x1)*p2'(x2)*p3(x3)  |
%                 \ p1(x1)*p2(x2)*p3'(x3) /
nBfElU = prod(spU.nDofElDir);
nBfElV = prod(spV.nDofElDir);
gradU = zeros(nDim, nBfElU, nQNEl, nEl);
gradV = zeros(nDim, nBfElV, nQNEl, nEl);
for iDim = 1:nDim
    
    indReshDerU = ones(1,3*nDim);
    indReshDerU([iDim, nDim+iDim, 2*nDim+iDim]) = [spU.nDofElDir(iDim),mesh.nQNElDir(iDim),mesh.nElDir(iDim)]; 
    valU = reshape(spU.bFDerDir{iDim},indReshDerU);
    
    indReshDerV = ones(1,3*nDim);
    indReshDerV([iDim, nDim+iDim, 2*nDim+iDim]) = [spV.nDofElDir(iDim),mesh.nQNElDir(iDim),mesh.nElDir(iDim)]; 
    valV = reshape(spV.bFDerDir{iDim},indReshDerV);
    
    for iNDim = setdiff(1:nDim,iDim)
        indReshU = ones(1,3*nDim);
        indReshU([iNDim, nDim+iNDim, 2*nDim+iNDim]) = [spU.nDofElDir(iNDim),mesh.nQNElDir(iNDim),mesh.nElDir(iNDim)]; 
        valU = bsxfun(@times, valU, reshape(spU.bFDir{iNDim},indReshU));
        
        indReshV = ones(1,3*nDim);
        indReshV([iNDim, nDim+iNDim, 2*nDim+iNDim]) = [spV.nDofElDir(iNDim),mesh.nQNElDir(iNDim),mesh.nElDir(iNDim)]; 
        valV = bsxfun(@times, valV, reshape(spV.bFDir{iNDim},indReshV));
    end
    
    gradU(iDim,:,:,:) = reshape(valU,nBfElU,nQNEl,nEl);
    gradV(iDim,:,:,:) = reshape(valV,nBfElV,nQNEl,nEl);
end

% LOOP through the elements of the mesh
rows = zeros(prod(mesh.nElDir)*nBfElU*nBfElV,1);
cols = zeros(prod(mesh.nElDir)*nBfElU*nBfElV,1);
vals = zeros(prod(mesh.nElDir)*nBfElU*nBfElV,1);

loc2GlobU = spU.loc2Glob;
loc2GlobV = spV.loc2Glob;

nCounter = 0;
for iEl = 1:prod(mesh.nElDir)
    gradU_iEl = reshape(gradU(:,:,:,iEl),1,nDim,nBfElU,nQNEl);
    gradV_iEl = reshape(gradV(:,:,:,iEl),1,nDim,nBfElV,nQNEl);
    JacInv_iEl = reshape(JacInv(:,:,:,iEl),nDim,nDim,1,nQNEl);
    
    % See remark in geo_2d.m or geo_3d.m
    JacInvGradU_iEl = reshape(sum(bsxfun(@times,JacInv_iEl,gradU_iEl),2),nDim,nBfElU,1,nQNEl);
    JacInvGradV_iEl = reshape(sum(bsxfun(@times,JacInv_iEl,gradV_iEl),2),nDim,1,nBfElV,nQNEl);
    
    gradUGradV_iEl = bsxfun(@times,JacInvGradU_iEl,JacInvGradV_iEl);
    muGradUGradV_iEl = bsxfun(@times,gradUGradV_iEl,reshape(mu(:,iEl),1,1,1,nQNEl));
    qWDetJac_iEl = reshape(weights(:,iEl).*abs(detJac(:,iEl)),1,1,1,nQNEl);
    
    vals_iEl = sum(sum(bsxfun(@times,muGradUGradV_iEl,qWDetJac_iEl),4),1);
    
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
    error ('opGradUGradV: Wrong number of output arguments')
end

end

