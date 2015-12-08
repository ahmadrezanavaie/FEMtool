function [ varargout ] = spErrSol( space, mesh, geo, u_h, uEx, varargin  )
%SPERRSOL evaluates the error between the numerical solution u_h and a
%given exact soluton uEx in L^2 and H^1 norm in a given FE-space.
%
%
% -------------------------------------------------------------------------
%
%   errL2 = spErrSol( space, mesh, geo, u_h, uEx );
%   [ errL2, errH1 ] = spErrSol( space, mesh, geo, u_h, uEx, uExGrad );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% space             (see output of space_2d.m/space_3d.m)
% mesh              (see output of mesh_2d.m/mesh_3d.m)
% geo               (see output of geo_2d.m/geo_3d.m --> or define your own)
% u_h               (space.nDof x 1 double) Vector of dof-weights of the
%                   numerical solution
% uEx               (function handle) Representing the exact solution
% (optional)
% uExGrad           (function handle) Representing the gradient of the
%                   exact solution
%
%
% OUTPUT
% ------
% errL2             (1 x 1 double) L^2 error between the numerical and the
%                   exact solution
% (optional)
% errH1             (1 x 1 double) H^1 error between the numerical and the
%                   exact solution
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

nArgIn0 = -nargin('spErrSol')-1;

if nargin-nArgIn0 == 1
    compH1 = 1;
    uExGrad = varargin{1};
else
    compH1 = 0;
end

nDim = space.nDim;
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
uEx = reshape(uEx(geo.map(pQN)),nQNEl,nEl);

if compH1
    uExGrad = reshape(uExGrad(geo.map(pQN)),nDim,nQNEl,nEl);
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
nBfEl = prod(space.nDofElDir);

indResh = ones(1,3*nDim);
indResh([1, nDim+1, 2*nDim+1]) = [space.nDofElDir(1),mesh.nQNElDir(1),mesh.nElDir(1)]; 
bFun = reshape(space.bFDir{1},indResh);

for iDim = 2:nDim
    indResh = ones(1,3*nDim);
    indResh([iDim, nDim+iDim, 2*nDim+iDim]) = [space.nDofElDir(iDim),mesh.nQNElDir(iDim),mesh.nElDir(iDim)]; 
    bFun = bsxfun(@times, bFun, reshape(space.bFDir{iDim},indResh));
end
bFun = reshape(bFun,nBfEl,nQNEl,nEl);


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
nBfEl = prod(space.nDofElDir);
bFGrad = zeros(nDim, nBfEl, nQNEl, nEl);
for iDim = 1:nDim
    
    indReshDer = ones(1,3*nDim);
    indReshDer([iDim, nDim+iDim, 2*nDim+iDim]) = [space.nDofElDir(iDim),mesh.nQNElDir(iDim),mesh.nElDir(iDim)]; 
    val = reshape(space.bFDerDir{iDim},indReshDer);
    
    for iNDim = setdiff(1:nDim,iDim)
        indResh = ones(1,3*nDim);
        indResh([iNDim, nDim+iNDim, 2*nDim+iNDim]) = [space.nDofElDir(iNDim),mesh.nQNElDir(iNDim),mesh.nElDir(iNDim)]; 
        val = bsxfun(@times, val, reshape(space.bFDir{iNDim},indResh));
    end
    
    bFGrad(iDim,:,:,:) = reshape(val,nBfEl,nQNEl,nEl);
end

% REGROUP JACINV such that 
% size(...) = [#dim, #dim, #QN per el, #Elements]
JacInv = reshape(geo.JacInvT(geo.map(pQN)),nDim,nDim,nQNEl,nEl);
clear pQN qN qW;

loc2Glob = space.loc2Glob;

errL2 = 0;
for iEl = 1:prod(mesh.nElDir)
    bFunU_h_iEl = reshape(sum(bsxfun(@times,bFun(:,:,iEl),u_h(loc2Glob(:,iEl))),1),nQNEl,1);
    errL2 = errL2 + sum(((bFunU_h_iEl-uEx(:,iEl)).^2).*weights(:,iEl).*abs(detJac(:,iEl)),1);
end
varargout{1} = sqrt(errL2);

if compH1
    errH1s = 0;
    for iEl = 1:prod(mesh.nElDir)
        JIBFgrad_iEl = squeeze(sum(bsxfun(@times,reshape(JacInv(:,:,:,iEl),nDim,nDim,1,nQNEl),reshape(bFGrad(:,:,:,iEl),1,nDim,nBfEl,nQNEl)),2));
        bFGradU_h_iEl = reshape(sum(bsxfun(@times,JIBFgrad_iEl,reshape(u_h(loc2Glob(:,iEl)),1,nBfEl,1)),2),nDim,nQNEl);
        vals_iEl = reshape(sum((bFGradU_h_iEl-uExGrad(:,:,iEl)).^2,1),nQNEl,1);
        errH1s = errH1s + sum(vals_iEl.*weights(:,iEl).*abs(detJac(:,iEl)),1);
    end
    varargout{2} = sqrt(errL2 + errH1s);
end


end

