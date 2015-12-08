function [ evalU ] = spEval( space, mesh, geo, u_h, opt )
%SPEVAL computes the value or the derivatives, of a function given by its
%degrees of freedom, at the quadrature nodes defined in the mesh. We recall
%that in space_2d and space_3d the basis functions are defined as:
%
% 2d:
%   u(x) = p1(x1)*p2(x2)
%
% 3d:
%   u(x) = p1(x1)*p2(x2)*p3(x3),
%
% therefore the gradients are
%
% 2d:
%                 / p1'(x1)*p2(x2) \
%   grad(u(x)) = |                  |
%                 \ p1(x1)*p2'(x2) /
%
% 3d:
%                 / p1'(x1)*p2(x2)*p3(x3) \
%   grad(u(x)) = |  p1(x1)*p2'(x2)*p3(x3)  |
%                 \ p1(x1)*p2(x2)*p3'(x3) /.
%
%
% -------------------------------------------------------------------------
%
%   evalU = spGradSol( space, mesh, geo, u_h, option );
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
% opt               (string) Possible options are 'value' and 'gradient'
%
%
% OUTPUT
% ------
% evalU
%   'value'         (nBfEl x nQNEl x nEl double) Evaluated values of the
%                   function given by the degrees of freedom
%   'gradient'      (nDim x nBfEl x nQNEl x nEl double) Evaluatede
%                   gradient of the function given by the degrees of
%                   freedom
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

nDim = space.nDim;
nBfEl = prod(space.nDofElDir);
nQNEl = prod(mesh.nQNElDir);
nEl = mesh.nElTot;

% Compute the values
if strcmp(opt,'value')
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
    
    u_h = reshape(u_h(space.loc2Glob),nBfEl,1,nEl);
    evalU = bsxfun(@times,bFun,u_h);
    
% Compute the gradient
elseif strcmp(opt,'gradient')
    
    pQN = zeros(nDim,nQNEl*nEl);
    for iDim = 1:nDim
        indResh = ones(1,2*nDim);
        indResh([iDim, nDim + iDim]) = [mesh.nQNElDir(iDim), mesh.nElDir(iDim)];
        indRep = [mesh.nQNElDir', mesh.nElDir'];
        indRep([iDim, nDim + iDim]) = [1, 1];

        qN = reshape(mesh.qNDir{iDim},indResh);
        qN = repmat(qN,indRep);
        pQN(iDim,:) = reshape(qN,1,numel(qN));
    end
    JacInv = reshape(geo.JacInv(geo.map(pQN)),nDim,nDim,1,nQNEl,nEl);
    
    grad = zeros(nDim, nBfEl, nQNEl, nEl);
    for iDim = 1:nDim
        indReshDer = ones(1,3*nDim);
        indReshDer([iDim, nDim+iDim, 2*nDim+iDim]) = [space.nDofElDir(iDim),mesh.nQNElDir(iDim),mesh.nElDir(iDim)]; 
        val = reshape(space.bFDerDir{iDim},indReshDer);

        for iNDim = setdiff(1:nDim,iDim)
            indResh = ones(1,3*nDim);
            indResh([iNDim, nDim+iNDim, 2*nDim+iNDim]) = [space.nDofElDir(iNDim),mesh.nQNElDir(iNDim),mesh.nElDir(iNDim)]; 
            val = bsxfun(@times, val, reshape(space.bFDir{iNDim},indResh));
        end

        grad(iDim,:,:,:) = reshape(val,nBfEl,nQNEl,nEl);
    end
    grad = squeeze(sum( bsxfun(@times, JacInv, reshape(grad,1,nDim,nBfEl,nQNEl,nEl) ) ,2));

    u_h = reshape(u_h(space.loc2Glob),1,nBfEl,1,nEl);
    evalU = bsxfun(@times,grad,u_h);
else
    error('spEval: Wrong option (possible options are value and gradient')
end


end

