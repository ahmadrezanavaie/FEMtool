classdef space_2d
%SPACE2D generates the FE function space on the parametric domain, where
%the mesh elements are quadrilaterals. Depending on the order of the
%piecewise polynomial basis functions we may add some degrees of freedom to
%the vertices in the mesh in order to guarantee continuity accross the
%element boundaries. The local degrees of freedom are enumerated as
%follows:
%
%         order = 1                   order = 2
%         ---------                   ---------
%        ___________                 _____ _____
%       3           4               7     8     9
%       |           |               |           |
%       |           |     K_ref     4     5     6
%       |           |               |           |
%       1___________2               1_____2_____3
%
%
%         order = 3                   order = 4
%         ---------                   ---------
%        ___ ___ ___                 __ __ __ __
%      13  14  15   16             21 22 23 24  25
%       9  10  11   12             16 17 18 19  20
%       |           |     K_ref    11 12 13 14  15
%       5   6   7   8               6  7  8  9  10
%       1___2___3___4               1__2__3__4__5
%  p2      
%   ^        
%   |     
%  -|---> p1
%
% 
%   !!!!!!!!!!!!   CAUTION   !!!!!!!!!!!! 
%The piecewise polynomial basis functions of order r in a quadrilateral
%mesh element are of the form
%
%   f(x,y) = p(x)q(y),
%
%with deg(p),deg(q) <= r. Therefore the gradient of a basis function reads
%
%                   / p'(x)q(y) \
%   grad(f(x,y)) = |             | .
%                   \ p(x)q'(y) /
%
%In the class structure "space = space_2d(...)" we only store the
%polynomials and its derivatives evaluated at the quadrature nodes in the
%corresponding direction. The gradient evaluated at the tensor product
%quadrature points in an element el given by
%
%   el = el1 + (el2-1)*mesh.nElDir(1)
%
%can be obtained by
%
%   [X, Y] = meshgrid(space.bFDerDir{1}(i,:,el1),space.bFDir{2}(j,:,el2)))
%   df_ij/dx = X'.*Y';
%   [X, Y] = meshgrid(space.bFDir{1}(i,:,el1),space.bFDerDir{2}(j,:,el2)))
%   df_ij/dy = X'.*Y';
%
%
% -------------------------------------------------------------------------
%
%   space = space_2d( mesh, rBFDir );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% mesh              (see output of mesh_2d.m)
% rBFDir            (2 x 1 double) Order of the piecewise polynomial basis
%                   functions in each parametric direction
%
%
% OUTPUT
% ------
% space             (1 x 1 struct)
%   .nDim           (1 x 1 double) Dimension of the computational domain
%   .nDof           (1 x 1 double) Number of degrees of freedom
%   .nDofDir        (2 x 1 double) Number of degrees of freedom in each
%                   parametric direction
%   .rBFDir         (2 x 1 double) Order of the piecewise polynomial basis
%                   functions in each parametric direction
%   .nDofElDir      (2 x 1 double) Number of degrees of freedom per element
%                   in each parametric direction
%   .bFDir          (2 x 1 cell-array)
%       .bFDir{i}   (nDofElDir(i) x mesh.nQNElDir(i) x mesh.nElDir(i)
%                   double ) For each element in the different directions,
%                   this contains all the non-zero basis functions
%                   evaluated at all the quadrature nodes
%   .bFDerDir       (2 x 1 cell-array)
%       .bFDerDir{i}(nDofElDir(i) x mesh.nQNElDir(i) x mesh.nElDir(i)
%                   double ) For each element in the different directions,
%                   this contains all the derivatives of the non-zero basis
%                   functions evaluated at all the quadrature nodes
%   .loc2Glob       (prod(nDofElDir) x prod(mesh.nElDir) double) For each
%                   element this local-to-global map indicates the global
%                   dof indices for each non-zero basis function in a given
%                   element
%
%
% Remark:
% -------
% We construct basis functions f_ij(x,y) = p_i(x)q_j(y) with deg(p_i) = r
% and deg(q_j) = s such that
%
%   p_i(x) = (x-r_i1)(x-r_i2)*...*(x-r_ir)/c_i1
%   q_j(y) = (y-s_j1)(y-s_j2)*...*(y-s_js)/c_j2
%
% and c_i1 resp. c_j2 such that p_i(x_i) = 1 = q_j(y_j).
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
    
    properties
        type = 'Function Space on 2-d Domain';
        nDim = [];
        nDof = [];
        nDofDir = [];
        rBFDir = [];
        nDofElDir = [];
        bFDir = [];
        bFDerDir = [];
        loc2Glob = [];
    end
    
    methods
        function obj = space_2d( mesh, rBFDir )
            if rBFDir(1) == 0 || rBFDir(2) == 0
                warning(['Not sure that it works for constant basis functions ',...
                    'because of the definition of pts = linspace(...), ',...
                    'that places the point at the end and not in the center --> check this']);
            end

            % Space Dimensions
            obj.nDim = 2;

            % Polynomial Order and Degrees of Freedom
            obj.rBFDir = rBFDir;
            obj.nDofElDir = rBFDir+1;
            obj.nDof = prod(rBFDir.*mesh.nElDir + 1,1);
            obj.nDofDir = rBFDir.*mesh.nElDir + 1;

            % Basis Functions and Derivatives Evaluated at Quadrature Nodes
            bF = {zeros(obj.nDofElDir(1),mesh.nQNElDir(1),mesh.nElDir(1));...
                zeros(obj.nDofElDir(2),mesh.nQNElDir(2),mesh.nElDir(2))};
            bFDer = {zeros(obj.nDofElDir(1),mesh.nQNElDir(1),mesh.nElDir(1));...
                zeros(obj.nDofElDir(2),mesh.nQNElDir(2),mesh.nElDir(2))};

            indRoots1 = repmat((1:rBFDir(1)+1)',1,rBFDir(1)+1); indRoots1 = indRoots1(~eye(size(indRoots1)));
            indRoots1 = reshape(indRoots1,rBFDir(1),rBFDir(1)+1)';
            indRootsDer1 = repmat((1:rBFDir(1))',1,rBFDir(1)); indRootsDer1 = indRootsDer1(~eye(size(indRootsDer1)))';
            for iEl1 = 1:mesh.nElDir(1)
                vert = mesh.vertDir{1}(iEl1:iEl1+1);
                pts = linspace(vert(1),vert(2),rBFDir(1)+1)';

                rootsBF1 = pts(indRoots1);
                constBF1 = prod(bsxfun(@minus,pts,rootsBF1),2);

                evalFactors = bsxfun(@minus,mesh.qNDir{1}(:,iEl1)',reshape(rootsBF1,size(rootsBF1,1),1,size(rootsBF1,2)));
                bF{1}(:,:,iEl1) = bsxfun(@rdivide,prod(evalFactors,3),constBF1);

                rootsBFDer1 = reshape(rootsBF1(:,indRootsDer1),rBFDir(1)+1,1,rBFDir(1)-1,rBFDir(1));
                evalFactorsDer = bsxfun(@minus,mesh.qNDir{1}(:,iEl1)',rootsBFDer1);
                bFDer{1}(:,:,iEl1) = bsxfun(@rdivide,sum(prod(evalFactorsDer,3),4),constBF1);
            end
            indRoots2 = repmat((1:rBFDir(2)+1)',1,rBFDir(2)+1); indRoots2 = indRoots2(~eye(size(indRoots2)));
            indRoots2 = reshape(indRoots2,rBFDir(2),rBFDir(2)+1)';
            indRootsDer2 = repmat((1:rBFDir(2))',1,rBFDir(2)); indRootsDer2 = indRootsDer2(~eye(size(indRootsDer2)))';
            for iEl2 = 1:mesh.nElDir(2)
                vert = mesh.vertDir{2}(iEl2:iEl2+1);
                pts = linspace(vert(1),vert(2),rBFDir(2)+1)';

                rootsBF2 = pts(indRoots2);
                constBF2 = prod(bsxfun(@minus,pts,rootsBF2),2);

                evalFactors = bsxfun(@minus,mesh.qNDir{2}(:,iEl2)',reshape(rootsBF2,size(rootsBF2,1),1,size(rootsBF2,2)));
                bF{2}(:,:,iEl2) = bsxfun(@rdivide,prod(evalFactors,3),constBF2);

                rootsBFDer2 = reshape(rootsBF2(:,indRootsDer2),rBFDir(2)+1,1,rBFDir(2)-1,rBFDir(2));
                evalFactorsDer = bsxfun(@minus,mesh.qNDir{2}(:,iEl2)',rootsBFDer2);
                bFDer{2}(:,:,iEl2) = bsxfun(@rdivide,sum(prod(evalFactorsDer,3),4),constBF2);
            end
            obj.bFDir = bF;
            obj.bFDerDir = bFDer;

            indBf1 = bsxfun(@plus,(1:obj.nDofElDir(1))',0:(obj.nDofElDir(1)-1):((obj.nDofElDir(1)-1)*(mesh.nElDir(1)-1)));
            indBf1 = reshape(indBf1,obj.nDofElDir(1),1,mesh.nElDir(1),1);
            indBf2 = bsxfun(@plus,(1:obj.nDofElDir(2))',0:(obj.nDofElDir(2)-1):((obj.nDofElDir(2)-1)*(mesh.nElDir(2)-1)));
            indBf2 = reshape(indBf2,1,obj.nDofElDir(2),1,mesh.nElDir(2));
            obj.loc2Glob = reshape(bsxfun(@plus, indBf1, (indBf2-1)*(mesh.nElDir(1)*(obj.nDofElDir(1)-1)+1) ),...
                prod(obj.nDofElDir),prod(mesh.nElDir));
        end
    end
    
end

