classdef mesh_2d
%MESH_2D generates a mesh of uniform quadrilateral elements of size dx x dy
%in a rectangular PARAMETRIC domain (see geo_2d).
%
%
%                        ( origin(1)+size(1), origin(2)+size(2) )
%                   ___________________
%                  |   |   |   |   |   |      Local Mesh Vertices:
%                  |   |   |   |...| n |           
%                  |___|___|___|___|___|          3 ______ 4
%       p2         |   |   |   |   |   |           |      |
%        ^         |   |   |   |   |   |           | el_i |
%        |         |___|___|___|___|___|           |______|
%       -|---> p1  |   |   |   |   |   |          1        2
%                  |   |   |   |   |   |
%                  |___|___|___|___|___|
%                  |   |   |   |   |   |
%                  | 1 | 2 | 3 |...|   | dy = size(2)/nElY
%                  |___|___|___|___|___|
%                                   dx = size(1)/nElX
%     ( origin(1), origin(2) )
%
%
% -------------------------------------------------------------------------
%
%   mesh = mesh_2d( geo, nElDir, nQRDir);
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% geo           (see output of geo_2d.m --> or define your own)
% nElDir        (2 x 1 double) Number of mesh elements in each parametric
%               direction
% nQRDir        (2 x 1 double) Order of tensor product quadrature rule (QR)
%               in each parametric direction
%
%
% OUTPUT
% ------
% mesh              (1 x 1 struct)
%   .paramDom       (a copy of geo.paramDom)
%   .nElDir         (2 x 1 double) Number of elements in each direction
%   .nElTot         (1 x 1 double) Total number of elements
%   .dElDir         (2 x 1 double) Element size in each direction
%   .nVertDir       (2 x 1 double) Number of vertices in each direction
%   .vertDir        (2 x 1 cell-array) 
%       .vertDir{i} (1 x nVertDir(i) double) Contains the unique
%                   i-coordinates of the mesh vertices
%   .nQNElDir       (2 x 1 double) Number of quadrature nodes per element
%                   in each parametric direction
%   .qNDir          (1 x 2 cell-array)
%       .qNDir{i}   (nQNElDir(i) x nElDir(i) double) Contains the quadr.
%                   nodes in each element in parametric direction i.
%   .qWDir          (1 x 2 cell-array) 
%       .qWDir{i}   (nQNElDir(i) x nElDir(i) double) Contains the quadr.
%                   weights in each element in parametric direction i.
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


    
    properties
        type = '2-d Mesh';
        paramDom = [];
        nElDir = [];
        nElTot = [];
        dElDir = [];
        nVertTot = [];
        nVertDir = [];
        vertDir = [];
        nQNElDir = [];
        qNDir = [];
        qWDir = [];
    end
    
    methods
        function obj = mesh_2d( geo, nElDir, nQRDir)
            % Parametric Domain
            origin = geo.paramDom.origin;
            sizeDir = geo.paramDom.sizeDir;
            obj.paramDom = geo.paramDom;

            % Elements
            obj.nElDir = nElDir;
            nEl = nElDir(1)*nElDir(2);
            obj.nElTot = nEl;
            dEl = sizeDir./nElDir;
            obj.dElDir = dEl;

            % Vertices
            obj.nVertTot = (nElDir(1)+1)*(nElDir(2)+1);
            obj.nVertDir = [(nElDir(1)+1);(nElDir(2)+1)];
            obj.vertDir = {origin(1):dEl(1):(origin(1)+sizeDir(1));...
                origin(2):dEl(2):(origin(2)+sizeDir(2))};

            % Quadrature Rule
            obj.nQNElDir = nQRDir;
            [p1QR, w1QR] = lgwt(nQRDir(1),-1,1); p1QR = flipud(p1QR); w1QR = flipud(w1QR);
            [p2QR, w2QR] = lgwt(nQRDir(2),-1,1); p2QR = flipud(p2QR); w2QR = flipud(w2QR);

            qN = cell(2,1);
            qW = cell(2,1);

            cen1 = origin(1) + ((1:nElDir(1))-1/2)*dEl(1);
            cen2 = origin(2) + ((1:nElDir(2))-1/2)*dEl(2);
            qN{1} = bsxfun(@plus,p1QR*dEl(1)/2,cen1);
            qN{2} = bsxfun(@plus,p2QR*dEl(2)/2,cen2);
            qW{1} = repmat(w1QR*(dEl(1)/2),1,length(cen1));
            qW{2} = repmat(w2QR*(dEl(2)/2),1,length(cen2));
            obj.qNDir = qN;
            obj.qWDir = qW;

        end
    end
    
end

