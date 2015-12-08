function [ centr ] = mesh_elCentrCoord( mesh, geo )
%MESH_ELCENTRCOORD Is a method of the class mesh_2d and computes the center
%coordinates of the quadrilateral mesh elements in the PHYSICAL DOMAIN.
%
%
% -------------------------------------------------------------------------
%
%   centr = mesh_elCentrCoord( mesh, geo );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% mesh          (see output mesh_2d.m)
% geo           (see output geo_2d.m)
%
%
% OUTPUT
% ------
% centr         (2 x mesh.nElTot double) Physical center coordinates
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

xx = origin(1) + (((1:mesh.nElDir(1))-0.5)*mesh.dElDir(1));
yy = origin(2) + (((1:mesh.nElDir(2))-0.5)*mesh.dElDir(2));

xx = reshape(repmat(xx',1,mesh.nElDir(2)),1,mesh.nElTot);
yy = reshape(repmat(yy,mesh.nElDir(1),1),1,mesh.nElTot);

centr = geo.map([xx;yy]);

end

