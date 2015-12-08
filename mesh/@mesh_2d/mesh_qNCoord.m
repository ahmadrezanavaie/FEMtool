function [ qNCoord ] = mesh_qNCoord( mesh, geo )
%MESH_QNCOORD Is a method of the class mesh_2d and computes the coordinates
%of the quadrature nodes in the physical domain.
%
%
% -------------------------------------------------------------------------
%
%   qNCoord = mesh_qNCoord( mesh, geo );
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
% qNCoord       (2 x prod(mesh.nQNElDir) x mesh.nElTot double) Physical
%               coordinates of the QN
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

xx = reshape(mesh.qNDir{1},mesh.nQNElDir(1),1,mesh.nElDir(1),1);
yy = reshape(mesh.qNDir{2},1,mesh.nQNElDir(2),1,mesh.nElDir(2));

xx = reshape(repmat(xx,1,mesh.nQNElDir(2),1,mesh.nElDir(2)),1,prod(mesh.nQNElDir)*prod(mesh.nElDir));
yy = reshape(repmat(yy,mesh.nQNElDir(1),1,mesh.nElDir(1),1),1,prod(mesh.nQNElDir)*prod(mesh.nElDir));

qNCoord = reshape(geo.map([xx;yy]),2,prod(mesh.nQNElDir),prod(mesh.nElDir));


end