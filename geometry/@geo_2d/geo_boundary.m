function [ geo_bndry ] = geo_boundary( geo, iSide )
%GEO_BOUNDARY is a method of the class geo_2d aiming to define a 1D
%boundary domain \Gamma of the computational domain \Omega defined by
%geo_2d. If F: P -> \Omega is the bijective parametrization map from the
%parameter domain P = [a,b] x [c,d] into the computational domain \Omega = 
%F(P), we identify the 4 different boundary parts as
%
%   \Gamma_1 = F( a, [c,d] )
%   \Gamma_2 = F( b, [c,d] )
%   \Gamma_3 = F( [a,b], c )
%   \Gamma_4 = F( [a,b], d ).
%
%We define a parameterization function gTilde(x,iSide) such that:
%
%   gTilde(x,1) = [ a; x ]
%   gTilde(x,2) = [ b; x ]
%   gTilde(x,3) = [ x; c ]
%   gTilde(x,4) = [ x; d ].
%
%
%However, it is not necessary to define the same type of boundary condition
%on the entire boundary \Gamma_i (it can be splitted into parts -> see
%mesh_2d/mesh_boundary).
%
%
% -------------------------------------------------------------------------
%
%   geo_bndry = geo_boundary( geo, iSide );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% geo           (see output of geo_2d.m)
% iSide         (1 x 1 integer) Indicating the identification number of the
%               boundary
%
%
% OUTPUT
% ------
% geo_bndry     (1 x 1 struct)
%   .map        (function handle) Definition of the parameterization
%               function F : I -> \Gamma_i, where I = [a,b] or I = [c,d].
%   .detJac     (function handle) Definition of the derivative of F
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

paramOrigin = geo.paramDom.origin;
paramSizeDir = geo.paramDom.sizeDir;

geo_bndry.map = @(x)( geo.map(gTilde(x,iSide,paramOrigin,paramSizeDir)) );
geo_bndry.detJac = @(x)( detJac( geo.Jac(gTilde(x,iSide,paramOrigin,paramSizeDir)), iSide ) );

geo_bndry.type = '2-D Boundary Domain';
end

function der = detJac( evalJac, iSide)

if iSide == 1 || iSide == 2
    der = sqrt(squeeze(sum(evalJac(:,2,:).^2,1)))';
elseif iSide == 3 || iSide == 4
    der = sqrt(squeeze(sum(evalJac(:,1,:).^2,1)))';
else
    error('geo_2d/geo_boundary: Wrong side number index')
end

end

function xExp = gTilde(x,iSide,paramOrigin,paramSizeDir)
if iSide == 1
    xExp = [0*x + paramOrigin(1); x];
elseif iSide == 2
    xExp = [0*x + paramOrigin(1) + paramSizeDir(1); x];
elseif iSide == 3
    xExp = [x; 0*x + paramOrigin(2)];
elseif iSide == 4
    xExp = [x; 0*x + paramOrigin(2) + paramSizeDir(2)];
else
    error('geo_2d/geo_boundary: Wrong side number index')
end

end

