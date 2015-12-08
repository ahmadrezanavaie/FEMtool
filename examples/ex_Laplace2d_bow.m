%% Example Laplace 2-d in a bow
%Solves the Poisson problem in a bow-shaped computational domain with
%Dirichlet boundary conditions. The numerical solution is compared to a
%given exact solution.  
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

close all
clear all

%% Geometry
geo = geo_2d('bow');
prbl.geo = geo;

%% Simulation Parameters
meth.nElDir = [ 40; 40];     % Mesh Size
meth.nQRDir = [ 3; 3];       % Order of Quadrature Rule
meth.rBFDir = [ 1; 2];       % Order of FEM Basis Functions

%% Exact Solution
uEx = @(x)( sin( x(1,:)*pi ).*cos( x(2,:)*pi ) );
uExGrad = @(x)( [ pi*cos( x(1,:)*pi ).*cos( x(2,:)*pi );...
    -pi*sin( x(1,:)*pi ).*sin( x(2,:)*pi ) ] );

%% Problem Parameters
% Forcing Term and Diffusion Coefficient
prbl.f = @(x)( (2*pi.^2)*uEx(x) );
prbl.mu = @(x)( ones(1,size(x,2)) );

% Boundary Conditions
prbl.BC.neum_sides = [];

prbl.BC.dir_sides = [ 1, 2, 3, 4];
prbl.BC.dir_lim =[ 0, 0, 0, 0; 1, 1, 1, 1];
prbl.BC.dir_fun{1} = @(x)( uEx(x) );
prbl.BC.dir_fun{2} = @(x)( uEx(x) );
prbl.BC.dir_fun{3} = @(x)( uEx(x) );
prbl.BC.dir_fun{4} = @(x)( uEx(x) );

%% Solve Problem
[sol, mesh, space] = solve_Laplace2d( prbl, meth);

%% Compare to Exact Solution
[errL2, errH1] = spErrSol( space, mesh, geo, sol, uEx, uExGrad );
fprintf('\n   ||u_h - u||_{L2} = %g\n   ||u_h - u||_{H1} = %g\n\n',errL2,errH1)

%% Plot Solution
paramO = geo.paramDom.origin;
paramS = geo.paramDom.sizeDir;
nDofDir = space.nDofDir;
[X, Y] = meshgrid(linspace(paramO(1),paramO(1)+paramS(1),nDofDir(1)),...
    linspace(paramO(2),paramO(2)+paramS(2),nDofDir(2)));
sizeResh = size(X');
Z = geo.map([reshape(X',1,prod(sizeResh));...
    reshape(Y',1,prod(sizeResh))]);
figure
MX = reshape(Z(1,:),sizeResh);
MY = reshape(Z(2,:),sizeResh);
surf(MX(1:meth.rBFDir(1):end,1:meth.rBFDir(2):end),MY(1:meth.rBFDir(1):end,1:meth.rBFDir(2):end),zeros(meth.nElDir'+1))
title('Mesh')
view([0,90])
figure
surf(reshape(Z(1,:),sizeResh),reshape(Z(2,:),sizeResh),reshape(sol,nDofDir(1),nDofDir(2)),'EdgeColor','none')
colorbar
title('Numerical Solution u_h')

