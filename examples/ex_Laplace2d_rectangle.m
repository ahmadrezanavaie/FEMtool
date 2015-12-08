%% Example Laplace 2-d in a rectangle
%Solves the Poisson problem in a rectangle with mixed boundary conditions.
%The numerical solution is compared to a given exact solution.
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

% close all
% clear all

%% Geometry
origin = [ 10; 20];
sizeDom = [ 2; 5];
geo = geo_2d( 'rectangle', origin, sizeDom);
prbl.geo = geo;

%% Simulation Parameters
meth.nElDir = [ 40; 40];     % Mesh Size
meth.nQRDir = [ 2; 2];       % Order of Quadrature Rule
meth.rBFDir = [ 1; 1];       % Order of FEM Basis Functions

%% Exact Solution
origin = geo.map(geo.paramDom.origin);
sizeDom = geo.map(geo.paramDom.origin + geo.paramDom.sizeDir)-origin;
a = 0.7/sizeDom(1); b = 2*pi/sizeDom(2);
uEx = @(x)( exp( -a*(x(1,:)-origin(1)) ).*sin( b*(x(2,:)-origin(2)) ) );
uExGrad = @(x)( [-a*uEx(x); b*exp( -a*(x(1,:)-origin(1)) ).*cos( b*(x(2,:)-origin(2)) ) ] );

%% Problem Parameters
% Forcing Term and Diffusion Coefficient
prbl.f = @(x)( (b^2-a^2)*uEx(x) );
prbl.mu = @(x)( ones(1,size(x,2)) );

% Boundary Conditions
prbl.BC.neum_sides = [ 1, 3, 4];
prbl.BC.neum_lim = [ 0, 0, 0; 0.45, 0.5, 1];
prbl.BC.neum_fun{1} = @(x)(  a*sin( b*(x(2,:)-origin(2)) ) );
prbl.BC.neum_fun{2} = @(x)( -b*exp(-a*(x(1,:)-origin(1))) );
prbl.BC.neum_fun{3} = @(x)(  b*exp(-a*(x(1,:)-origin(1))) );

prbl.BC.dir_sides = [ 1, 2, 3];
prbl.BC.dir_lim = [ 0.45, 0, 0.5; 1, 1, 1];
prbl.BC.dir_fun{1} = @(x)( uEx(x) );
prbl.BC.dir_fun{2} = @(x)( uEx(x) );
prbl.BC.dir_fun{3} = @(x)( uEx(x) );

%% Solve Problem
[u_h, mesh, space] = solve_Laplace2d( prbl, meth);

%% Compare to Exact Solution
[errL2, errH1] = spErrSol( space, mesh, geo, u_h, uEx, uExGrad );
fprintf('\n   ||u_h - u||_{L2} = %g\n   ||u_h - u||_{H1} = %g\n\n',errL2,errH1)

%% Plot Solution
% % % paramO = geo.paramDom.origin;
% % % paramS = geo.paramDom.sizeDir;
% % % nDofDir = space.nDofDir;
% % % [X, Y] = meshgrid(linspace(paramO(1),paramO(1)+paramS(1),nDofDir(1)),...
% % %     linspace(paramO(2),paramO(2)+paramS(2),nDofDir(2)));
% % % sizeResh = size(X');
% % % Z = geo.map([reshape(X',1,prod(sizeResh));...
% % %     reshape(Y',1,prod(sizeResh))]);
figure
MX = reshape(Z(1,:),sizeResh);
MY = reshape(Z(2,:),sizeResh);
surf(MX(1:meth.rBFDir(1):end,1:meth.rBFDir(1):end),MY(1:meth.rBFDir(2):end,1:meth.rBFDir(2):end),zeros(meth.nElDir'+1))
title('Mesh')
view([0,90])
figure
surf(reshape(Z(1,:),sizeResh),reshape(Z(2,:),sizeResh),reshape(u_h,nDofDir(1),nDofDir(2)),'EdgeColor','none')
colorbar
title('Numerical Solution u_h')
