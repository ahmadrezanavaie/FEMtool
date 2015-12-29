%% Example Heat Equation 2-d in a rectangle
%Solves the heat equation in a rectangle with mixed boundary conditions and
%zero reaction term. The initial solution is taken similar as in the figure
%on Wikipedia (https://en.wikipedia.org/wiki/Heat_equation). We use an
%adaptive time selection algorithm together with the implicit Radau IIA
%time integration method.
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
% clear all
ticTot = tic;

%% Geometry
origin = [ 0; 0];
sizeDom = [ 10; 10];
geo = geo_2d( 'rectangle', origin, sizeDom);
prbl.geo = geo;

%% Simulation Parameters
% Mesh and Space
meth.nElDir = [ 100; 100];
meth.nQRDir = [ 4; 4]; 
meth.rBFDir = [ 1; 1];
meth.lMax = 10000;
% (Initial) Time Step
meth.dt = 1/200;
% Adaptive Time Step
meth.adapt_dt.flag = 1;
meth.adapt_dt.facmax = 2.0;
meth.adapt_dt.facmin = 0.1;
meth.adapt_dt.fac = 0.9;
meth.adapt_dt.tol = 1e-1;
% Time Integration Scheme
meth.timeScheme = 'GL4';
meth.eta = 2/13; % Damping factor for RKC2 method
% Discretize Initial Solution
meth.discIn = 'IP';

%% Problem Parameters
% Time
prbl.Tmax = 20;
% Forcing Term and Diffusion Coefficient and Initial Solution
prbl.f = @(x)( 0*x(1,:) );
prbl.mu = @(x)( ones(1,size(x,2)) );
prbl.sigma = @(x)( ones(1,size(x,2)) );
prbl.u_0 = @(x)( 2*( x(1,:)<=7 & x(1,:)>=3 & x(2,:)<=7 & x(2,:)>=3 &...
    (bsxfun(@minus,x(1,:),3).^2 + bsxfun(@minus,x(2,:),5).^2 >= 1 ) ) );

% Boundary Conditions
prbl.BC.neum_sides = [ 1, 3, 4];
prbl.BC.neum_lim = [ 0, 0, 0; 1, 1, 1];
prbl.BC.neum_fun{1} = @(x)( 0*x(1,:) );
prbl.BC.neum_fun{2} = @(x)( 0*x(1,:) );
prbl.BC.neum_fun{3} = @(x)( 0*x(1,:) );

prbl.BC.dir_sides = 2;
prbl.BC.dir_lim = [ 0; 1];
prbl.BC.dir_fun{1} = @(x)( 0*x(1,:) );


%% Solve Problem
output.screen_flag = 1;
output.vtk_flag = 1;
output.vtk_name = ['VTK/output/rect2_', meth.timeScheme, '_'];
output.vtk_freq = 1;
output.lmA_flag = 0;
output.smA_flag = 0;

fprintf('Completed: ');
[u_h, data] = solve_Heat2d( prbl, meth, output);
fprintf('\n');

fprintf('\n              h/H:   1/%g\n',min(meth.nElDir))
fprintf(' Time int. scheme:   %s\n',meth.timeScheme)
fprintf('       Time steps:   %i (%i)\n',length(find(data.dtHist(2,:))),size(data.dtHist,2))
fprintf('     Time elapsed:   %g [s]\n\n',toc(ticTot))


