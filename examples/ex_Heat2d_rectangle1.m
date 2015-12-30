%% Example Heat Equation 2-d in a rectangle
%Solves the heat equation in a rectangle with mixed boundary conditions.
%The numerical solution is compared to a given exact solution from which
%the boundary- and the initial conditions as well as the reaction term are
%derived. In this example we use a fixed time step and the implicit
%Crank-Nicolson time integration scheme.
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
ticTot = tic;

%% Geometry
origin = [ 0; 0];
sizeDom = [ 10; 2];
geo = geo_2d( 'rectangle', origin, sizeDom);
prbl.geo = geo;

%% Simulation Parameters
% Mesh and Space
meth.nElDir = [ 100; 20];
meth.nQRDir = [ 4; 4]; 
meth.rBFDir = [ 1; 1];
meth.lMax = 10000;
% (Initial) Time Step
meth.dt = 1/25;
% Adaptive Time Step
meth.adapt_dt.flag = 0;
meth.adapt_dt.facmax = 2.0;
meth.adapt_dt.facmin = 0.1;
meth.adapt_dt.fac = 0.9;
meth.adapt_dt.tol = 1e-1;
% Time Integration Scheme
meth.timeScheme = 'RIIA';
% meth.eta = 2/13; % Damping factor for RKC2 method
% Discretize Initial Solution
meth.discIn = 'IP';

%% Exact Solution
origin = geo.map(geo.paramDom.origin);
sizeDom = geo.map(geo.paramDom.origin + geo.paramDom.sizeDir)-origin;
a = 0.7/sizeDom(1); b = 2*pi/sizeDom(2); c = 0.1;
uEx = @(x,t)( exp(-c*t)*exp( -a*(x(1,:)-origin(1)) ).*sin( b*(x(2,:)-origin(2)) ) );
uExGrad = @(x,t)( [-a*uEx(x,t); b*exp(-c*t)*exp( -a*(x(1,:)-origin(1)) ).*cos( b*(x(2,:)-origin(2)) ) ] );

%% Problem Parameters
% Time
prbl.Tmax = 5;
% Forcing Term and Diffusion Coefficient and Initial Solution
prbl.f = @(x,t)( (-c+b^2-a^2)*uEx(x,t) );
prbl.mu = @(x)( ones(1,size(x,2)) );
prbl.sigma = @(x)( ones(1,size(x,2)) );
prbl.u_0 = @(x)( uEx(x,0) );

% Boundary Conditions
prbl.BC.neum_sides = [ ];

prbl.BC.dir_sides = [ 1, 2, 3, 4];
prbl.BC.dir_lim = [ 0, 0, 0, 0; 1, 1, 1, 1];
prbl.BC.dir_fun{1} = @(x,t)( uEx(x,t) );
prbl.BC.dir_fun{2} = @(x,t)( uEx(x,t) );
prbl.BC.dir_fun{3} = @(x,t)( uEx(x,t) );
prbl.BC.dir_fun{4} = @(x,t)( uEx(x,t) );

% prbl.BC.neum_sides = [ 1, 2, 3, 4 ];
% prbl.BC.neum_lim = [ 0, 0, 0, 0; 1, 1, 1, 1];
% prbl.BC.neum_fun{1} = @(x,t)(  a*uEx(x,t) );
% prbl.BC.neum_fun{2} = @(x,t)( -a*uEx(x,t) );
% prbl.BC.neum_fun{3} = @(x,t)( -b*exp(-c*t)*exp( -a*(x(1,:)-origin(1)) ).*cos( b*(x(2,:)-origin(2)) ) );
% prbl.BC.neum_fun{4} = @(x,t)(  b*exp(-c*t)*exp( -a*(x(1,:)-origin(1)) ).*cos( b*(x(2,:)-origin(2)) ) );
% 
% prbl.BC.dir_sides = [];


%% Solve Problem
output.screen_flag = 1;
output.vtk_flag = 0;
output.vtk_name = ['VTK/output/rect1_', meth.timeScheme, '_'];
output.vtk_freq = 1;
output.lmA_flag = 0;
output.smA_flag = 1;

fprintf('Completed: ');
[u_h, data, mesh, space] = solve_Heat2d( prbl, meth, output);
fprintf('\n')
dt = data.dtHist(1,data.dtHist(2,:)>0);
tn = cumsum([0, dt]);

%% Write Reference Solution VTK file
% if output.vtk_flag
%     refName = 'VTK/output/ref_';
%     sigma = @(x)( 0*x(1,:) + 1 );
%     M = opUV(  space, space, mesh, sigma, geo  );
%     for i = 1:output.vtk_freq:size(u_h,2)
%         fun = @(x)( uEx(x,tn(i)) );
%         matlab2vtk( M\opFV(space,mesh,fun,geo), space, geo, i-1, refName)
%     end
% end


%% Compare to Exact Solution
for it = 1:length(dt)
    [errL2_i(it), errH1_i(it)] = spErrSol( space, mesh, geo, u_h(:,it+1), @(x)(uEx(x,tn(it+1))), @(x)(uExGrad(x,tn(it+1))) );
end
errT_totL2 = errL2_i(end).^2 + 2*data.smA*sum(dt.*(errH1_i.^2));
fprintf('\n   ||u_h(T)-u(T)||_{L2}^2 + 2*alpha*sum(dt*||u_h(tn)-u(tn)||_{V}^2) = %g\n',errT_totL2)
errT_totH1 = errH1_i(end).^2 + 2*data.smA*sum(dt.*(errH1_i.^2));
fprintf('   ||u_h(T)-u(T)||_{V}^2  + 2*alpha*sum(dt*||u_h(tn)-u(tn)||_{V}^2) = %g\n',errT_totH1)
errT_L2 = errL2_i(end);
fprintf('   ||u_h(T)-u(T)||_{L^2}                                            = %g\n',errT_L2)
errT_H1 = errH1_i(end);
fprintf('   ||u_h(T)-u(T)||_{H^1}                                            = %g\n\n',errT_H1)
fprintf(' Time elapsed:   %g [s]\n\n',toc(ticTot))


