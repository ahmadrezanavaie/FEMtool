function [ ] = test_Heat2d_stream( )
%TEST_HEAT2D_STREAM tests the numerical solution against its analytical one
%given in the book
%
%"Quantitative Hydrogeology" by Ghislain de Marsily (1986)
%
%in section 8.5. We solve the problem in two dimensions but independant of
%the second one.
%
%The problem consist of an aquifer initially in equilibrium next to a
%stream. At time t = 0 the aquifer receives a uniform and instantaneous
%recharge h_0 all over its surface. The exact solution is given on page
%198-201.
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

%% Geometry
L = 100;
origin = [ 0; 0];
sizeDom = [ L; 10];
geo = geo_2d( 'rectangle', origin, sizeDom);
prbl.geo = geo;

%% Simulation Parameters
% Mesh and Space
meth.nElDir = [ 40; 10];
meth.nQRDir = [ 4; 4]; 
meth.rBFDir = [ 1; 1];
meth.lMax = 1000;
% (Initial) Time Step
meth.dt = 1/10;
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

%% Problem Parameters
T = 3.5*1e-4;
S = 1.2*1e-5;
h_0 = 0.1;
% Time
prbl.Tmax = 20;
% Reaction Term and Diffusion Coefficient and Initial Solution
prbl.f = @(x)( 0*x(1,:) );
muFact = T/S;
prbl.mu = @(x)( muFact*ones(1,size(x,2)) );
prbl.sigma = @(x)( ones(1,size(x,2)) );
prbl.u_0 = @(x)( 0*x(1,:) + h_0 );

% Boundary Conditions
prbl.BC.neum_sides = [ 2, 3, 4];
prbl.BC.neum_lim = [ 0, 0, 0; 1, 1, 1];
prbl.BC.neum_fun{1} = @(x)( 0*x(1,:) );
prbl.BC.neum_fun{2} = @(x)( 0*x(1,:) );
prbl.BC.neum_fun{3} = @(x)( 0*x(1,:) );

prbl.BC.dir_sides = 1;
prbl.BC.dir_lim = [ 0; 1];
prbl.BC.dir_fun{1} = @(x)( 0*x(1,:) );


%% Solve Problem
output.screen_flag = 1;
output.vtk_flag = 0;
output.vtk_name = ['VTK/output/stream_', meth.timeScheme, '_'];
output.vtk_freq = 1;
output.lmA_flag = 0;
output.smA_flag = 0;

fprintf('Completed: ')
[u_h, data, ~, space] = solve_Heat2d( prbl, meth, output);
fprintf('\n')
dt = data.dtHist(1,data.dtHist(2,:)>0);
tn = cumsum([0, dt]);

xDof = space_dofCoord(space,geo);
xDof = xDof(1,1:space.nDofDir(1));

timeStps = [51, 101, 151, 201];
figure
for iTme = 1:4
    uEx = h_0*(1-erfc(xDof*sqrt(S/(4*T*tn(timeStps(iTme))))));
    subplot(2,2,iTme)
    plot(xDof,uEx,'k','LineWidth',2)
    hold on
    plot(xDof,u_h(1:space.nDofDir(1),timeStps(iTme)),'r--','LineWidth',2)
    hold off
    title(sprintf('Time tn = %g [s]',round(10*tn(timeStps(iTme)))/10))
    legend('Ex. Sol.','u_h(x)','Location','SouthEast')
    xlabel('\bf{x}')
end

end

