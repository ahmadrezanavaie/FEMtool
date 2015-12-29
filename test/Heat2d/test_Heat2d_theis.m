function [ ] = test_Heat2d_theis( )
%TEST_HEAT2D_THEIS Solves the heat equation in a two dimensional square
%domain with zero initial solution, mixed boundary conditions and a
%negative point-source term f. We use an adaptive time selection algorithm
%together with the implicit Radau IIA time integration method.
%
%The theoretical drawdown in an infinite domain s = s(r,t) = h(r,t) - h_0
%is given by the THEIS EQUATION and reads
%
%   s = Q/(4*pi*T)*E_1(u), where u(r,t) = r^2*S/(4*T*t)
%
%where
%   Q is the pumping rate [m^3/s],
%   T is the transmissivity [m^2/s],
%   S is the storativity [-],
%   r is the distance from the well, and
%   E_1(.) denotes the exponentional integral.
%
%In this function we aim to compare the numerical soluton in a finite
%domain to the theoretical one.
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

%% Geometry
origin = [ 0; 0];
sizeDom = [ 40; 40];
geo = geo_2d( 'rectangle', origin, sizeDom);
prbl.geo = geo;

%% Simulation Parameters
% Mesh and Space
meth.nElDir = [ 200; 200];
meth.nQRDir = [ 4; 4]; 
meth.rBFDir = [ 1; 1];
meth.lMax = 1000;
% (Initial) Time Step
meth.dt = 1/10;
% Adaptive Time Step
meth.adapt_dt.flag = 1;
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
litrPerMin = 100;
Q = litrPerMin/60000;

xWell = origin + 0.5*sizeDom;

volBasisFun = 1*prod(sizeDom)/prod(meth.nElDir);
reaction = Q/(volBasisFun*S);

% Time
prbl.Tmax = 1000;
% Reaction Term and Diffusion Coefficient and Initial Solution
prbl.f = @(x)( -reaction*( (abs(x(1,:)-xWell(1)) <= 0.5*sizeDom(1)/meth.nElDir(1)) & (abs(x(2,:)-xWell(2)) <= 0.5*sizeDom(2)/meth.nElDir(2)) ) );
muFact = T/S;
prbl.mu = @(x)( muFact*ones(1,size(x,2)) );
prbl.sigma = @(x)( ones(1,size(x,2)) );
prbl.u_0 = @(x)( 0*x(1,:) );

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
output.vtk_flag = 0;
output.vtk_name = ['VTK/output/theis_', meth.timeScheme, '_'];
output.vtk_freq = 1;
output.lmA_flag = 0;
output.smA_flag = 0;

fprintf('Completed: ')
[u_h, data, ~, space] = solve_Heat2d( prbl, meth, output);
fprintf('\n')
dt = data.dtHist(1,data.dtHist(2,:)>0);
tn = cumsum([0, dt]);

%% Plot Theis Curves
alpha_0 = 0;
alpha_1 = pi/2;
r = linspace(0,20,100); r = r(2:end);

xDof = reshape(space_dofCoord(space,geo),2,1,space.nDof);

xObsP_0 = bsxfun(@plus,bsxfun(@times,r,[cos(alpha_0);sin(alpha_0)]),xWell);
xObsP_1 = bsxfun(@plus,bsxfun(@times,r,[cos(alpha_1);sin(alpha_1)]),xWell);
[~,dofObsP_0] = min(reshape(sum(bsxfun(@minus,xObsP_0,xDof).^2,1),size(xObsP_0,2),space.nDof),[],2);
[~,dofObsP_1] = min(reshape(sum(bsxfun(@minus,xObsP_1,xDof).^2,1),size(xObsP_1,2),space.nDof),[],2);
% unique(.) without sorting
[dof_sp_0, ind_sp_0] = sort(dofObsP_0); UVp_0(ind_sp_0) = ([1; diff(dof_sp_0)] ~= 0); dofObsP_0 = dofObsP_0(UVp_0);
[dof_sp_1, ind_sp_1] = sort(dofObsP_1); UVp_1(ind_sp_1) = ([1; diff(dof_sp_1)] ~= 0); dofObsP_1 = dofObsP_1(UVp_1);
r_0 = r(UVp_0);
r_1 = r(UVp_1);

xObsM_0 = bsxfun(@plus,bsxfun(@times,-r_0,[cos(alpha_0);sin(alpha_0)]),xWell);
[~,dofObsM_0] = min(reshape(sum(bsxfun(@minus,xObsM_0,xDof).^2,1),size(xObsM_0,2),space.nDof),[],2);
xObsM_1 = bsxfun(@plus,bsxfun(@times,-r_1,[cos(alpha_1);sin(alpha_1)]),xWell);
[~,dofObsM_1] = min(reshape(sum(bsxfun(@minus,xObsM_1,xDof).^2,1),size(xObsM_1,2),space.nDof),[],2);

figure
indTn = 10:13;
tplot = tn(indTn);
for iT = 1:4
    refVals_1 = (Q/(4*pi*T))*expint((S*r_1.^2)./(4*T*tplot(iT)));
    subplot(2,2,iT)
    plot(r_1,-refVals_1,'k','LineWidth',2)
    hold on
    plot(r_0,u_h(dofObsP_0,indTn(iT)),'b','LineWidth',2)
    plot(r_1,u_h(dofObsP_1,indTn(iT)),'r','LineWidth',2)
    
    plot(-r_1,-refVals_1,'k','LineWidth',2)
    plot(-r_0,u_h(dofObsM_0,indTn(iT)),'b','LineWidth',2)
    plot(-r_1,u_h(dofObsM_1,indTn(iT)),'r','LineWidth',2)
    hold off
    title(sprintf('Time tn = %g [s]',round(10*tplot(iT))/10))
    legend('Theis','h-h_0 (0)','h-h_0 (\pi/2)','Location','SouthEast')
end

end

