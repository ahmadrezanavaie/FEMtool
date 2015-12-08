%% Convergence Test
%Tests the convergence rate of the numerical solution under h-refinement
%against the theoretical result reading
%
%   ||u-u_h||_{H^s} ~ h^{r-s+1},
%
%where ||.||_{H^s} denotes the usual norm on the Hilbert space H^s and r is
%the polynomial degree of the FE basis functions. 
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
nElDir = [8, 16, 32, 64; 8, 16, 32, 64];
meth.nQRDir = [6; 6];
rBFDir = [1, 2, 3, 4;1, 2, 3, 4];

%% Exact Solution
uEx = @(x)( sin( x(1,:)*pi ).*cos( x(2,:)*pi ) );
uExGrad = @(x)( [ pi*cos( x(1,:)*pi ).*cos( x(2,:)*pi );...
    -pi*sin( x(1,:)*pi ).*sin( x(2,:)*pi ) ] );

%% Problem Parameters
% Forcing Term and Diffusion Coefficient
prbl.f = @(x)( (2*pi.^2)*uEx(x) );
prbl.mu = @(x)( ones(1,size(x,2)) );

% Boundary Conditions
prbl.BC.neum_sides = [ 1, 4];
prbl.BC.neum_lim = [ 0, 0; 0.75, 1];
prbl.BC.neum_fun{1} = @(x)( -(pi*cos( x(1,:)*pi ).*cos( x(2,:)*pi ))./sqrt(1+pi^2*(cos(pi*x(2,:)).^2)) -...
    (pi*sin( x(1,:)*pi ).*sin( x(2,:)*pi )).*(pi*cos(pi*x(2,:)))./sqrt(1+pi^2*(cos(pi*x(2,:)).^2))  );
prbl.BC.neum_fun{2} = @(x)( -pi*sin( x(1,:)*pi ).*sin( x(2,:)*pi )  );

prbl.BC.dir_sides = [1, 2, 3];
prbl.BC.dir_lim = [ 0.75, 0, 0; 1, 1, 1];
prbl.BC.dir_fun{1} = @(x)( uEx(x) );
prbl.BC.dir_fun{2} = @(x)( uEx(x) );
prbl.BC.dir_fun{3} = @(x)( uEx(x) );

%% Solve Problems and Compute Errors
errL2 = zeros(size(rBFDir,2),size(nElDir,2));
errH1 = zeros(size(rBFDir,2),size(nElDir,2));
fprintf('\nh-refinement\n------------\n')
for iOrder = 1:size(rBFDir,2)
    fprintf('   r = %i, nQR = %i\n',min(rBFDir(:,iOrder)),min(meth.nQRDir))
    for iH = 1:size(nElDir,2)
        fprintf('   - h/H = 1/%i, nDof = %i\n',min(nElDir(:,iH)),prod(rBFDir(:,iOrder).*nElDir(:,iH)+1))
        
        % Solve Problem
        meth.nElDir = nElDir(:,iH);
        meth.rBFDir = rBFDir(:,iOrder);
        [u_h, mesh, space] = solve_Laplace2d( prbl, meth);

        % Compute Errors
        [errL2(iOrder,iH), errH1(iOrder,iH)] = spErrSol( space, mesh, geo, u_h, uEx, uExGrad );
        
%         save errL2_Laplace2D errL2
%         save errH1_Laplace2D errH1
    end
end

%% Plot Convergence
fontSize = 12;

refValX = [2,2.2];
refValY = 1./[2^4, 2^(4.2)];

fact = [5.5, 2.6, 2, 1.5];
figure
semilogy(errL2(1,:),'LineWidth',2)
hold on
semilogy(errL2(2,:),'LineWidth',2)
semilogy(errL2(3,:),'LineWidth',2)
semilogy(errL2(4,:),'LineWidth',2)

semilogy(refValX,fact(1)*(refValY).^(rBFDir(1,1)+1),'k')
semilogy(refValX,fact(2)*(refValY).^(rBFDir(1,2)+1),'k')
semilogy(refValX,fact(3)*(refValY).^(rBFDir(1,3)+1),'k')
semilogy(refValX,fact(4)*(refValY).^(rBFDir(1,4)+1),'k')

text((2*refValX(1)+refValX(2))/3,0.01,'2')
text((2*refValX(1)+refValX(2))/3,0.00028,'3')
text((2*refValX(1)+refValX(2))/3,0.000013,'4')
text((2*refValX(1)+refValX(2))/3,0.0000005,'5')
hold off
title('L^2 Error')
xlabel('h/H')
legend('r = 1','r = 2','r = 3','r = 4')
set(gca,'XTick',1:size(nElDir(1,:),2))
set(gca,'XTickLabel',{'1/8','1/16','1/32','1/64'})
set(gca,'FontSize',fontSize)

fact = [12, 18, 20, 20];
figure
semilogy(errH1(1,:),'LineWidth',2)
hold on
semilogy(errH1(2,:),'LineWidth',2)
semilogy(errH1(3,:),'LineWidth',2)
semilogy(errH1(4,:),'LineWidth',2)

semilogy(refValX,fact(1)*(refValY).^(rBFDir(1,1)),'k')
semilogy(refValX,fact(2)*(refValY).^(rBFDir(1,2)),'k')
semilogy(refValX,fact(3)*(refValY).^(rBFDir(1,3)),'k')
semilogy(refValX,fact(4)*(refValY).^(rBFDir(1,4)),'k')

text((2*refValX(1)+refValX(2))/3,0.45,'1')
text((2*refValX(1)+refValX(2))/3,0.04,'2')
text((2*refValX(1)+refValX(2))/3,0.0025,'3')
text((2*refValX(1)+refValX(2))/3,0.00015,'4')
hold off
title('H^1 Error')
xlabel('h/H')
legend('r = 1','r = 2','r = 3','r = 4')
set(gca,'XTick',1:size(nElDir(1,:),2))
set(gca,'XTickLabel',{'1/8','1/16','1/32','1/64'})
set(gca,'FontSize',fontSize)

