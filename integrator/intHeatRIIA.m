function [ u_1, varargout ] = intHeatRIIA( u_0, F_0, f, FE, RK, op, BC)
%INTHEATRIIA is the numerical integrator of one time step by the Radau IIA
%Runge Kutta method of the heat problem:
%
% find u : \Omega x ]0,T] --> R s.t.
%
%              du/dt - grad( \mu div(u) ) = f(x,t),     in \Omega x ]0,T]
%                                  u(.,0) = u_0,        in \Omega
%                                  u(x,t) = g(x,t),     on \Gamma_D x ]0,T]
%                          \mu du(x,t)/dn = h(x,t),     on \Gamma_N x ]0,T]
%
%where \Omega is a compact subset of R^2 parametrizable from a parametric
%rectangular domain.
%
%The Butcher Tableau for the Radau IIA method reads:
%
%   1/3 | 5/12  -1/12
%     1 |  3/4    1/4
%   -----------------------
%       |  3/4    1/4
%
%For minimizing the computational costs we compute a low order
%approximation by using a time integration scheme defined by the following
%Butcher tableau:
%
%     0 |  0     0      0
%   1/3 |  0    5/12  -1/12
%     1 |  0    3/4    1/4
%   -----------------------
%       | b_1   b_2    b_3
%
%where b_1, b_2 and b_3 are chosen such that the method is of order 2. It
%is easy to see that the quantities computed in the Radau IIA method can be
%reused for the lower order approximations. The additional effort consists
%of computing the right hand side at time t_0. 
%
%
% -------------------------------------------------------------------------
%
%   u_1 = intHeatRIIA( u_0, F_0, f, FE, RK, op, BC);
%   [u_1, outData] = intHeatRIIA( u_0, F_0, f, FE, RK, op, BC);
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% u_0               (nDof x 1 double) Dof-weights of the solution at t_0
% F_0               (nDof x 1 double) Reaction term plus boundary condition
%                   weights at t_0
% f                 (function handle) Function handle of the reaction term
% FE                (struct)
%   .geo            (see geo_2d.m)
%   .mesh           (see mesh_2d.m)
%   .space          (see space_2d.m)
% RK                (struct)
%   .t              (1 x 1 double) Time t_0
%   .dt             (1 x 1 double) Time step dt: t_1 = t_0 + dt
% op                (struct)
%   .A              (nDof x nDof double) Stiffness matrix
%   .M              (nDof x nDof double) Mass matrix
%   (for constant reaction term)
%   .f              (nDof x 1 double) Reaction weights
% BC                (struct)
%	.dir_sides      (1 x nDir double) Dirichlet boundary indices
%	.dir_lim        (2 x nDir double) Dirichlet limits
%	.dir_fun        (1 x nDir cell array) Dirichlet function handles
%	.neum_sides     (1 x nNeum double) Neumann boundary indices
% |-(if .neum_sides non-empty)
% | .neum_lim       (2 x nNeum double) Neumann limits
% |-.neum_fun       (1 x nNeum cell array) Neumann function handles
%   .dir            (see solve_Heat2d:bndry_info)
%   .neum           (see solve_Heat2d:bndry_info)
%
%
% OUTPUT
% ------
% u_1               (nDof x 1 double) Dof-weights of the solution at t_1
% (optional)
% outData           (struct)
%   .F_1            (nDof x 1 double) Reaction term plus boundary condition
%                   weights at t_1
%   .u_1p           (nDof x 1 double) Dof-weights of a lower order solution
%                   at t_1
%   .orderp         (1 x 1 integer) Lower order
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

h = RK.dt;

dirDof = BC.dir.dof;
intDof = setdiff(1:FE.space.nDof,BC.dir.dof);
nIntDof = length(intDof);

BA = [5/12, -1/12; 3/4, 1/4];
Bb = [3/4, 1/4];
Bc = [1/3, 1];

t_13 = RK.t+h*Bc(1);
t_1 = RK.t+h;

F_13 = rhsF( t_13, f, op, FE ) +...
     rhsDir( t_13, op.A, FE.space, BC ) +...
     rhsNeum( t_13, FE.space, BC );

[rhsDir_1, dirWeights_1] = rhsDir( t_1, op.A, FE.space, BC );
F_1 = rhsF( t_1, f, op, FE ) +...
     rhsDir_1 + rhsNeum( t_1, FE.space, BC );
 
Mint = op.M(intDof,intDof);
Aint = op.A(intDof,intDof);

Mat = [Mint + h*BA(1,1)*Aint, h*BA(1,2)*Aint;...
    h*BA(2,1)*Aint, Mint+h*BA(2,2)*Aint];
B = [-Aint*u_0(intDof) + F_13(intDof);...
    -Aint*u_0(intDof) + F_1(intDof)];

K = Mat\B;

u_1 = 0*u_0;
u_1(intDof) = u_0(intDof) + h*(Bb(1)*K(1:nIntDof)+Bb(2)*K((nIntDof+1):end));
u_1(dirDof) = dirWeights_1;

if nargout == 2
    % Define b_i such that the order conditions up to order 2 are satisfied
    Bbp1 = 1/3; % b_1 for the low order approximation
    
    order2Mat = [ 1, 1; 1/3, 1];
    Bbp = order2Mat\[1-Bbp1;1/2];
    Bbp = [Bbp1; Bbp]; % [b_1; b_2; b_3] for the low order approximation
    
    % Compute low order approximation
    K2p = K(1:nIntDof);
    K3p = K((nIntDof+1):end);
    
    u_1p = 0*u_0;
    u_1p(intDof) = u_0(intDof) + h*Bbp(1)*( Mint\(-Aint*u_0(intDof) + F_0(intDof)) ) +...
        h*Bbp(2)*K2p + h*Bbp(3)*K3p;
    u_1p(dirDof) = dirWeights_1;

    outData.u_1p = u_1p;
    outData.F_1 = F_1;
    outData.orderp = 2;
    varargout{1} = outData;
end

end

