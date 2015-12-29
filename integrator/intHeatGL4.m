function [ u_1, varargout ] = intHeatGL4( u_0, F_0, f, FE, RK, op, BC)
%INTHEATGL4 is the numerical integrator of one time step by the 4th order
%Gauss-Legendre Runge Kutta method of the heat problem:
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
%The Butcher Tableau for the 4th order Gauss-Legendre method reads:
%
%   1/2-a |  1/4   1/4-a
%   1/2+a | 1/4+a   1/4             with a = sqrt(3)/6
%   -----------------------
%         |  1/2    1/2
%
%For minimizing the computational costs we compute a low order
%approximation by using a time integration scheme defined by the following
%Butcher tableau:
%
%      0  |   0      0      0      0
%   1/2-a |   0     1/4   1/4-a    0
%   1/2+a |   0    1/4+a   1/4     0        with a = sqrt(3)/6
%      1  |   0     1/2    1/2     0
%   --------------------------------
%         |  b_1    b_2    b_3    b_4
%
%where b_1, b_2, b_3 and b_4 are chosen such that the method is of order 3.
%It is easy to see that the quantities computed in the 4th order
%Gauss-Legendre method can be reused for the lower order approximations.
%The additional effort consists of computing the right hand side at time
%t_1 (re-use it as F_0 in the next time step).
%
%
% -------------------------------------------------------------------------
%
%   u_1 = intHeatGL4( u_0, F_0, f, FE, RK, op, BC);
%   [u_1, outData] = intHeatGL4( u_0, F_0, f, FE, RK, op, BC);
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

a = sqrt(3)/6;
BA = [1/4, 1/4-a; 1/4+a, 1/4];
Bb = [1/2, 1/2];
Bc = [1/2-a, 1/2+a];

t_12ma = RK.t+h*Bc(1);
t_12pa = RK.t+h*Bc(2);

F_12ma = rhsF( t_12ma, f, op, FE ) +...
     rhsDir( t_12ma, op.A, FE.space, BC ) +...
     rhsNeum( t_12ma, FE.space, BC );

F_12pa = rhsF( t_12pa, f, op, FE ) +...
     rhsDir( t_12pa, op.A, FE.space, BC ) +...
     rhsNeum( t_12pa, FE.space, BC );
 
Mint = op.M(intDof,intDof);
Aint = op.A(intDof,intDof);

Mat = [Mint + h*BA(1,1)*Aint, h*BA(1,2)*Aint;...
    h*BA(2,1)*Aint, Mint+h*BA(2,2)*Aint];
B = [-Aint*u_0(intDof) + F_12ma(intDof);...
    -Aint*u_0(intDof) + F_12pa(intDof)];

K = Mat\B;

[rhsDir_1, dirWeights_1] = rhsDir( RK.t+h, op.A, FE.space, BC );
u_1 = 0*u_0;
u_1(intDof) = u_0(intDof) + h*(Bb(1)*K(1:nIntDof)+Bb(2)*K((nIntDof+1):end));
u_1(dirDof) = dirWeights_1;

if nargout == 2
    % Define b_i such that the order conditions up to order 3 are satisfied
    Bbp1 = 0.25; % b_1 for the low order approximation
    
    order3Mat = [ 1, 1, 1;...
        0.5-a, 0.5+a, 1;...
        (0.5-a)^2, (0.5+a)^2, 1;...
        0.25*(0.5-a)+(0.25-a)*(0.5+a), (0.25+a)*(0.5-a)+0.25*(0.5+a), 0.5];
    Bbp = (order3Mat'*order3Mat)\(order3Mat'*[1-Bbp1;1/2;1/3;1/6]);
    Bbp = [Bbp1; Bbp]; % [b_1; b_2; b_3; b_4] for the low order approximation
    
    % Compute low order approximation
    F_1 = rhsF( RK.t+h, f, op, FE ) +...
        rhsDir_1 + rhsNeum( RK.t+h, FE.space, BC );
    
    K2p = K(1:nIntDof);
    K3p = K((nIntDof+1):end);
    
    u_1p = 0*u_0;
    u_1p(intDof) = u_0(intDof) +...
        h*Bbp(1)*( Mint\(-Aint*u_0(intDof) + F_0(intDof)) ) +...
        h*Bbp(2)*K2p + h*Bbp(3)*K3p +...
        h*Bbp(4)*( Mint\(-(Aint*( u_0(intDof) + 0.5*h*(K2p + K3p))) + F_1(intDof)) ) ;
    u_1p(dirDof) = dirWeights_1;

    outData.u_1p = u_1p;
    outData.F_1 = F_1;
    outData.orderp = 3;
    varargout{1} = outData;
end

end

