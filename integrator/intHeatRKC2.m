function [ u_1, varargout ] = intHeatRKC2( u_0, F_0, f, FE, RKC, op, BC)
%INTHEATRKC2 is the numerical integrator of one time step by the second
%order DAMPED and explicit RUNGE KUTTA CHEBYSHEV method of the heat
%problem:
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
%See for example:
%   B.P. Sommejer, L.F. Shampine, J.G. Verwer. RKC: An explicit solver for
%   parabolic PDEs. Journal of Computational and Applied Mathematics 88
%   (1997) 315-326.
%
%
% -------------------------------------------------------------------------
%
%   u_1 = intHeatRKC2( u_0, F_0, f, FE, RK, op, BC);
%   [u_1, outData] = intHeatRKC2( u_0, F_0, f, FE, RK, op, BC);
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
% RKC               (struct)
%   .t              (1 x 1 double) Time t_0
%   .dt             (1 x 1 double) Time step dt: t_1 = t_0 + dt
%   .s              (1 x 1 integer) Stage number s
%   .eta            (1 x 1 double) Damping coefficient eta
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

t_0 = RKC.t; dt = RKC.dt;
dirDof = BC.dir.dof;
intDof = setdiff(1:FE.space.nDof,BC.dir.dof);

[mu, nu, mup, gammap, ct] = coeffRKC2(RKC.s, RKC.eta);

Aint = op.A(intDof,intDof);
Mlumped = sum(op.M,2);

R_0 = 0*u_0;
R_0(intDof) = (-Aint*u_0(intDof) + F_0(intDof))./Mlumped(intDof);

[rhsDir_jm1, wDir_jm1] = rhsDir( t_0 + ct(2)*dt, op, FE.space, BC );
g_0 = u_0; g_jm2 = g_0;
g_jm1 = 0*g_0; g_jm1(dirDof) = wDir_jm1;
g_jm1(intDof) = g_0(intDof) + mup(2)*dt*R_0(intDof);
g_j = g_jm1;

F_jm1 = rhsF( t_0+ct(2)*dt, f, op, FE ) + rhsDir_jm1 + rhsNeum( t_0+ct(2)*dt, FE.space, BC );
R_jm1 = 0*R_0;
R_jm1(intDof) = (-Aint*g_jm1(intDof) + F_jm1(intDof))./Mlumped(intDof);

for jj = 2:RKC.s
    [rhsDir_j, wDir_j] = rhsDir( t_0 + ct(jj+1)*dt, op, FE.space, BC );
    g_j = 0*g_jm1; g_j(dirDof) = wDir_j;
    g_j(intDof) = (1-mu(jj+1)-nu(jj+1))*g_0(intDof) + mu(jj+1)*g_jm1(intDof) + ...
        nu(jj+1)*g_jm2(intDof) + mup(jj+1)*dt*R_jm1(intDof) + gammap(jj+1)*dt*R_0(intDof);
    
    g_jm2 = g_jm1;
    g_jm1 = g_j;
    
    F_jm1 = rhsF( t_0+ct(jj+1)*dt, f, op, FE ) + rhsDir_j + rhsNeum( t_0+ct(jj+1)*dt, FE.space, BC );
    R_jm1(intDof) = (-Aint*g_j(intDof) + F_jm1(intDof))./Mlumped(intDof);
end

u_1 = g_j;

if nargout == 2
    outData.F_1 = F_jm1;
%     outData.u_1p = intHeatEI( u_0, F_0, f, FE, RK, op, BC, rhsDir_j, wDir_j, F_jm1);
    outData.u_1p = intHeatEE( u_0, F_0, f, FE, RKC, op, BC, rhsDir_j, wDir_j);
    outData.orderp = 1;
    varargout{1} = outData;
end

end

