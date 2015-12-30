function [ u_1, varargout ] = intHeatEI( u_0, F_0, f, FE, RK, op, BC, varargin)
%INTHEATEI is the numerical integrator of one time step by the EULER
%IMPLICIT Runge Kutta method of the heat problem: 
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
%
% -------------------------------------------------------------------------
%
%   u_1 = intHeatEI( u_0, F_0, f, FE, RK, op, BC, ...);
%   [u_1, outData] = intHeatEI( u_0, F_0, f, FE, RK, op, BC, ...);
%   [u_1, ...] = intHeatEI( u_0, F_0, f, FE, RK, op, BC);
%   [u_1, ...] = intHeatEI( u_0, F_0, f, FE, RK, op, BC, rhsDir_1, wDir_1, F_1);
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
% (optional)
% rhsDir_1          (nDof x 1) Right hand side weights associated to the
%                   Dirichlet boundary condition at t_1 (on the "internal"
%                   nodes it is given by -A(int,dir)*wDir_1 )
% wDir_1            (nDirDof x 1) Dirichlet projection weights at t_1
% F_1               (nDof x 1 double) Reaction term plus boundary condition
%                   weights at t_1
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

nArgIn0 = -nargin('intHeatEi')-1;

t_1 = RK.t+RK.dt;
dirDof = BC.dir.dof;
intDof = setdiff(1:FE.space.nDof,BC.dir.dof);

if nargin - nArgIn0 == 3
    rhsDir_1 = varargin{1};
    wDir_1 = varargin{2};
    F_1 = varargin{3};
else
    [rhsDir_1, wDir_1] = rhsDir( t_1, op, FE.space, BC );
    F_1 = rhsF( t_1, f, op, FE ) + rhsDir_1 + rhsNeum( t_1, FE.space, BC );
end
 
Mint = op.M(intDof,intDof);
Aint = op.A(intDof,intDof);

u_1 = 0*u_0;
u_1(intDof) = (Mint + RK.dt*Aint)\(Mint*u_0(intDof) + RK.dt*F_1(intDof));
u_1(dirDof) = wDir_1;

if nargout == 2
    outData.F_1 = F_1;
    outData.u_1p = intHeatEE( u_0, F_0, f, FE, RK, op, BC, rhsDir_1, wDir_1);
    outData.orderp = 1;
    varargout{1} = outData;
end

end