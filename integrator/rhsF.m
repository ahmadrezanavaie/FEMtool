function [ rhs ] = rhsF( t, f, op, FE )
%RHSF computes the right hand side weights of the reaction term. If it is
%stationary, previously computed weights are re-used.
%
%
% -------------------------------------------------------------------------
%
%   rhs = rhsF( t, f, op, BC );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% t                 (1 x 1 double) Time t
% f                 (function handle) Function handle of the reaction term
% op                (struct)
%   (for constant reaction term)
%   .f              (nDof x 1 double) Reaction weights
% FE                (struct)
%   .geo            (see geo_2d.m)
%   .mesh           (see mesh_2d.m)
%   .space          (see space_2d.m)
%
%
% OUTPUT
% ------
% rhs               ( nDof x 1 double) Right hand side weights associated
%                   to the reaction term
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

if nargin(f) == 2
    fun = @(x)( f(x,t) );
    rhs = opFV(FE.space,FE.mesh,fun,FE.geo);
else
    rhs = op.f;
end

end

