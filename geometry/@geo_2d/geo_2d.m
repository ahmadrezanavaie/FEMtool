classdef geo_2d
%GEO_2D is a class that defines a computational domain \Omega trough a
%bijective and smooth mapping function F starting from a rectangular
%parametric domain [a,b]x[c,d] into the computational domain \Omega:
%
%                         (b,d)
%        ___________________
%       |                   |
%       |                   |
%       |                   |    F
%       |                   |   --->    \Omega \subset \mathbb{R}^2
%       |                   |
% p2    |                   |
% ^     |___________________|
% |   (a,c)
% |
%-|---> p1 
%
%
% -------------------------------------------------------------------------
%
%   geo = geo_2d( );
%   geo = geo_2d( 'geometry' );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% 'geometry' can be a user-specified string that refers to its own geometry
% defined in the constructor.
%
%
% OUTPUT
% ------
% geo           (geo_2d object)
%   .paramDom   (1 x 1 struct)
%       .origin (2 x 1 double) Defines the origin of the parametric domain
%       .sizeDir(2 x 1 double) Defines the size of the parametric domain
%               in each parametric direction
%   .map        (function handle) Definition of the parameterization
%               function F
%   .Jac        (function handle) Definition of the Jacobian matrix of F
%   .detJac     (function handle) Definition of the determinant of the
%               Jacobian matrix of F
%   .JacInvT    (function handle) Definition of the Jacobian matrix of the
%               inverse map (for the gradient transformation)
%
%
% Remarks:
% --------
% 1) The field .Jac is a function R^d -> R^(dxd) defined in the PARAMETRIC
% DOMAIN. If x is of size (d x N) then .Jac(x) is of size (d x d x N).
%
% 2) The field .detJac is a function R^d -> R defined in the PARAMETRIC
% DOMAIN and denotes the determinant of the Jacobian matrix of F. If x is
% of size (d x N) then .detJac(x) is of size (1 x N).
%
% 3) The field .JacInvT is a function R^d -> R^(dxd) defined in the
% PHYSICAL SPACE. If we denote by JF^{-1} the Jacobian matrix of the
% INVERSE OF F, then JacInvT represents JF^{-T}. This is useful
% because if grad_x (resp. grad_p) denotes the gradient in the physical
% (resp. parameter) domain, then
%
%                grad_x = JF^{-T}grad_p
%
% and
%
%     (grad_x U)^T grad_x V = (grad_p U)^T JF^{-1}*(JF^{-1})^T grad_p V.
%
% If x is of size (d x N) then .JacInvT(x) is of size (d x d x N).
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
    
    properties
        type = '2-d Geometry';
        paramDom = [];
        map = [];
        Jac = [];
        detJac = [];
        JacInvT = [];
    end
    
    methods
        function obj = geo_2d(varargin)
            if nargin == 0 || strcmp(varargin{1},'rectangle')
                % Parametric Domain
                obj.paramDom.origin = [0; 0];    % = [a;c]
                obj.paramDom.sizeDir = [1; 1];   % = [b-a;d-c]

                % Computational Domain
                if nargin == 1
                    origin = [0;0];
                    sizeDir = [100;100];
                elseif nargin == 3
                    origin = varargin{2};
                    sizeDir = varargin{3};
                else
                    error('geo_2d: Wrong number of input arguments for rectangle geometry')
                end
                obj.map = @(x)( bsxfun(@plus,bsxfun(@times,x,sizeDir./obj.paramDom.sizeDir),...
                    origin-obj.paramDom.origin.*sizeDir./obj.paramDom.sizeDir) );

                % JACOBIAN of MAP F
                obj.Jac = @(x)( repmat([sizeDir(1)/obj.paramDom.sizeDir(1), 0; 0, sizeDir(2)/obj.paramDom.sizeDir(2)],1,1,size(x,2)) );
                % DETERMINANT of JACOBIAN of MAP F
                obj.detJac = @(x)( repmat((sizeDir(1)*sizeDir(2))./(obj.paramDom.sizeDir(1)*obj.paramDom.sizeDir(2)),1,size(x,2)) );

                % TRANSPOSED JACOBIAN of INVERSE MAP of F (for the gradient change)
                % This function has to be defined over the PHYSICAL SPACE, and NOT over
                % the parameter space (as geo.detJac and geo.map)
                obj.JacInvT = @(x)( repmat([obj.paramDom.sizeDir(1)/sizeDir(1), 0; 0, obj.paramDom.sizeDir(2)/sizeDir(2)],1,1,size(x,2)) );

            elseif strcmp(varargin{1},'paral')
                a = 1;

                % Parametric Domain
                obj.paramDom.origin = [0; 0];    % = [a;c]
                obj.paramDom.sizeDir = [1; 1];   % = [b-a;d-c]

                % Computational Domain
                obj.map = @(x)( [x(1,:) + a*x(2,:); x(2,:)] );

                % JACOBIAN of MAP F
                obj.Jac = @(x)( repmat([1,a;0,1],1,1,size(x,2)) );
                % DETERMINANT of JACOBIAN of MAP F
                obj.detJac = @(x)( ones(1,size(x,2)) );

                % TRANSPOSED JACOBIAN of INVERSE MAP of F (for the gradient change)
                % This function has to be defined over the PHYSICAL SPACE, and NOT over
                % the parameter space (as geo.detJac and geo.map)
                obj.JacInvT = @(x)( repmat([1,0;-a,1],1,1,size(x,2)) );

            elseif strcmp(varargin{1},'bow')
                % Parametric Domain
                obj.paramDom.origin = [0; 0];    % = [a;c]
                obj.paramDom.sizeDir = [1; 1];   % = [b-a;d-c]

                % Computational Domain
                obj.map = @(x)( [x(1,:) + sin(x(2,:).*pi); x(2,:)] );

                % JACOBIAN of MAP F
                obj.Jac = @(x)( bsxfun( @plus, eye(2), bsxfun(@times,[0,1;0,0],reshape(pi*cos(pi*x(2,:)),1,1,size(x,2))) ) );
                % DETERMINANT of JACOBIAN of MAP F
                obj.detJac = @(x)( ones(1,size(x,2)) );

                % TRANSPOSED JACOBIAN of INVERSE MAP of F (for the gradient change)
                % This function has to be defined over the PHYSICAL SPACE, and NOT over
                % the parameter space (as geo.detJac and geo.map)
                obj.JacInvT = @(x)( bsxfun( @plus, eye(2), bsxfun(@times,[0,0;1,0],reshape(-pi*cos(pi*x(2,:)),1,1,size(x,2))) ) );

            elseif strcmp(varargin{1},'banana')
                % Parametric Domain
                obj.paramDom.origin = [0; 0];    % = [a;c]
                obj.paramDom.sizeDir = [1; 1];   % = [b-a;d-c]

                % Computational Domain
                obj.map = @(x)( [x(1,:) + sin(x(2,:).*(pi/2)); x(2,:)] );

                % JACOBIAN of MAP F
                obj.Jac = @(x)( bsxfun( @plus, eye(2), bsxfun(@times,[0,1;0,0],reshape((pi/2)*cos((pi/2)*x(2,:)),1,1,size(x,2))) ) );
                % DETERMINANT of JACOBIAN of MAP F
                obj.detJac = @(x)( ones(1,size(x,2)) );

                % TRANSPOSED JACOBIAN of INVERSE MAP of F (for the gradient change)
                % This function has to be defined over the PHYSICAL SPACE, and NOT over
                % the parameter space (as geo.detJac and geo.map)
                obj.JacInvT = @(x)( bsxfun( @plus, eye(2), bsxfun(@times,[0,0;1,0],reshape((-pi/2)*cos((pi/2)*x(2,:)),1,1,size(x,2))) ) );

            elseif strcmp(varargin{1},'trapezoid')
                % Parametric Domain
                obj.paramDom.origin = [0; 0];    % = [a;c]
                obj.paramDom.sizeDir = [1; 1];   % = [b-a;d-c]

                % Computational Domain
                obj.map = @(x)( [ (x(1,:)-0.5).*(3-2*x(2,:)); x(2,:) ] );

                % JACOBIAN of MAP F
                obj.Jac = @(x)( bsxfun(@plus, bsxfun(@times,[1,0;0,0],reshape(3-2*x(2,:),1,1,size(x,2))),...
                    bsxfun( @plus, repmat([0,0;0,1],1,1,size(x,2)),...
                    bsxfun(@times,[0,1;0,0],reshape(-2*(x(1,:)-0.5),1,1,size(x,2))) ) ) );
                % DETERMINANT of JACOBIAN of MAP F
                obj.detJac = @(x)( 3-2*x(2,:) );

                % JTRANSPOSED JACOBIAN of INVERSE MAP of F (for the gradient change)
                % This function has to be defined over the PHYSICAL SPACE, and NOT over
                % the parameter space (as geo.detJac and geo.map)
                obj.JacInvT = @(x)( bsxfun(@plus, bsxfun(@times,[1,0;0,0],reshape(1./(3-2*x(2,:)),1,1,size(x,2))),...
                    bsxfun( @plus, repmat([0,0;0,1],1,1,size(x,2)),...
                    bsxfun(@times,[0,0;1,0],reshape((2*x(1,:))./((3-2*x(2,:)).^2),1,1,size(x,2))) ) ) );

            else
                error('geo_2d: Geometry %s is not defined yet',varargin{1})
            end
        end
    end
    
end

