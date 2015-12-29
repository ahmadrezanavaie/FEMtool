function [mu, nu, mup, gammap, ct] = coeffRKC2( s, eta )
%COEFFRKC2 computes the coefficients for the second order damped Runge
%Kutta Chebyshev (RKC) method for a given stage number s and damping 
%coefficient eta. If s=1 the RKC method is equivalent to the Euler explicit
%method and therefore, the coefficients mu, nu, and gammap are all zeros.
%
%See for example:
%   B.P. Sommejer, L.F. Shampine, J.G. Verwer. RKC: An explicit solver for
%   parabolic PDEs. Journal of Computational and Applied Mathematics 88
%   (1997) 315-326.
%
%
% -------------------------------------------------------------------------
%
%   [mu, nu, mup, gammap, ct] = coeffRKC2( s, eta );
%
% -------------------------------------------------------------------------
%
%
% INPUT
% -----
% s                 (1 x 1 integer) Stage number s
% eta               (1 x 1 double) Damping coefficient eta
%
%
% OUTPUT
% ------
% mu, nu,
% mup, gammap, ct   (1 x s+1 double) Coefficients of the stabilized Runge
%                   Kutta Chebyshev method
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

mu = zeros(1,s+1); nu = mu; mup = mu; gammap = mu; ct = mu;

if s <= 0
    error('coeffRKC2: Bad stage number s = %i',s)
elseif s == 1;
    ct(2) = 1;
    mup(2) = 1;
else
    w0 = 1+eta/(s.^2);
    [Tw0, dTw0, ddTw0] = cheb(s,w0);
    w1 = dTw0(s+1)/ddTw0(s+1);
    
    ct(1) = 0; ct(end) = 1;
%     ct(3:(end-1)) = ((2:(s-1)).^2-1)./(s.^2-1);
%     ct(2) = ct(3)/4;
    ct(3:(end-1)) = (dTw0(end)/ddTw0(end))*(ddTw0(3:(end-1))./dTw0(3:(end-1)));
    ct(2) = ct(3)/dTw0(3);
    
    b = zeros(1,s+1);
    b(3:end) = ddTw0(3:end)./(dTw0(3:end).^2);
    b(2) = b(3); b(1) = b(3);
    
    mup(2) = b(2)*w1;
    mup(3:end) = (b(3:end)./(b(2:(end-1))))*2*w1;
    
    mu(3:end) = (b(3:end)./(b(2:(end-1))))*2*w0;
    
    nu(3:end) = -(b(3:end)./(b(1:(end-2))));
    
    gammap(3:end) = -(1-b(2:(end-1)).*Tw0(2:(end-1))).*mup(3:end);
    
end

end

function [T, dT, ddT] = cheb(s,x)
%CHEB computes the Chebyshev polynomial of order <= s as well as its first
%and second derivative at x

T = zeros(1,s+1); dT = T; ddT = T;

T(1) = 1; T(2) = x;
dT(2) = 1;

for ii = 2:s
    T(ii+1) = 2*x*T(ii) - T(ii-1);
    dT(ii+1) = 2*T(ii) + 2*x*dT(ii) - dT(ii-1);
    ddT(ii+1) = 4*dT(ii) + 2*x*ddT(ii) - ddT(ii-1);
end

end

