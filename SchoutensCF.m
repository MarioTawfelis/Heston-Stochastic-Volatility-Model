function y = SchoutensCF(phi,kappa,theta,rho,sigma,t,S,K,r,q,v0,Pnum)

% This function provides the characteristic function as proposed by Schoutens et al. (2004)
% ----------------------------------------------------------------------------------------s
% Inputs:
%   phi = Integration variable
%   kappa  = Heston parameter: mean reversion speed
%   theta  = Heston parameter: mean reversion level
%   rho    = Heston parameter: correlation
%   sigma  = Heston parameter: volatility of vol
%   v0     = Heston parameter: initial variance
%   S = Market parameter: Spot price
%   K = Market parameter: Strike
%   T = Market parameter: Time to maturity
%   r = Market parameter: Risk free rate
%   q = Market parameter: Dividend rate
%   Pnum: 1 or 2 (for the probabilities)
% -------------------------------------------------------
% Outputs:
%   y: Real part of the characteristic funtion
% -------------------------------------------------------
% Dependencies:
%   This is an independent function
% -------------------------------------------------------

% Log of the stock price
x = log(S);

% Forward
F = x + (r-q)*t;

% Auxillary variables
a = kappa*theta;

if Pnum==1
    u = 0.5;
    b = kappa - sigma*rho;
else
    u = -0.5;
    b = kappa;
end

d = sqrt((sigma*rho*1i*phi - b)^2 - sigma^2*(2*u*1i*phi - phi^2));
g1 = (b - sigma*rho*1i*phi + d) / (b - sigma*rho*1i*phi - d);
g2 = 1/g1;

D = (b - sigma*rho*1i*phi - d)/sigma^2*((1-exp(-d*t))/(1-g2*exp(-d*t)));
G = (1 - g2*exp(-d*t))/(1-g2);
C = F*1i*phi + a/sigma^2*((b - sigma*rho*1i*phi - d)*t - 2*log(G));

% The characteristic function
f = exp(C + D.*v0 + 1i.*phi.*x);

% Return the real part of the integrand
y = real((exp(-1i*phi*log(K/S))*f) / (1i*phi));
end

