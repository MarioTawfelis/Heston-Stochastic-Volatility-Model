function y = CuiCF(phi,kappa,theta,rho,sigma,t,S,K,r,q,v0,Pnum)

% This function provides the characteristic function as proposed by Cui et al. (2017)
% ----------------------------------------------------------------------------------
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

xi = kappa - sigma*rho*1i*phi;
d = sqrt((sigma*rho*1i*phi - b)^2 - sigma^2*(2*u*1i*phi - phi^2));

A1 = (phi^2 + 1i*phi)*sinh(d*t/2);
A2 = (d/v0)*cosh(d*t/2) + (xi/v0)*sinh(d*t/2);
A = A1 / A2;
B = (d*exp(kappa*t/2)) / (v0*A2);


% Characteristic function as per equation (14)
f = exp((1i*phi*F) - ((a*rho*t*1i*phi) / sigma) - A) * B^(2*a / sigma^2);

% Return the real part of the integrand
y = real((exp(-1i*phi*log(K/S))*f) / (1i*phi));
end