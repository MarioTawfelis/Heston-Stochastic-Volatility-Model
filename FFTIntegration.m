function [CallValue, PutValue] = FFTIntegration(S,K,r,q,T,kappa,theta,rho,sigma,v0,x,w,u,du,CF,method)

% This function provides numerical methods for integration
% -------------------------------------------------------
% Inputs:
%   S = Market parameter: Spot price
%   K = Market parameter: Strike
%   T = Market parameter: Time to maturity
%   r = Market parameter: Risk free rate
%   q = Market parameter: Dividend rate
%   kappa  = Heston parameter: mean reversion speed
%   theta  = Heston parameter: mean reversion level
%   rho    = Heston parameter: correlation
%   sigma  = Heston parameter: volatility of vol
%   v0     = Heston parameter: initial variance
%   x = Integration parameter: Gauss Laguerre abscissas
%   w = Integration parameter: Gauss Laguerre weights
%   CF: 0 = Schoutens et al (2004) Characteristic function
%       1 = Cui et al (2017) Characteristic function
% -------------------------------------------------------
% Outputs:
%   CallValue: Heston call price
%   PutValue: Heston put price
% -------------------------------------------------------
% Dependencies:
%   SchoutensCF(phi,kappa,theta,rho,sigma,tau,S,K,r,q,v0,Pnum)
%   CuiCF(phi,kappa,theta,rho,sigma,tau,S,K,r,q,v0,Pnum)
% -------------------------------------------------------

if CF == 0      
    % Schoutens et al (2004)
    if strcmp(method, 'quad')
        % Initialise results vectors to improve performance
        int1 = zeros(length(x),1);
        int2 = zeros(length(x),1);
        
        for k=1:length(x)
            int1(k) = w(k)*SchoutensCF(x(k),kappa,theta,rho,sigma,T,S,K,r,q,v0, 1);    % P1
            int2(k) = w(k)*SchoutensCF(x(k),kappa,theta,rho,sigma,T,S,K,r,q,v0, 2);    % P2
        end
        
        % The integrals
        I1 = sum(int1);
        I2 = sum(int2);
    elseif strcmp(method, 'trapz')
        % Initialise results vectors to improve performance
        int1 = zeros(length(u),1);
        int2 = zeros(length(u),1);

        for k=1:length(u)
            int1(k) = SchoutensCF(u(k),kappa,theta,rho,sigma,T,S,K,r,q,v0, 1);    % P1
            int2(k) = SchoutensCF(u(k),kappa,theta,rho,sigma,T,S,K,r,q,v0, 2);    % P2
        end
        
        % The integrals
        I1 = trapz(int1)*du;
        I2 = trapz(int2)*du;
    end
else
    % Cui et al (2017)
    if strcmp(method, 'quad')
        % Initialise results vectors to improve performance
        int1 = zeros(length(x),1);
        int2 = zeros(length(x),1);
        
        for k=1:length(x)
            int1(k) = w(k)*CuiCF(x(k),kappa,theta,rho,sigma,T,S,K,r,q,v0, 1);    % P1
            int2(k) = w(k)*CuiCF(x(k),kappa,theta,rho,sigma,T,S,K,r,q,v0, 2);    % P2
        end
        
        % The integrals
        I1 = sum(int1);
        I2 = sum(int2);
    elseif strcmp(method, 'trapz')
        % Initialise results vectors to improve performance
        int1 = zeros(length(u),1);
        int2 = zeros(length(u),1);

        for k=1:length(u)
            int1(k) = CuiCF(u(k),kappa,theta,rho,sigma,T,S,K,r,q,v0, 1);    % P1
            int2(k) = CuiCF(u(k),kappa,theta,rho,sigma,T,S,K,r,q,v0, 2);    % P2
        end
        
        % The integrals
        I1 = trapz(int1)*du;
        I2 = trapz(int2)*du;
    end
end

% Discounted payoff
% P1(?; K, ?) and P2(?; K, ?) are solutions to certain pricing PDEs (Heston, 1993, Eq. (12))
P1 = 1/2 + 1/pi*I1;
P2 = 1/2 + 1/pi*I2;

% Call price
CallValue = S*P1 - K*exp(-r*T)*P2;

% Put price using Put-Call Parity
PutValue = CallValue - S + K*exp(-r*T);
end