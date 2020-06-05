%% COMP0043 - Numerical Methods for Finance
%   Heston Stochastic Volatility Method
%   Numerical Methods - Fourier

%% Initialisation

% Time grid
T = 1;              % Maturity
t = 0;            % Time now                       

% Market parameters as per Table 1 in Cui et al. (2017)
S0 = 1;            % Initial stock price
K = 1.1;          % Strik price
r = 0.02;         % Interest rate
q = 0;            % Dividend rate
tau = T-t;        % Time to maturity

% Model parameters as per Table 1 in Cui et al. (2017)
sigma = 0.25;               % Volatility of Volatility
kappa = 3;                  % Mean-reversion rate
theta = 0.1;                % Long-term variance
v0 = 0.08;                  % Initial variance
rho = -0.8;                 % Correlation between the BMs                         

% Check Feller condition is satisfied
assert((2*kappa * theta > sigma^2),"Feller condition not satisfied!")

%% Integration setup 
tic
% Weights and abscissas for Gauss Laguerre integration
[x, w] = GenerateGaussLaguerre(32);

% Integration grid for integration using MATLAB's trapz				 
Lu = 1e-10;           % Lower limit
du = 0.01;            % Increment
Uu = 100;             % Upper limit
u = [Lu:du:Uu]; % Grid
NMInitialisation = toc;

%% Pricing
%% Schoutens et al. (2004)
tic

% Integration using the Gauss-Laguerre - Schoutens et al (2004)
[SchoutensCall1, SchoutensPut1] = ...
        FFTIntegration(S0,K,r,q,tau,kappa,theta,rho,sigma,v0,x,w,u,du,0,'quad');
    
% Integration using trapz - Schoutens et al (2004)
[SchoutensCall2, SchoutensPut2] = ...
        FFTIntegration(S0,K,r,q,tau,kappa,theta,rho,sigma,v0,x,w,u,du,0,'trapz');
    
%% Cui et  al. (2017)

% Integration using the Gauss-Laguerre - Cui et al (2017)
[CuiCall1, CuiPut1] = ...
    FFTIntegration(S0,K,r,q,tau,kappa,theta,rho,sigma,v0,x,w,u,du,1,'quad');

% Integration using trapz - Cui et al (2017)
[CuiCall2, CuiPut2] = ...
    FFTIntegration(S0,K,r,q,tau,kappa,theta,rho,sigma,v0,x,w,u,du,1,'trapz');

NMSimulation = toc;

%% Output

disp(' ')
fprintf('Characteristic function           Call        Put\n')
fprintf('----------------------------------------------------\n')
fprintf('Schoutens et al (2004):\n')
fprintf('\t\tGauss-Laguerre    %5.4f      %5.4f \n', SchoutensCall1,SchoutensPut1);
fprintf('\t\tTrapz             %5.4f      %5.4f \n', SchoutensCall2,SchoutensPut2);
fprintf('Cui et al (2017):\n')
fprintf('\t\tGauss-Laguerre    %5.4f      %5.4f \n', CuiCall1,CuiPut1);
fprintf('\t\tTrapz             %5.4f      %5.4f \n', CuiCall2,CuiPut2);
disp(' ')
fprintf('Schoutens-Cui price error: \t  %1.4f      %1.4f\n',SchoutensCall1-CuiCall1,SchoutensPut1-CuiPut1);
disp(' ')
fprintf("Numerical methods initialisation execution time %3.1f\n:",NMInitialisation)
fprintf("Numerical methods simulation execution time %3.1f\n:",NMSimulation)
fprintf("Numerical methods total execution time %3.1f\n:",NMInitialisation+NMSimulation)
disp("Numerical methods completed.")