%% COMP0043 - Numerical Methods for Finance
%   Heston Stochastic Volatility Method

%% Initialisation
% Monte Carlo parameters
nsamples = 10000; % Number of time steps/samples
nblocks = 20000;  % Number of paths

% Time grid
T = 1;                     % Maturity
t = 0;                      % Time now
dt = T/nsamples;            % Time step

% Market parameters as per Table 1 in Cui et al. (2017)
S0 = 1;                   % Initial stock price
K = 1.1;                    % Strik price
r = 0.02;                   % Interest rate
q = 0;                      % Dividend rate
tau = T-t;                  % Time to maturitay

% Model parameters as per Table 1 in Cui et al. (2017)
sigma = 0.25;               % Volatility of Volatility
kappa = 3;                  % Mean-reversion rate
theta = 0.1;                % Long-term variance
v0 = 0.08;                  % Initial variance
rho = -0.8;                 % Correlation between the BMs                         

% Check Feller condition is satisfied
assert((2*kappa * theta > sigma^2),"Feller condition not satisfied!")

%% Monte Carlo
% Initialisation
tic
S = [S0*ones(nsamples,1) zeros(nsamples,nblocks)]; % Stock paths initialised with S0
v = [v0*ones(nsamples,1) zeros(nsamples,nblocks)]; % Variance paths initialised with v0

call = zeros(nblocks,1);
put = zeros(nblocks,1);

% Generate correlated normal random variables
W1 = randn(nsamples,nblocks);
W2 = rho * W1 + sqrt(1-rho^2) * randn(nsamples,nblocks);

MCInitialisation = toc;

%% Simulation
tic
for i = 1:nblocks  
    
    % Geometric Brownian Motion - Price with stochastic volatility goverened by v(:,:)
    S(:,i+1) = S(:,i) .*  (1 + r*dt + sqrt(v(:,i)*dt) .*  W1(:,i)); 
    
    % Cox-Ignersoll-Ross square root process - Stochastic volatility
    v(:,i+1) = v(:,i) + kappa.*(theta - v(:,i)).*dt + sigma .* sqrt(v(:,i)*dt) .* W2(:,i);
    v(:,i+1) = max(v(:,i+1),zeros(nsamples,1));
    
    % Discounted payoff
    call(i) = exp(-r*T) * mean(max(S(:,i+1)-K,0));
    put(i) = exp(-r*T) * mean(max(K-S(:,i+1),0));
            
end
MCSimulation = toc;

%% Pricing
% Discounted expected payoff
MC_call = mean(call);
MC_put = mean(put);

% Standard deviation of discounted expected payoff
MC_call_stdev = sqrt(var(call)/nblocks);
MC_put_stdev = sqrt(var(put)/nblocks);

%% Numerical Methods 
% Integration setup 
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
tic
% Schoutens et al. (2004)
% Integration using the Gauss-Laguerre - Schoutens et al (2004)
[SchoutensCall1, SchoutensPut1] = ...
        FFTIntegration(S0,K,r,q,tau,kappa,theta,rho,sigma,v0,x,w,u,du,0,'quad');
    
% Integration using trapz - Schoutens et al (2004)
[SchoutensCall2, SchoutensPut2] = ...
        FFTIntegration(S0,K,r,q,tau,kappa,theta,rho,sigma,v0,x,w,u,du,0,'trapz');
    
% Cui et  al. (2017)
% Integration using the Gauss-Laguerre - Cui et al (2017)
[CuiCall1, CuiPut1] = ...
    FFTIntegration(S0,K,r,q,tau,kappa,theta,rho,sigma,v0,x,w,u,du,1,'quad');

% Integration using trapz - Cui et al (2017)
[CuiCall2, CuiPut2] = ...
    FFTIntegration(S0,K,r,q,tau,kappa,theta,rho,sigma,v0,x,w,u,du,1,'trapz');

NMSimulation = toc;

%% Output
disp(' ')
fprintf('Method                               Call        Put\n')
fprintf('------------------------------------------------------\n')
fprintf('Monte Carlo                         %5.4f      %5.4f   \n', MC_call,MC_put);
fprintf('with standard deviation of          %5.4f      %5.4f   \n', MC_call_stdev,MC_put_stdev);
disp(' ')
fprintf('Schoutens et al. (2004):\n')
fprintf('          Gauss-Laguerre            %5.4f      %5.4f \n', SchoutensCall1,SchoutensPut1);
fprintf('          Trapz                     %5.4f      %5.4f \n', SchoutensCall2,SchoutensPut2);
fprintf('Cui et al. (2017):\n')
fprintf('          Gauss-Laguerre            %5.4f      %5.4f \n', CuiCall1,CuiPut1);
fprintf('          Trapz                     %5.4f      %5.4f \n', CuiCall2,CuiPut2);
fprintf('------------------------------------------------------\n')
disp(' ')
fprintf('Absolute Error                       Call        Put\n')
fprintf('------------------------------------------------------\n')
fprintf('Schoutens-Cui:                      %1.4f      %1.4f\n', SchoutensCall1-CuiCall1,SchoutensPut1-CuiPut1);
fprintf('Schoutens-Monte Carlo               %5.4f      %5.4f\n', SchoutensCall1-MC_call,SchoutensPut1-MC_put);
fprintf('Monte Carlo-Cui:                    %1.4f      %1.4f\n', MC_call-CuiCall1,MC_put-CuiPut1);
disp(' ')