%% COMP0043 - Numerical Methods for Finance
%   Heston Stochastic Volatility Method
%   Monte Carlo

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

%% Output
disp('')
fprintf('Method          Call      Put   \n')
fprintf('--------------------------------\n')
fprintf('Monte Carlo     %5.4f   %5.4f   \n', MC_call,MC_put);
disp('')
fprintf("Monte Carlo initialisation execution time %3.1f:\n",MCInitialisation)
fprintf("Monte Carlo simulation execution time %3.1f:\n",MCSimulation)
fprintf("Monte Carlo total execution time %3.1f seconds:\n",MCInitialisation+MCSimulation)
disp("Monte Carlo completed.")
disp('')