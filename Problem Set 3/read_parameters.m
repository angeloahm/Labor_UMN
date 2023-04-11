%% Read parameters
% This code has all the model/numerical parameters to solve Shimer (2005)

%% Model parameters (and gridpoints)

n = 1000;                     %gridpoints for y

p_star = 1;                 %long-run productivity 
s      = 0.1;               %long-run separation rate
z0     = 0.4;               %value of unemployment
r      = 0.012;             %interest rate
gamma  = 0.004;             %mean of the OU-process 
sigma  = 0.0165;            %sd of the OU-process
lambda = n*gamma;           %poisson shock
c      = 0.213;             %cost of keep an open vacancy  
mu     = 1.355;             %parameter of the matching function
alpha  = 0.72;              %parameter of the matching function
beta0  = 0.72;              %bargaining power

%Source of shocks
%bargain = 0;

f = @(t) mu*t^(1-alpha);            %matching function
q = @(t) f(t)/t;                    %vacancy mathing rate
du = @(s, u, t) s*(1-u) - f(t)*u;   %motion for u
%% Numerical parameters 

options = optimset(@fsolve);
options.Display = 'off';

eps         = 1e-7;                   %tolerance
max_iter    = 1000;                   %maximum iterations
print       = 50;                     %print           

Delta = sigma/(lambda^(1/2));               %step size

% Grid for the latent state
y_grid = linspace(-n*Delta, n*Delta, 2*n+1);

P = zeros(length(y_grid), length(y_grid));                        %allocate transition matrix

%Simulation parameters
T = 1212;               %number of periods (quarters) 
W = 10000;              %number of simulations

%% Setup shocks
% Fill the transition matrix
for i=1:length(y_grid)
    y = y_grid(i);
        
    if i==1                     %lower reflecting barrier
        P(i, i+1)=1;            
    elseif i==length(y_grid)    %upper reflecting barrier
        P(i,i-1)=1;             
    else
        P(i,i-1) = (1/2)*(1+y/(n*Delta));
        P(i,i+1) = (1/2)*(1-y/(n*Delta));
    end
        
end

p_grid = z0 + exp(y_grid)*(p_star-z0);        %productivity grid
%s_grid = exp(y_grid)*s_star;                %separation grid



 



