
function out = budget_balance(taxes)
% This code computes the government budget balance
global dist

read_parameters;

%Define taxes
tau=taxes;

%Solve model for a given tau
fprintf('----------- Solving Model -----------\n')
solve_model;


%Compute stationary distribution
fprintf('----------- Running Simulation -----------\n')
simulation;

%Compute budget balance
w_matrix = w_grid' * ones(1,n_h);
h_matrix = ones(n_w,1) * h_grid;

help = (1-policy_employed) .* dist;
unemp = sum( help(:) );

revenue = sum( tau*b.*(1-policy_employed(:)).*dist(:) ) + sum(tau*policy_employed(:).*dist(:).*w_matrix(:).*h_matrix(:));
cost = sum( b.*(1-policy_employed(:)).*dist(:) );



%Output for function:
out = revenue-cost;



fprintf('tau:               %.4f \n', tau);
fprintf('Excess revenue:    %.4f \n', out);
fprintf('Unemp. rate:       %.4f \n', unemp);
