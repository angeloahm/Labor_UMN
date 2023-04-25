function out = g(ui)

% Read parameters of the model
read_parameters;

% Solve model with a specific value for the unemployment insurance
b_grid = ui * w_grid;
solve_model


%Simulate the model 
simulate

%Define the objective
u_target = 0.065;
h = avg_u - u_target;

out = (1/T) * sum(h(:)^2);


end