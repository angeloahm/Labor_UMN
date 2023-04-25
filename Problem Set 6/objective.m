function g = objective(ui)

% Read parameters of the model
read_parameters;

% Solve model with a specific value for the unemployment insurance
b_grid = ui * w_grid;
solve_model


%Simulate the model 
simulate

%Define the objective
u_target = 0.065;
u_sim = 1-mean(e_sim);
h = u_sim - u_target;

g = (1/T) * sum(h(:).^2);

fprintf('Objective function is:                %0.8f \n', g)
write = [ui g];
dlmwrite('result.csv',write,'delimiter',',','-append');


end