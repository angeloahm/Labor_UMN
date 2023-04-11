function out = objective_calibration(guess)
%This function computes the distance between the model moments and data
%targets


read_parameters

%Set changing parameters if in calibration mode:
c  = guess(1);
mu = guess(2);



% Turn bargaining shocks on 
bargain=1;

% Solve model with varying beta
solve_model

% Define the avg. 
wage_bar = 0.9340;

out = max(abs(wages(:) - wage_bar)); 

fprintf('Distance is %.4f for parameters (c,mu)=(%.4f, %.4f) \n', out, c, mu)

write = [guess' out];
dlmwrite('result.csv',write,'delimiter',',','-append');
