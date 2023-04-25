%% Problem Set 6 - 8581
% Author: Angelo Avelar Hermeto Mendes
% This code estimates the model in Problem Set 5. 
% We look for a value of unemployment benefit b that minimizes the distance
% of the average unemployment and the target 6.5%

lb = 0.3;     %lower bound of values for b
ub = 1;   %upper bound of values for b

myoptions = optimoptions('patternsearch','MaxTime',3600*8);
[solution, fval] = patternsearch(@objective,0.75,[],[],[],[],lb,ub,[],myoptions);


objective(solution);

% Evaluate the derivative
g0 = g(0.83);
g1 = g(0.83 + 0.001);

diff = (g1 - g0) / 0.001;
disp(diff);



