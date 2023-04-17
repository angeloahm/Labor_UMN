
%% Problem Set 5
% Angelo Avelar Hermeto Mendes

clear 
 

% Read model/numeric parameters
read_parameters

% Solve the model 
tic
solve_model
toc

% Plots
[W,Z] = meshgrid(z_grid, w_grid);
figure(1)
surf(W, Z, policy_effort)
xlabel('Aggregate shock')
ylabel('Benefit')
zlabel('Search Policy')
saveas(gcf, 'search.png')

figure(2)
surf(W, Z, policy_wage)
xlabel('Aggregate shock')
ylabel('Benefit')
zlabel('Wage Posting Policy')
saveas(gcf, 'wage.png')


figure(3)
plot(w_grid, p(theta(:,7)), LineWidth=2, Color='black')
hold on 
plot(w_grid, p(theta(:,15)), LineStyle='--', LineWidth=2, Color='red')
legend(['z=1'; 'z_H'])
title('Job Finding Rate (p(\theta))')
xlabel('Posted Wage')
ylabel('Job Finding Rate')
saveas(gcf, 'p.png')

figure(4)
plot(w_grid, theta(:,7), LineWidth=2, Color='black')
hold on 
plot(w_grid, theta(:,15), LineStyle='--', LineWidth=2, Color='red')
legend(['z=1'; 'z_H'])
title('Market Tightness (\theta)')
xlabel('Posted Wage')
ylabel('Market Tightness (\theta)')
saveas(gcf, 'theta.png')


% Run simulation and compute moments
tic
simulate
toc












