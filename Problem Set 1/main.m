% Code to solve a simplified version of Ljungqvist and Sargent (1998)
% Angelo Mendes 
% University of Minnesota

clc
clear
close all

%Start timer
tic


read_parameters

if eq==1
    %Find equilibrium
    fprintf('=========== FINDING EQ. TAXES ===========\n')
    tau_star = bisection(@budget_balance, 0.015, 0.02, 100, 0.001);
    %Solve model again with the equilibrium tau
    tau=tau_star;
    solve_model;
    simulation;

else
    %Solve model for a specific tau
    tau = tau_guess;
    solve_model;
    simulation;
    
end


figure(1)
plot(h_grid, policy_effort)
title('Policy Effort')
xlabel('Human capital')
ylabel('s(h)')
saveas(gcf,'effort_policy.png')

figure(2)
plot(h_grid, w_res)
title('Reservation wage')
xlabel('Human capital')
ylabel('Reservation wage')
saveas(gcf,'res_wage.png')

figure(3)
contourf(w_grid,h_grid,policy_employed',[0 1])
xlabel('w')
ylabel('h')
title('Employment')
colorbar
saveas(gcf,'employment_heatmap.png')



% tau_grid = linspace(0,1,10);
% vector_excess = ones(1, 10);
% 
% for t=1:10
%     tau = tau_grid(t);
%     out = budget_balance(tau);
%     vector_excess(t) = out;
% end
% 
% figure(1)
% plot(tau_grid, vector_excess,tau_grid, zeros(1,length(tau_grid)))



% End timer
toc











