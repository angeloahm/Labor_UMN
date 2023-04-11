%% Simulation 
% This code creates a function to that takes the solution of the model (the vector theta) and simulate it for T periods
function [p_sim, theta_sim, f_sim, v_sim, u_sim] = simulate_model(T, theta)

% Read numerical/model parameters 
read_parameters

% Solve wage equation to find the vector of theta
%solve_model


%% Simulation preliminaries 

fprintf('=============================\n')
fprintf('Simulating the model...\n')
fprintf('=============================\n')


% Initial values for some variables
quarter = 0;
time    = 0;
u       = 0.05;

arrivaltime = exprnd(1/lambda, 1000);       %time between shocks 
totaltime   = cumsum(arrivaltime);          %total periods simulated (measured in quarters)
nshocks     = sum(totaltime(:) < T);        %total number of shocks simulated

shocks      = rand(1, nshocks);             

% Allocate simulated vectors
p_sim       = zeros(1,T);
theta_sim   = zeros(1,T);
f_sim       = zeros(1,T);
u_sim       = zeros(1,T);
v_sim       = zeros(1,T);

%% Run simulation
idx = find(y_grid==0);   %start simulation with y=0 


for n=1:nshocks

    % Time of the next Poisson shock
    time = time + arrivaltime(n);
    
    % Update unemployment according to the flow equation
    u = u + du(s, u, theta(idx));

    if shocks(n)<P(idx,idx-1)                   %prob of going down
        idx = idx-1;
    else                                   %prob of going up
        idx = idx+1;                    
    end


    if time+arrivaltime(n+1) >= quarter+1       %If next shocks is in another quarter...

        % Update quarter
        quarter = quarter+1;

        % Save variables for this quarter
        p_sim(quarter)      = p_grid(idx);  
        theta_sim(quarter)  = theta(idx);
        f_sim(quarter)      = f(theta(idx));
        v_sim(quarter)      = theta(idx)*u_sim(quarter);
        u_sim(quarter+1)    = u;

        if u>1
            keyboard
        end

        %End after 1212 periods
        if quarter>T
            break
        end

    end     %end quarter update

end         %end shocks


% Burn first 1,000 "quarters"
p_sim       = p_sim(1000:T);
theta_sim   = theta_sim(1000:T);
f_sim       = f_sim(1000:T);
v_sim       = v_sim(1000:T);
u_sim       = u_sim(1000:T);









