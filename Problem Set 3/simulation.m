% Simulation 
% This code simulates the model for T periods and uses a HP filter to get
% the results from Table 3 in the paper. 

fprintf('=============================\n')
fprintf('Simulating the model...\n')
fprintf('=============================\n')

%% Simulation preliminaries 

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

% Use HP filter
[p_trend, p_cycle] = hpfilter(p_sim, Smoothing=10^5);
[theta_trend, theta_cycle] = hpfilter(theta_sim, Smoothing=10^5);
[f_trend, f_cycle] = hpfilter(f_sim, Smoothing=10^5);
[v_trend, v_cycle] = hpfilter(v_sim, Smoothing=10^5);
[u_trend, u_cycle] = hpfilter(u_sim, Smoothing=10^5);

%% Compute moments

%Standard deviation of de-trended series
sd_p        = std(log(p_sim(:) ./ p_trend(:)));
sd_f        = std(log(f_sim(:) ./ f_trend(:)));
sd_theta    = std(log(theta_sim(:) ./ theta_trend(:)));
sd_v        = std(log(v_sim(:) ./ v_trend(:)));
sd_u        = std(log(u_sim(:) ./ u_trend(:)));

%Autocorrelation
ac_p = autocorr(log(p_sim(:) ./ p_trend(:)), NumLags=2);
ac_f = autocorr(log(f_sim(:) ./ f_trend(:)), NumLags=2);
ac_theta = autocorr(log(theta_sim(:) ./ theta_trend(:)), NumLags=2);
ac_v = autocorr(log(v_sim(:) ./ v_trend(:)), NumLags=2);
ac_u = autocorr(log(u_sim(:) ./ u_trend(:)), NumLags=2);







