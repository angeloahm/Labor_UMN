% This code runs the simulation of the model using the effort policy and 
% the reservation wage obtained from 'solve_model.m'


shock = rand(N,T);                                             %generate random numbers
%index_initial_skill_distribution = randi(n_h, [N,1]);         %initial guess for the distribution of skills (index)
%index_initial_wages_distribution = randi(n_w, [N,1]);         %initial guess for the distribution of wages (index)


% We want to keep track the employment status, human capital, wage.
% We generate an initial skill distribution and employment status
index_simulated_skill = ones(N, T);             %I'm starting all agets with the smallest human capital level
index_simulated_wages = zeros(N, T);

simulated_employment = zeros(N, T);
simulated_effort = zeros(N,T);


% Resolving uncertainty 
fired = (shock<lambda);
dead = (shock<alpha);

for t=2:T       %for all periods
    for n=1:N   %for all agents
        
        if dead(n,t-1)==1     %if the agent dies, replace by the initial state
            simulated_employment(n,t)    = simulated_employment(n,1);
            index_simulated_skill(n,t)   = index_simulated_skill(n,1);
            index_simulated_wages(n,t)   = index_simulated_wages(n,1);
        else    %if alive
            if simulated_employment(n,t-1)==1  %if employed
                    simulated_employment(n,t) = 1-fired(n,t-1);     %0 if fired, 1 o/w
                    % human capital dynamics (fired vs. not fired)
                    index_simulated_skill(n,t) = fired(n,t-1)*max(index_simulated_skill(n,t-1)-psi_F, 1) + ...
                                                 (1-fired(n,t-1))*min(index_simulated_skill(n,t-1)+1,n_h);
                    index_simulated_wages(n,t) = (1-fired(n,t-1))*index_simulated_wages(n,t-1) + fired(n,t-1);
            else    %if unemployed
                    index_simulated_skill(n,t) = max(index_simulated_skill(n,t-1)-psi_U, 1);      %loose h continuously
                    simulated_effort(n,t) = policy_effort(index_simulated_skill(n,t));            %apply h in the effort policy 
                    match = (shock(n,t)<pi(simulated_effort(n,t)));                               %indicator of getting an offer
                    
                    % Draw wage from F
                    i_wage_drawn = sum(shock(n,t) >= F)+1;
                               
                    %match and get a good draw
                    simulated_employment(n,t) = match*(i_wage_drawn >= iw_res(index_simulated_skill(n,t)));                                                                                                            %at the policy_index
                    index_simulated_wages(n,t) = simulated_employment(n,t)*i_wage_drawn + (1-simulated_employment(n,t));
            end %end if employed  
        end %end if dead

    end %for all agents
end     %for all periods


% Exclude first 100 obs
index_simulated_skill = index_simulated_skill(:,101:T)  ;
index_simulated_wages = index_simulated_wages(:,101:T)  ;
simulated_employment  = simulated_employment(:,101:T)   ;
 
% We now use the last periods of the NxT matrix to compute the stationary
% distribution
dist = zeros(n_w, n_h);

for i_w=1:n_w
    for i_h=1:n_h
        logical = (index_simulated_wages==i_w) .* (index_simulated_skill==i_h);
        dist(i_w, i_h) = sum( logical(:) ) / (N*length(index_simulated_skill(1,:)));
    end
end




