%% Simulation 
% This code simulates the model using eq. outcomes from 'solve_model.m'

fprintf('=============================\n')
fprintf('Simulating the model...\n')
fprintf('=============================\n')

% Allocate arrays 
iw_sim = ones(N, T);     %keep track of indices!
ib_sim = ones(N, T);     %keep track of indices!
s_sim = ones(N, T);
u_sim = ones(N, T);    
%u_sim(1:N/2,T) = 1;

shock = rand(N,T);                      %generate random numbers

mc = dtmc(P);                           %create a markov chain from P
iz_sim = simulate(mc, T);                %simulate markov chain (indexes)

for t=1:T-1
    for n=1:N
        
        s_sim(n,t) = policy_effort(iw_sim(n,t), iz_sim(t));
        
        %Unemployed agent at t
        if u_sim(n,t) == 1
            if shock(n,t)<lambda(s_sim(n,t))*p(theta(iw_sim(n,t), iz_sim(t))) %matches a firm and gets a job
                u_sim(n,t+1) = 0;                                            %becomes employed
                iw_sim(n,t+1) = idx_wage(iw_sim(n,t), iz_sim(t));             %wage from tomorrow on
                ib_sim(n,t+1) = ib_sim(n,t+1);                                           
            else                                                            %does not match a firm
                u_sim(n,t+1) = 1;   %remains unemployed
                iw_sim(n,t+1) = iw_sim(n,t+1);
                ib_sim(n,t+1) = (shock(n,t)<0.1)*1 + (shock(n,t)>0.1)*ib_sim(n,t);    %keep unemployment benefit w.p 0.9
            end     %end matching case
        %Employed agent at t
        else
            s_sim(n,t) = policy_effort(iw_sim(n,t), iz_sim(t));
            if shock(n,t)<delta(iw_sim(n,t), iz_sim(t+1))                     %hit by separation shock
                u_sim(n,t+1) = 1;   %becomes unemployed
                iw_sim(n,t+1) = iw_sim(n,t);
                ib_sim(n,t+1) = (iw_sim(n,t)>1)*floor(iw_sim(n,t)/2) + (iw_sim(n,t)==1);    %first UI is w/2
            else
                u_sim(n,t+1) = 0;           %remains employed
                iw_sim(n,t+1) = iw_sim(n,t);  %keep same wage
                ib_sim(n,t+1) = ib_sim(n,t);
            end     %end separation shock case
        end         %end employment status

        


    end     %end loop N
end         %end loop T


iw_sim = iw_sim(:, burn:T);
iz_sim = iz_sim(burn:T);
u_sim = u_sim(:, burn:T);
ib_sim = ib_sim(:, burn:T);

% Print statistics
z_sim=z_grid(iz_sim(:));
fprintf('Mean unemployment rate:                  %.4f \n', mean(u_sim(:)));
fprintf('SD unemployment rate/SD productivity:    %.4f \n', std(u_sim(:))/std(z_sim(:)) );


