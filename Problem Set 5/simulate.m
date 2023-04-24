%% Simulation 
% This code simulates the model using eq. outcomes from 'solve_model.m'

fprintf('=============================\n')
fprintf('Simulating the model...\n')
fprintf('=============================\n')

rng(17)

% Allocate arrays 
iw_sim = zeros(N, T);       %keep track of indices! (all agents start unemployed: iw=0)
ib_sim = ones(N, T);        %keep track of indices!
s_sim = zeros(N, T);
e_sim = zeros(N, T);    

%generate random numbers
shock_l = rand(N,T);
shock_ui = rand(N,T);
shock_b = rand(N,1);
shock_p = rand(N,T);
shock_delta = rand(N,T);
shock_z = rand(1,T);

% iz_sim(1) = 8;
% for t=2:T
%     [~, iz_sim(t)] = max(shock_z(t)<=cumsum(P(iz_sim(t-1),:)));
% end

mc = dtmc(P);                           %create a markov chain from P
iz_sim = simulate(mc, T);               %simulate markov chain (indexes)

%Initial states for the benefit are uniformly distributed
b_pdf=ones(n_w,1)*(1/n_w);
b_cdf=cumsum(b_pdf);

for n=1:N

    [~,ib_sim(n,1)]=max(shock_b(n,1)<=b_cdf);               %set initial stage for benefits
    s_sim(n,1) = policy_effort(ib_sim(n,1), iz_sim(1));

    for t=2:T
        
        %Unemployed agent at t
        if e_sim(n,t-1) == 0
           
            ib_sim(n,t) = (shock_ui(n,t)<=0.1)*1 + (shock_ui(n,t)>0.1)*ib_sim(n,t-1);      %keep unemployment benefit w.p 0.9
            %matches a firm
            if shock_l(n,t)<=lambda(s_sim(n,t-1))                                   
               %Note that the probability of a meeting depends of the
               %effort made according to (b,z) and just at t+1 I get to
               %find if I meet or not the firm...
                
                %gets a job
                if shock_p(n,t)<=p( theta( idx_wage(ib_sim(n,t), iz_sim(t)), iz_sim(t) ) )            
                    e_sim(n,t)  = 1;                                                   %becomes employed
                    iw_sim(n,t) = idx_wage(ib_sim(n,t), iz_sim(t));                    %wage from tomorrow on
                    ib_sim(n,t) = 0;
                    s_sim(n,t) = 0;                                                    %will not search since found a job

                %does not get a job
                else                                                                
                    e_sim(n,t) = 0;                                                 %remains unemployed
                    iw_sim(n,t) = 0;
                    s_sim(n,t) = policy_effort(ib_sim(n,t), iz_sim(t));             %search according to policy
                end %end p shock

            %does not match a firm
            else                                                                  %does not match a firm
                e_sim(n,t) = 0;                                                   %remains unemployed
                iw_sim(n,t) = 0;
                s_sim(n,t) = policy_effort(ib_sim(n,t), iz_sim(t));               %search policy today
            end     %end lambda shock

        %Employed agent at t
        else
            
            %hit by separation shock
            if shock_delta(n,t)<=delta(w_grid(iw_sim(n,t-1)), z_grid(iz_sim(t)))
                e_sim(n,t) = 0;                                                                     %becomes unemployed
                ib_sim(n,t) = iw_sim(n,t-1);%(iw_sim(n,t-1)>1)*floor(iw_sim(n,t-1)/2) + (iw_sim(n,t-1)==1);%benefit indexed according to wage today!
                iw_sim(n,t) = 0;               
                s_sim(n,t) = policy_effort(ib_sim(n,t), iz_sim(t));                                   %employed agent does not exhert any effort

            %does not hit by separation shock
            else
                e_sim(n,t) = 1;               %remains employed 
                iw_sim(n,t) = iw_sim(n,t-1);  %keep same wage
                ib_sim(n,t) = 0;
                s_sim(n,t) = 0;             %employed agent does not exhert any effort
            end     %end separation shock case
        end         %end employment status


    end     %end loop T
end         %end loop N


iw_sim = iw_sim(:, burn:T);
iz_sim = iz_sim(burn:T);
e_sim = e_sim(:, burn:T);
ib_sim = ib_sim(:, burn:T);

% Print statistics
z_sim=z_grid(iz_sim(:));
fprintf('Mean unemployment rate:                  %.4f \n', 1-mean(mean(e_sim)));
fprintf('SD unemployment rate/SD productivity:    %.4f \n', std((1-e_sim(:)))/std(z_sim(:)) );


