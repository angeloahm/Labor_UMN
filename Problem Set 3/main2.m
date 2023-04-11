%% Problem Set 3 (8581) - Shimer (2005)
% Angelo Mendes
% University of Minnesota

clear all 

%% Table 3
%Read parameters 
read_parameters

bargain=0;  %no bargain shocks
leisure=0;  %no leisure shocks

%Solve the model to find theta
solve_model

% Simulate the model 
sd_p_vec        = zeros(1,W);
sd_theta_vec    = zeros(1,W);
sd_f_vec        = zeros(1,W);
sd_v_vec        = zeros(1,W);
sd_u_vec        = zeros(1,W);
ac_p_vec        = zeros(1,W);
ac_theta_vec    = zeros(1,W);
ac_f_vec        = zeros(1,W);
ac_v_vec        = zeros(1,W);
ac_u_vec        = zeros(1,W);
correl          = zeros(5,5,W);
elasticity_vec  = zeros(1,W);

for w=1:W
    
    [p_sim, theta_sim, f_sim, v_sim, u_sim] = simulate_model(1212, theta);

    %Compute wages 
    wages_sim = (1-beta)*z + beta*(p_sim + c*theta_sim );
    coef_reg = cov(log(wages_sim), log(p_sim)) / var(p_sim);
    elasticity = coef_reg(1,2);
    elasticity_vec(w) = elasticity;
    
    % Use HP filter
    [p_trend, p_cycle] = hpfilter(p_sim, Smoothing=10^5);
    [theta_trend, theta_cycle] = hpfilter(theta_sim, Smoothing=10^5);
    [f_trend, f_cycle] = hpfilter(f_sim, Smoothing=10^5);
    [v_trend, v_cycle] = hpfilter(v_sim, Smoothing=10^5);
    [u_trend, u_cycle] = hpfilter(u_sim, Smoothing=10^5);

    % Compute moments
    % Standard deviation of de-trended series
    sd_p        = std(log(p_sim(:) ./ p_trend(:)));
    sd_theta    = std(log(theta_sim(:) ./ theta_trend(:)));
    sd_f        = std(log(f_sim(:) ./ f_trend(:)));
    sd_v        = std(log(v_sim(:) ./ v_trend(:)));
    sd_u        = std(log(u_sim(:) ./ u_trend(:)));

    % Autocorrelation
    ac_p = autocorr(log(p_sim(:) ./ p_trend(:)), NumLags=10);
    ac_f = autocorr(log(f_sim(:) ./ f_trend(:)), NumLags=10);
    ac_theta = autocorr(log(theta_sim(:) ./ theta_trend(:)), NumLags=10);
    ac_v = autocorr(log(v_sim(:) ./ v_trend(:)), NumLags=10);
    ac_u = autocorr(log(u_sim(:) ./ u_trend(:)), NumLags=10);

    %Compute correlation matrix
    matrix = [log(u_sim(:) ./ u_trend(:))'; log(v_sim(:) ./ v_trend(:))'; log(theta_sim(:) ./ theta_trend(:))'; ...
              log(f_sim(:) ./ f_trend(:))'; log(p_sim(:) ./ p_trend(:))'];
    
    corr=corrcoef(matrix');

    sd_p_vec(w) = sd_p;
    sd_theta_vec(w) = sd_theta;
    sd_f_vec(w)     = sd_f;
    sd_v_vec(w)     = sd_v;
    sd_u_vec(w)     = sd_u;
    ac_p_vec(w)     = ac_p(2);
    ac_theta_vec(w) = ac_theta(2);
    ac_f_vec(w)     = ac_f(2);
    ac_v_vec(w)     = ac_v(2);
    ac_u_vec(w)     = ac_u(2);
    
    for i=1:5
        for j=1:5
            correl(i,j,w) = corr(i,j);
        end
    end


end


%Print results
fprintf('        Table 3 - Labor Productivity Shocks         \n')
disp([mean(sd_u_vec(:)), mean(sd_v_vec(:)), mean(sd_theta_vec(:)), mean(sd_f_vec(:)), mean(sd_p_vec(:))])
disp([std(sd_u_vec(:)), std(sd_v_vec(:)), std(sd_theta_vec(:)), std(sd_f_vec(:)), std(sd_p_vec(:))])
disp([mean(ac_u_vec(:)) mean(ac_v_vec(:)) mean(ac_theta_vec(:)) mean(ac_f_vec(:)) mean(ac_p_vec(:))])
disp([std(ac_u_vec(:)), std(ac_v_vec(:)), std(ac_theta_vec(:)), std(ac_f_vec(:)), std(ac_p_vec(:))])

fprintf('        Correlation Matrix - Labor Shocks         \n')
disp([mean(correl(1,1,:)), mean(correl(1,2,:)), mean(correl(1,3,:)), mean(correl(1,4,:)), mean(correl(1,5,:))])
disp([std(correl(1,1,:)), std(correl(1,2,:)), std(correl(1,3,:)), std(correl(1,4,:)), std(correl(1,5,:))])
disp([mean(correl(2,1,:)), mean(correl(2,2,:)), mean(correl(2,3,:)), mean(correl(2,4,:)), mean(correl(2,5,:))])
disp([std(correl(2,1,:)), std(correl(2,2,:)), std(correl(2,3,:)), std(correl(2,4,:)), std(correl(2,5,:))])
disp([mean(correl(3,1,:)), mean(correl(3,2,:)), mean(correl(3,3,:)), mean(correl(3,4,:)), mean(correl(3,5,:))])
disp([std(correl(3,1,:)), std(correl(3,2,:)), std(correl(3,3,:)), std(correl(3,4,:)), std(correl(3,5,:))])
disp([mean(correl(4,1,:)), mean(correl(4,2,:)), mean(correl(4,3,:)), mean(correl(4,4,:)), mean(correl(4,5,:))])
disp([std(correl(4,1,:)), std(correl(4,2,:)), std(correl(4,3,:)), std(correl(4,4,:)), std(correl(4,5,:))])
disp([mean(correl(5,1,:)), mean(correl(5,2,:)), mean(correl(5,3,:)), mean(correl(5,4,:)), mean(correl(5,5,:))])
disp([std(correl(5,1,:)), std(correl(5,2,:)), std(correl(5,3,:)), std(correl(5,4,:)), std(correl(5,5,:))])


%% Some plots
% Plot distribution over elasticity
[dist, supp] = ksdensity(elasticity_vec(:));
figure(1)
plot(supp, dist, LineWidth=2)
title('Elasticity of wages wrt productivity')
xlabel('Elasticity')
saveas(gcf, 'elasticity.png')

% Plot Beveridge Curve
figure(2)
scatter(u_sim, v_sim)
title('Beveridge Curve')
xlabel('u')
ylabel('v')
saveas(gcf, 'beveridge.png')


%% Bargaining shocks

%Find new values for (c,mu)
guess = [c; mu];

%out = objective_calibration(guess);

[parameters, out] = fminsearch(@objective_calibration,myguess);


%% Leisure shocks
bargain=0;
leisure=1;

read_parameters;
solve_model;

% Simulate the model 
sd_p_vec        = zeros(1,W);
sd_theta_vec    = zeros(1,W);
sd_f_vec        = zeros(1,W);
sd_v_vec        = zeros(1,W);
sd_u_vec        = zeros(1,W);
ac_p_vec        = zeros(1,W);
ac_theta_vec    = zeros(1,W);
ac_f_vec        = zeros(1,W);
ac_v_vec        = zeros(1,W);
ac_u_vec        = zeros(1,W);
correl          = zeros(5,5,W);
elasticity_vec  = zeros(1,W);

for w=1:W
    
    [p_sim, theta_sim, f_sim, v_sim, u_sim] = simulate_model(1212, theta);

    %Compute wages 
    wages_sim = (1-beta)*(p_sim-0.0001) + beta*(p_sim + c*theta_sim );
    coef_reg = cov(log(wages_sim), log(p_sim)) / var(p_sim);
    elasticity = coef_reg(1,2);
    elasticity_vec(w) = elasticity;
    
    % Use HP filter
    [p_trend, p_cycle] = hpfilter(p_sim, Smoothing=10^5);
    [theta_trend, theta_cycle] = hpfilter(theta_sim, Smoothing=10^5);
    [f_trend, f_cycle] = hpfilter(f_sim, Smoothing=10^5);
    [v_trend, v_cycle] = hpfilter(v_sim, Smoothing=10^5);
    [u_trend, u_cycle] = hpfilter(u_sim, Smoothing=10^5);

    % Compute moments
    % Standard deviation of de-trended series
    sd_p        = std(log(p_sim(:) ./ p_trend(:)));
    sd_theta    = std(log(theta_sim(:) ./ theta_trend(:)));
    sd_f        = std(log(f_sim(:) ./ f_trend(:)));
    sd_v        = std(log(v_sim(:) ./ v_trend(:)));
    sd_u        = std(log(u_sim(:) ./ u_trend(:)));

    % Autocorrelation
    ac_p = autocorr(log(p_sim(:) ./ p_trend(:)), NumLags=10);
    ac_f = autocorr(log(f_sim(:) ./ f_trend(:)), NumLags=10);
    ac_theta = autocorr(log(theta_sim(:) ./ theta_trend(:)), NumLags=10);
    ac_v = autocorr(log(v_sim(:) ./ v_trend(:)), NumLags=10);
    ac_u = autocorr(log(u_sim(:) ./ u_trend(:)), NumLags=10);

    %Compute correlation matrix
    matrix = [log(u_sim(:) ./ u_trend(:))'; log(v_sim(:) ./ v_trend(:))'; log(theta_sim(:) ./ theta_trend(:))'; ...
              log(f_sim(:) ./ f_trend(:))'; log(p_sim(:) ./ p_trend(:))'];
    
    corr=corrcoef(matrix');

    sd_p_vec(w) = sd_p;
    sd_theta_vec(w) = sd_theta;
    sd_f_vec(w)     = sd_f;
    sd_v_vec(w)     = sd_v;
    sd_u_vec(w)     = sd_u;
    ac_p_vec(w)     = ac_p(2);
    ac_theta_vec(w) = ac_theta(2);
    ac_f_vec(w)     = ac_f(2);
    ac_v_vec(w)     = ac_v(2);
    ac_u_vec(w)     = ac_u(2);
    
    for i=1:5
        for j=1:5
            correl(i,j,w) = corr(i,j);
        end
    end


end


%Print results
fprintf('        Leisure Shocks         \n')
disp([mean(sd_u_vec(:)), mean(sd_v_vec(:)), mean(sd_theta_vec(:)), mean(sd_f_vec(:)), mean(sd_p_vec(:))])
disp([std(sd_u_vec(:)), std(sd_v_vec(:)), std(sd_theta_vec(:)), std(sd_f_vec(:)), std(sd_p_vec(:))])
disp([mean(ac_u_vec(:)) mean(ac_v_vec(:)) mean(ac_theta_vec(:)) mean(ac_f_vec(:)) mean(ac_p_vec(:))])
disp([std(ac_u_vec(:)), std(ac_v_vec(:)), std(ac_theta_vec(:)), std(ac_f_vec(:)), std(ac_p_vec(:))])

% Plot Beveridge Curve
figure(3)
scatter(u_sim, v_sim)
title('Beveridge Curve')
xlabel('u')
ylabel('v')
saveas(gcf, 'beveridge_leisure.png')

% Plot distribution over elasticity
[dist, supp] = ksdensity(elasticity_vec(:));
figure(4)
plot(supp, dist, LineWidth=2)
title('Elasticity of wages wrt productivity')
xlabel('Elasticity')
saveas(gcf, 'elasticity_leisure.png')

fprintf('        Correlation Matrix - Leisure Shocks         \n')
disp([mean(correl(1,1,:)), mean(correl(1,2,:)), mean(correl(1,3,:)), mean(correl(1,4,:)), mean(correl(1,5,:))])
disp([std(correl(1,1,:)), std(correl(1,2,:)), std(correl(1,3,:)), std(correl(1,4,:)), std(correl(1,5,:))])
disp([mean(correl(2,1,:)), mean(correl(2,2,:)), mean(correl(2,3,:)), mean(correl(2,4,:)), mean(correl(2,5,:))])
disp([std(correl(2,1,:)), std(correl(2,2,:)), std(correl(2,3,:)), std(correl(2,4,:)), std(correl(2,5,:))])
disp([mean(correl(3,1,:)), mean(correl(3,2,:)), mean(correl(3,3,:)), mean(correl(3,4,:)), mean(correl(3,5,:))])
disp([std(correl(3,1,:)), std(correl(3,2,:)), std(correl(3,3,:)), std(correl(3,4,:)), std(correl(3,5,:))])
disp([mean(correl(4,1,:)), mean(correl(4,2,:)), mean(correl(4,3,:)), mean(correl(4,4,:)), mean(correl(4,5,:))])
disp([std(correl(4,1,:)), std(correl(4,2,:)), std(correl(4,3,:)), std(correl(4,4,:)), std(correl(4,5,:))])
disp([mean(correl(5,1,:)), mean(correl(5,2,:)), mean(correl(5,3,:)), mean(correl(5,4,:)), mean(correl(5,5,:))])
disp([std(correl(5,1,:)), std(correl(5,2,:)), std(correl(5,3,:)), std(correl(5,4,:)), std(correl(5,5,:))])